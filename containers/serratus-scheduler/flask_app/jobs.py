import csv
from itertools import islice
from datetime import datetime

from flask import Blueprint, jsonify, request, abort, render_template
from prometheus_client import Counter

from . import db
from .metrics import REGISTRY

bp = Blueprint("jobs", __name__, url_prefix="/jobs")

done_items_counter = Counter(
    "scheduler_done_items", "Items successfully processed", ["kind"], registry=REGISTRY
)


@bp.route("/add_sra_run_info/<filename>", methods=["POST"])
def add_sra_runinfo(filename):
    """Read the CSV into the DB.

    TODO: Use the filename part of the URL to prevent the same file from
    being uploaded multiple times.
    """
    insert_count = 0
    session = db.get_session()

    csv_data = csv.DictReader(map(lambda b: b.decode(), request.stream))
    while True:
        # Using bulk_insert_mappings, using the recipe recommended here:
        # https://docs.sqlalchemy.org/en/13/faq/performance.html
        mappings = [
            dict(state="new", sra_run_info=line) for line in islice(csv_data, 10000)
        ]

        if not mappings:
            break

        session.bulk_insert_mappings(db.Accession, mappings)
        insert_count += len(mappings)
        print("Reading data... {}".format(insert_count))

    session.commit()

    total_count = session.query(db.Accession).count()
    return jsonify({"inserted_rows": insert_count, "total_rows": total_count,})


@bp.route("/", methods=["GET"])
def show_jobs():
    session = db.get_session()
    accs = session.query(db.Accession).all()
    blocks = (
        session.query(db.Block, db.Accession)
        .filter(db.Block.acc_id == db.Accession.acc_id)
        .all()
    )
    return render_template("job_list.html", accs=accs, blocks=blocks)


@bp.route("/split/", methods=["POST"])
@bp.route("/dl/", methods=["POST"])
def start_split_job():
    """Get a job id and mark it as "running" in the DB
    returns: a job JSON file"""
    session = db.get_session()

    # get an item where state = new and update its state
    acc = (
        session.query(db.Accession)
        .filter_by(state="new")
        .with_for_update(skip_locked=True)
        .first()
    )

    if acc is None:
        return jsonify({"action": "wait"})

    acc.state = "splitting"
    acc.split_start_time = datetime.now()
    acc.split_end_time = None
    acc.split_worker = request.args.get("worker_id")

    session.add(acc)
    response = acc.to_dict()
    response["id"] = acc.acc_id
    session.commit()

    response["action"] = "process"
    response["split_args"] = db.get_config_val("DL_ARGS")

    # Send the response as JSON
    return jsonify(response)


@bp.route("/split/<acc_id>", methods=["POST"])
@bp.route("/dl/<acc_id>", methods=["POST"])
def finish_split_job(acc_id):
    """Mark a split job as finished, setting off N alignment jobs.

    param acc_id: The acc_id, as returned by start_split_job
    param state: Resulting state, see sched_states.dia
    param N_paired: Number of paired blocks the split resulted in
    param N_unpaired: Number of unpaired blocks the split resulted in

    Examples:

    Successful split into 10 blocks:
        http://scheduler/jobs/split/1?state=split_done&N_paired=10

    Failed split (bug in code or bad data):
        http://scheduler/jobs/split/1?state=split_err

    Terminated node (please redo this work):
        http://scheduler/jobs/split/1?state=new
    """
    state = request.args.get("state")
    try:
        n_unpaired = int(request.args.get("N_unpaired", 0))
        n_paired = int(request.args.get("N_paired", 0))
    except TypeError:
        abort(400)

    if state == "terminated":
        state = "new"

    if state == "done":
        state = "split_done"

    # Update the accessions table
    if state not in ("new", "split_err", "split_done"):
        raise ValueError("Invalid State {}".format(state))

    session = db.get_session()
    acc = (
        session.query(db.Accession)
        .filter_by(acc_id=int(acc_id))
        .with_for_update()
        .one()
    )

    if acc.state != "splitting":
        abort(400)
    acc.state = state
    acc.split_end_time = datetime.now()

    if state in ("new", "split_err"):
        # Not ready to process more
        session.commit()
        return jsonify({"result": "success"})

    acc.contains_paired = n_paired > 0
    acc.contains_unpaired = n_unpaired > 0

    # AB: if there is paired-read data then ignore unpaired data
    blocks = n_paired or n_unpaired

    # Insert N align jobs into the alignment table.
    for i in range(int(blocks)):
        block = db.Block(state="new", acc_id=acc.acc_id, n=i)
        session.add(block)

    acc.blocks = blocks
    session.add(acc)
    session.commit()

    done_items_counter.labels("dl").inc()

    return jsonify({"result": "success", "inserted_rows": blocks,})


@bp.route("/align/", methods=["POST"])
def start_align_job():
    """Get a job id and mark it as "running" in the DB

    returns: a job JSON file"""
    session = db.get_session()
    # get an item where state = new and update its state
    block = (
        session.query(db.Block)
        .filter(db.Block.state == "new")
        .with_for_update(skip_locked=True)
        .first()
    )

    if block is None:
        # TODO: If we got no results, then could be:
        #    * Waiting for split nodes: hang tight
        #    * No more work: shutdown
        return jsonify({"action": "wait"})

    block.state = "aligning"
    block.align_start_time = datetime.now()
    block.align_end_time = None
    block.align_worker = request.args.get("worker_id")
    session.add(block)
    session.commit()

    response = block.to_dict()

    # Query acc *without* for_update, since we're just using it to get extra
    # information about the block.
    acc = session.query(db.Accession).filter_by(acc_id=block.acc_id).one()

    response.update(acc.to_dict())
    response["id"] = block.block_id
    response["align_args"] = db.get_config_val("ALIGN_ARGS")
    response["genome"] = db.get_config_val("GENOME")
    response["action"] = "process"

    # Send the response as JSON
    return jsonify(response)


@bp.route("/align/<block_id>", methods=["POST"])
def finish_align_job(block_id):
    """Finished job, block_id is the same parameter from the start_job,
    state is one of (new, done, fail)"""
    state = request.args.get("state")

    if state == "terminated":
        state = "new"

    if state not in db.BLOCK_STATES:
        abort(400)

    session = db.get_session()
    block = session.query(db.Block).filter_by(block_id=block_id).with_for_update().one()
    block.state = state
    block.align_end_time = datetime.now()
    session.commit()

    done_items_counter.labels("align").inc()

    return jsonify({"result": "success"})


@bp.route("/merge/", methods=["POST"])
def start_merge_job():
    session = db.get_session()
    # Exclude accs with blocks in "new", "fail", or in progress state.
    exclude_accs = (
        session.query(db.Block.acc_id)
        .distinct()
        .filter(db.Block.state.in_(("new", "fail", "aligning")))
        .subquery()
    )

    acc = (
        session.query(db.Accession)
        .filter_by(state="split_done")
        .filter(~(db.Accession.acc_id.in_(exclude_accs)))
        .with_for_update(skip_locked=True)
        .first()
    )

    if acc is None:
        # TODO: Think about this behaviour
        return jsonify({"action": "wait"})

    acc.state = "merging"
    acc.merge_start_time = datetime.now()
    acc.merge_end_time = None
    acc.merge_worker = request.args.get("worker_id")
    session.add(acc)
    session.commit()

    response = acc.to_dict()
    response["id"] = acc.acc_id
    response["action"] = "process"
    response["genome"] = db.get_config_val("GENOME")

    # Send the response as JSON
    return jsonify(response)


@bp.route("/merge/<acc_id>", methods=["POST"])
def finish_merge_job(acc_id):
    state = request.args.get("state")

    if state == "terminated":
        state = "split_done"

    if state == "done":
        state = "merge_done"

    if state not in ("split_done", "merge_err", "merge_done"):
        abort(400)

    session = db.get_session()
    acc = (
        session.query(db.Accession)
        .filter_by(acc_id=int(acc_id))
        .with_for_update()
        .one()
    )

    if acc.state != "merging":
        abort(400)

    acc.state = state
    acc.merge_end_time = datetime.now()
    session.commit()

    done_items_counter.labels("merge").inc()

    return jsonify({"result": "success"})
