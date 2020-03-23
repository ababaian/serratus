from flask import Blueprint, jsonify, request, abort, render_template
from sqlalchemy import func
from . import db

bp = Blueprint('jobs', __name__, url_prefix='/jobs')

@bp.route('/reset_db', methods=['POST'])
def reset_db():
    db.init_db()
    return jsonify('Database cleared')

@bp.route('/add_sra_run_info/<filename>', methods=['POST'])
def add_sra_runinfo(filename):
    ## Read the CSV into the DB
    import csv, io
    insert_count = 0
    session = db.get_session()

    csv_data = io.StringIO(request.data.decode())
    for line in csv.DictReader(csv_data):
        insert_count += 1
        acc = db.Accession(state='new', sra_run_info=line)
        session.add(acc)

    total_count = session.query(db.Accession).count()
    session.commit()

    return jsonify({
        'inserted_rows': insert_count,
        'total_rows': total_count,
    })

@bp.route('/split/', methods=['POST'])
def start_split_job():
    """Get a job id and mark it as "running" in the DB
    returns: a job JSON file"""
    session = db.get_session()

    # get an item where state = new and update its state
    acc = session.query(db.Accession)\
        .filter_by(state='new')\
        .first()

    if acc is None:
        # TODO: Think about this behaviour
        return jsonify({'action': 'wait'})

    acc.state = 'splitting'
    session.add(acc)
    response = acc.to_dict()
    session.commit()
    response['action'] = 'process'

    response['split_args'] = ""

    # Send the response as JSON
    return jsonify(response)

@bp.route('/', methods=['GET'])
def show_jobs():
    session = db.get_session()
    accs = session.query(db.Accession).all()
    blocks = session.query(db.Block, db.Accession)\
        .filter(db.Block.acc_id == db.Accession.acc_id)\
        .all()
    return render_template('job_list.html', accs=accs, blocks=blocks)


@bp.route('/split/<acc_id>', methods=['POST'])
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
    state = request.args.get('state')
    try:
        n_unpaired = int(request.args.get('N_unpaired', 0))
        n_paired = int(request.args.get('N_paired', 0))
    except TypeError:
        abort(400)

    if state == "terminated":
        state = "new"

    # Update the accessions table
    if state not in ('new', 'split_err', 'split_done'):
        abort(400)

    session = db.get_session()
    acc = session.query(db.Accession)\
        .filter_by(acc_id=int(acc_id))\
        .one()

    if acc.state != 'splitting':
        abort(400)
    acc.state = state

    if state in ('new', 'split_err'):
        # Not ready to process more
        session.commit()
        return jsonify({'result': 'success'})

    acc.contains_paired = n_paired > 0
    acc.contains_unpaired = n_unpaired > 0

    # AB: if there is paired-read data then ignore unpaired data
    number_split = n_paired or n_unpaired

    # Insert N align jobs into the alignment table.
    for i in range(int(number_split)):
        block = db.Block(state='new', acc_id=acc.acc_id, n=i)
        session.add(block)

    row_count = session.query(db.Block).count()
    session.commit()

    return jsonify({
        'result': 'success',
        'inserted_rows': number_split,
        'total_rows': row_count,
    })

@bp.route('/align/', methods=['POST'])
def start_align_job():
    """Get a job id and mark it as "running" in the DB

    returns: a job JSON file"""
    session = db.get_session()
    # get an item where state = new and update its state
    query = session.query(db.Block, db.Accession)\
        .filter(db.Block.state == 'new')\
        .filter(db.Block.acc_id == db.Accession.acc_id)\
        .first()

    if query is None:
        # TODO: If we got no results, then could be:
        #    * Waiting for split nodes: hang tight
        #    * No more work: shutdown
        return jsonify({'action': 'wait'})

    block, acc = query

    block.state = 'aligning'
    session.add(block)

    response = block.to_dict()
    response.update(acc.to_dict())
    session.commit()

    # TODO Move these into the database
    response['align_args'] = "--very-sensitive-local"
    response['genome'] = "hgr1"
    response['action'] = "process"

    # Send the response as JSON
    return jsonify(response)

@bp.route('/align/<block_id>', methods=['POST'])
def finish_align_job(block_id):
    """Finished job, block_id is the same parameter from the start_job,
    state is one of (new, done, fail)"""
    state = request.args.get('state')

    if state == "terminated":
        state = "new"

    if state not in db.BLOCK_STATES:
        abort(400)

    session = db.get_session()
    block = session.query(db.Block).filter_by(block_id=block_id).one()
    block.state = state

    # Check the other blocks---are there any waiting or running still?
    state_counts_q = session.query(db.Block.state,
                                   func.count(db.Block.block_id))\
        .filter(db.Block.acc_id == block.acc_id)\
        .group_by(db.Block.state)\
        .all()

    state_counts = {k: 0 for k in db.BLOCK_STATES}
    state_counts.update(state_counts_q)

    if state_counts['new'] > 0 or state_counts['aligning'] > 0:
        # This is not the last block
        session.commit()
        return jsonify({
            'result': 'success'
        })
    elif state_counts['fail'] > 0:
        # :(  TODO What to do here?  Ideally we move Accession into a failed
        # state when a single align fails, so this should be no-op?
        session.commit()
        return jsonify({
            'result': ':('
        })
    else:
        # All blocks are done.  Move Accession into the merge_wait state.
        accession = session.query(db.Accession)\
            .filter(db.Accession.acc_id == block.acc_id)\
            .one()
        accession.state = 'merge_wait'
        session.commit()
        return jsonify({
            'result': 'success',
        })

@bp.route('/merge/', methods=['POST'])
def start_merge_job():
    session = db.get_session()
    # Get an Acc where all blocks have been aligned.
    acc = session.query(db.Accession)\
        .filter_by(state='merge_wait')\
        .first()

    if acc is None:
        # TODO: Think about this behaviour
        return jsonify({'action': 'wait'})

    acc.state = 'merging'
    session.add(acc)
    response = acc.to_dict()
    session.commit()
    response['action'] = 'process'

    response['merge_args'] = ""

    # Send the response as JSON
    return jsonify(response)

@bp.route('/merge/<acc_id>', methods=['POST'])
def finish_merge_job(acc_id):
    state = request.args.get('state')

    if state == "terminated":
        state = "merge_wait"

    if state not in ('merge_wait', 'merge_err', 'merge_done'):
        abort(400)

    session = db.get_session()
    acc = session.query(db.Accession)\
        .filter_by(acc_id=int(acc_id))\
        .one()

    if acc.state != 'merging':
        abort(400)

    acc.state = state
    session.commit()
    return jsonify({
        'result': 'success'
    })

