from flask import Blueprint, jsonify, request, abort
from . import db

bp = Blueprint('jobs', __name__, url_prefix='/jobs')

@bp.route('/reset_db', methods=['POST'])
def reset_db():
    db.init_db()
    return jsonify('Database cleared')

@bp.route('/add_sra_run_info', methods=['POST'])
def add_sra_runinfo():
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
        return jsonify({'action': 'shutdown'})

    acc.state = 'splitting'
    session.add(acc)
    response = acc.to_dict()
    session.commit()
    response['action'] = 'process'

    # Send the response as JSON
    return jsonify(response)

@bp.route('/split/<acc_id>', methods=['POST'])
def finish_split_job(acc_id):
    """Mark a split job as finished, setting off N alignment jobs.

    param acc_id: The acc_id, as returned by start_split_job
    param status: Resulting status, see sched_states.dia
    param N_paired: Number of paired blocks the split resulted in
    param N_unpaired: Number of unpaired blocks the split resulted in

    Examples:

    Successful split into 10 blocks:
        http://scheduler/jobs/split/1?status=split_done&N_paired=10

    Failed split (bug in code or bad data):
        http://scheduler/jobs/split/1?status=split_err

    Terminated node (please redo this work):
        http://scheduler/jobs/split/1?status=new
    """
    status = request.args.get('status')
    try:
        number_split = int(request.args.get('N_unpaired', 0))
        number_split = int(request.args.get('N_paired', 0))
    except TypeError:
        abort(400)

    # Update the accessions table
    if status not in ('new', 'split_err', 'split_done'):
        abort(400)

    session = db.get_session()
    acc = session.query(db.Accession)\
        .filter_by(acc_id=int(acc_id))\
        .one()

    if acc.state != 'splitting':
        abort(400)
    acc.state = status

    if status in ('new', 'split_err'):
        # Not ready to process more
        session.commit()
        return jsonify({'result': 'success'})

    # Insert N align jobs into the alignment table.
    for i in range(int(number_split)):
        chunk = db.Block(state='new', acc_id=acc.acc_id, n=i)
        session.add(chunk)
    session.commit()

    return jsonify('inserted {} jobs'.format(number_split))

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
        return jsonify({'action': 'shutdown'})

    chunk, acc = query

    chunk.state = 'aligning'
    session.add(chunk)

    response = chunk.to_dict()
    response.update(acc.to_dict())
    session.commit()

    # Send the response as JSON
    return jsonify(response)

@bp.route('/align/<block_id>', methods=['POST'])
def finish_align_job():
    """Finished job, chunk_id is the same parameter from the start_job,
    state is one of (new, done, fail)"""
    chunk_id = request.args.get('chunk_id')
    state = request.args.get('state')

    if state not in db.CHUNK_STATES:
        abort(400)

    session = db.get_session()
    chunk = session.query(db.Block).filter(chunk_id=chunk_id).one()
    chunk.state = state

    return 'success'

@bp.route('/merge/', methods=['POST'])
def start_merge_job():
    pass


@bp.route('/merge/<acc_id>', methods=['POST'])
def finish_merge_job(acc_id):
    pass



