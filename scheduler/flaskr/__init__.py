"""Master job server"""
import os

from flask import Flask, request, current_app, jsonify

from . import db

def add_endpoints(app):
    @app.route('/start_split_job')
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
            return 'shutdown\n'

        acc.state = 'splitting'
        session.add(acc)

        response = acc.to_dict()

        session.commit()
        response['action'] = 'process'

        # Send the response as JSON
        return jsonify(response)

    @app.route('/finish_split_job')
    def finish_split_job():
        """Mark a split job as finished, setting off N alignment jobs.

        param job_id: The acc_id, as returned by start_split_job
        param status: Resulting status, see sched_states.dia
        param N: Number of blocks the file was split into

        Examples:

        Successful split into 10 blocks:
            http://scheduler/finish_split_job?job_id=1&status=split_done&N=10

        Failed split (bug in code or bad data):
            http://scheduler/finish_split_job?job_id=1&status=split_err

        Terminated node (please redo this work):
            http://scheduler/finish_split_job?job_id=1&status=new
        """
        job_id = request.args.get('job_id')
        status = request.args.get('status')
        number_split = request.args.get('N')

        # Update the accessions table
        if status not in ('new', 'split_err', 'split_done'):
            raise ValueError()

        session = db.get_session()
        acc = session.query(db.Accession)\
            .filter_by(acc_id=int(job_id))\
            .one()

        if acc.state != 'splitting':
            raise ValueError('Job was not in progress')
        acc.state = status

        if status in ('new', 'split_err'):
            # Not ready to process more
            return 'updated db'

        # Insert N align jobs into the alignment table.
        for i in range(int(number_split)):
            chunk = db.Chunk(state='new', acc_id=acc.acc_id, n=i)
            session.add(chunk)
        session.commit()

        return jsonify('inserted {} jobs'.format(number_split))

    @app.route('/start_align_job')
    def start_align_job():
        """Get a job id and mark it as "running" in the DB

        returns: a job JSON file"""
        session = db.get_session()
        # get an item where state = new and update its state
        query = session.query(db.Chunk, db.Accession)\
            .filter(db.Chunk.state == 'new')\
            .filter(db.Chunk.acc_id == db.Accession.acc_id)\
            .first()

        if query is None:
            # TODO: If we got no results, then could be:
            #    * Waiting for split nodes: hang tight
            #    * No more work: shutdown
            return 'nothing to do\n'

        chunk, acc = query

        chunk.state = 'aligning'
        session.add(chunk)

        response = chunk.to_dict()
        response.update(acc.to_dict())
        session.commit()

        # Send the response as JSON
        return jsonify(response)

    @app.route('/finish_align_job')
    def finish_align_job():
        """Finished job, chunk_id is the same parameter from the start_job,
        state is one of (new, done, fail)"""
        chunk_id = request.args.get('chunk_id')
        state = request.args.get('state')

        if state not in db.CHUNK_STATES:
            raise ValueError()

        session = db.get_session()
        chunk = session.query(db.Chunk).filter(chunk_id=chunk_id).one()
        chunk.state = state

        return 'success'

def create_app(test_config=None):
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        SECRET_KEY='dev',
        DATABASE=os.path.join(app.instance_path, 'flaskr.sqlite')
    )

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    add_endpoints(app)

    db.init_app(app)
    return app

