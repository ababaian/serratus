import sqlite3

import click
from flask import current_app, g
from flask.cli import with_appcontext

ACC_STATES = ('new', 'splitting', 'split_done', 'merge_wait', 'merging',
              'merge_done', 'split_err', 'merge_err')

FQ_STATES = ('new', 'aligning', 'done', 'fail')

def get_db():
    if 'db' not in g:
        print('Creating new connection')
        g.db = sqlite3.connect(
            current_app.config['DATABASE'],
            detect_types=sqlite3.PARSE_DECLTYPES
        )
        g.db.row_factory = sqlite3.Row

    return g.db

def acc_states():
    """Return accession states as a dictionary"""
    if 'acc_states' not in g:
        g.acc_states = {}
        c = get_db().cursor()
        c.execute("SELECT acc_state_id,name FROM acc_state_defn")
        for id, name in c.fetchall():
            g.acc_states[name] = id

    return g.acc_states

def fastq_states():
    """Return fastq states as a dictionary"""
    if 'fq_states' not in g:
        g.fq_states = {}

        c = get_db().cursor()
        c.execute("SELECT fastq_state_id,name FROM fastq_state_defn")
        for id, name in c.fetchall():
            g.fq_states[name] = id

    return g.fq_states

def close_db(e=None):
    db = g.pop('db', None)

    if db is not None:
        db.close()

def init_db(job_csv):
    db = get_db()

    with current_app.open_resource('schema.sql') as f:
        db.executescript(f.read().decode('utf8'))

    cursor = db.cursor()
    for state in ACC_STATES:
        cursor.execute('INSERT INTO acc_state_defn VALUES (null, ?)', (state,))

    for state in FQ_STATES:
        cursor.execute('INSERT INTO fastq_state_defn VALUES (null, ?)', (state,))

    states = acc_states()
    import csv
    with open(job_csv) as f:
        for line in csv.reader(f):
            cursor.execute('INSERT INTO acc VALUES (null, ?, ?, ?, ?, ?)',
                           [states['new']] + line)

    db.commit()

@click.command('init-db')
@click.argument('job_csv', type=str)
@with_appcontext
def init_db_command(job_csv):
    """Clear the existing data and create new tables."""
    init_db(job_csv)
    click.echo('Initialized the database.')

def init_app(app):
    app.teardown_appcontext(close_db)
    app.cli.add_command(init_db_command)

