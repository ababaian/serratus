import click
from flask import Blueprint, current_app, g, send_file, jsonify, request
from flask.cli import with_appcontext

from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, JSON, ForeignKey, Boolean, Enum, DateTime, String
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

ACC_STATES = ('new', 'splitting', 'split_done', 'merge_wait', 'merging',
              'merge_done', 'split_err', 'merge_err')

bp = Blueprint('db', __name__, url_prefix='/db')

class Dicter():
    """Give Base class option to output dicts"""
    __table__ = None
    def to_dict(self):
        ret = {}
        for col in self.__table__.columns:
            ret[col.name] = getattr(self, col.name)

        return ret

class Accession(Base, Dicter):
    __tablename__ = 'acc'

    acc_id = Column(Integer, primary_key=True)
    state = Column(Enum(name='state', *ACC_STATES))

    contains_paired = Column(Boolean)
    contains_unpaired = Column(Boolean)
    blocks = Column(Integer)
    sra_run_info = Column(JSON)

    split_start_time = Column(DateTime)
    split_end_time = Column(DateTime)
    split_worker = Column(String)

    merge_start_time = Column(DateTime)
    merge_end_time = Column(DateTime)
    merge_worker = Column(String)

BLOCK_STATES = ('new', 'aligning', 'done', 'fail')

class Block(Base, Dicter):
    __tablename__ = 'blocks'

    block_id = Column(Integer, primary_key=True)
    state = Column(Enum(name='state', *BLOCK_STATES))
    acc_id = Column(Integer, ForeignKey('acc.acc_id'))
    n = Column(Integer)

    align_start_time = Column(DateTime)
    align_end_time = Column(DateTime)
    align_worker = Column(String)

def get_engine(echo=False, engine=[]):
    if not engine:
        path = 'sqlite:///' + current_app.config['DATABASE']
        engine.append(create_engine(path, echo=echo))

    return engine[0]

def get_session():
    if 'session' not in g:
        g.session = sessionmaker(bind=get_engine())()
    return g.session

def teardown_session(e=None):
    session = g.pop('session', None)

    if session is not None:
        session.close()

def init_db():
    """Clear the existing data and create new tables."""
    engine = get_engine(echo=False)
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)

@bp.route('', methods=['GET'])
def dump_db_sqlite():
    """Get a full copy of the SQLite database file"""
    return send_file(current_app.config['DATABASE'],
                     mimetype="application/x-sqlite3",
                     as_attachment=True,
                     attachment_filename='dump.sqlite')

@bp.route('', methods=['PUT'])
def load_db_sqlite():
    """Replace the DB with a given SQLite file"""
    if not request.data:
        init_db()
        return jsonify('Database reset to initial state')

    # Do I have to delete and replace my engine object?
    with open(current_app.config['DATABASE'], 'wb') as f:
        f.write(request.data)

    return jsonify("Database replaced")

@click.command('init-db')
@with_appcontext
def init_db_command():
    init_db()
    click.echo('Initialized the database.')

def init_app(app):
    app.teardown_appcontext(teardown_session)
    app.cli.add_command(init_db_command)
