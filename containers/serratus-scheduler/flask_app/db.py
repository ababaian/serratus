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
    state = Column(Enum(name='acc_state', *ACC_STATES), index=True)

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
    state = Column(Enum(name='blk_state', *BLOCK_STATES), index=True)
    acc_id = Column(Integer, ForeignKey('acc.acc_id'))
    n = Column(Integer)

    align_start_time = Column(DateTime)
    align_end_time = Column(DateTime)
    align_worker = Column(String)

CONFIG_DEFAULT = dict(
    GENOME="cov2m",
    CLEAR_INTERVAL=600,
    SCALING_INTERVAL=300,
    DL_ARGS="",
    DL_SCALING_ENABLE=True,
    DL_SCALING_CONSTANT=.1,
    DL_SCALING_MAX=10,
    ALIGN_ARGS="--very-sensitive-local",
    ALIGN_SCALING_ENABLE=True,
    ALIGN_SCALING_CONSTANT=0.75,
    ALIGN_SCALING_MAX=10,
    MERGE_ARGS="",
    MERGE_SCALING_ENABLE=True,
    MERGE_SCALING_CONSTANT=.1,
    MERGE_SCALING_MAX=3,
    VIRTUAL_ASG_MAX_INCREASE=10,
    VIRTUAL_SCALING_INTERVAL=60,
)

class Config(Base):
    __tablename__ = "config"

    key = Column(String, primary_key=True)
    value = Column(JSON)

def get_engine(engine=[]):
    if not engine:
        path = 'postgresql://postgres@localhost:5432/postgres'
        engine.append(create_engine(path, echo=False))

    return engine[0]

def get_session():
    if 'session' not in g:
        g.session = sessionmaker(bind=get_engine())(expire_on_commit=False)
    return g.session


def teardown_session(e=None):
    session = g.pop('session', None)

    if session is not None:
        session.close()


def update_config(data, create=True):
    session = get_session()
    for key, value in data.items():
        row = session.query(Config).filter_by(key=key).one_or_none()
        if row is None:
            if create:
                row = Config(key=key)
            else:
                raise ValueError("Invalid key: {}".format(key))
        row.value = value
        session.add(row)
    session.commit()


def get_config():
    session = get_session()
    rows = session.query(Config).all()
    for row in rows:
        yield row.key, row.value

def get_config_val(key):
    session = get_session()
    row = session.query(Config).filter_by(key=key).one()
    return row.value


def init_db():
    """Clear the existing data and create new tables."""
    engine = get_engine()
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    update_config(CONFIG_DEFAULT, create=True)

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
