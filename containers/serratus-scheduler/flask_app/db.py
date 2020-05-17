from flask import Blueprint, current_app, g, send_file, jsonify, request

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
        path = current_app.config["DATABASE"]
        engine.append(create_engine(path, echo=False))

    return engine[0]

def get_session():
    if 'session' not in g:
        g.session = sessionmaker(bind=get_engine())(expire_on_commit=False)
    return g.session


def teardown_session(err=None):
    session = g.pop('session', None)

    if session is not None:
        session.close()


def create_config(data):
    """Load default configuration values into the DB, iff they don't exist there yet"""
    session = get_session()
    for key, value in data.items():
        row = session.query(Config)\
            .filter_by(key=key)\
            .with_for_update()\
            .one_or_none()

        if row is None:
            session.add(Config(key=key, value=value))
            session.commit()


def update_config(data):
    """Load updated config into the DB.  If keys don't exist, raise a ValueError()"""
    session = get_session()
    for key, value in data.items():
        row = session.query(Config)\
            .filter_by(key=key)\
            .one_or_none()
        if row is None:
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


def init_db(reset=False):
    """Clear the existing data and create new tables."""
    engine = get_engine()
    if reset:
        Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    create_config(CONFIG_DEFAULT)

@bp.route('', methods=['GET'])
def dump_db_sqlite():
    """Get a full copy of the SQLite database file"""
    # TODO: Convert -> CSV then load that into a newly created SQLite DB.
    raise NotImplementedError()
    #return send_file(database_file,
    #                 mimetype="application/x-sqlite3",
    #                 as_attachment=True,
    #                 attachment_filename='dump.sqlite')

@bp.route('', methods=['PUT'])
def load_db_sqlite():
    """Replace the DB with a given SQLite file"""
    init_db(reset=True)

    if not request.data:
        return jsonify('Database reset to initial state')

    raise NotImplementedError()
    # TODO: Convert -> CSV then load the CSV into postgres.
    #with open(current_app.config['DATABASE'], 'wb') as f:
    #    f.write(request.data)

    return jsonify("Database replaced")

def init_app(app):
    app.teardown_appcontext(teardown_session)
