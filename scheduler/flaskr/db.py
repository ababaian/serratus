import click
from flask import current_app, g
from flask.cli import with_appcontext

from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, ForeignKey, Boolean
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

ACC_STATES = ('new', 'splitting', 'split_done', 'merge_wait', 'merging',
              'merge_done', 'split_err', 'merge_err')

FQ_STATES = ('new', 'aligning', 'done', 'fail')

class AccessionState(Base):
    __tablename__ = 'acc_state_defn'

    acc_state_id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

class FastQState(Base):
    __tablename__ = 'fastq_state_defn'

    fastq_state_id = Column(Integer, primary_key=True)
    name = Column(String, unique=True)

class Accession(Base):
    __tablename__ = 'acc'

    acc_id = Column(Integer, primary_key=True)
    acc_state_id = Column(Integer, ForeignKey('acc_state_defn.acc_state_id'))

    acc = Column(String)
    sra_url = Column(String)
    split_cmd = Column(String)
    merge_cmd = Column(String)

class FastQ(Base):
    __tablename__ = 'fastq'

    fastq_id = Column(Integer, primary_key=True)
    fastq_state_id = Column(Integer,
                            ForeignKey('fastq_state_defn.fastq_state_id'))
    acc_id = Column(Integer, ForeignKey('acc.acc_id'))
    n = Column(Integer)
    align_cmd = Column(String)
    paired = Column(Boolean)

def get_engine(echo=True):
    path = 'sqlite:///' + current_app.config['DATABASE']
    print(path)
    return create_engine(path, echo=echo)

def get_session(**kwargs):
    if 'session' not in g:
        g.session = sessionmaker(bind=get_engine(**kwargs))()

    return g.session

def acc_states():
    """Return accession states as a dictionary"""
    if 'acc_states' not in g:
        g.acc_states = {}
        session = get_session()
        for row in session.query(AccessionState).all():
            g.acc_states[row.name] = row.acc_state_id

    return g.acc_states

def fastq_states():
    """Return fastq states as a dictionary"""
    if 'fq_states' not in g:
        g.fq_states = {}
        session = get_session()
        for row in session.query(FastQState).all():
            g.fq_states[row.name] = row.acc_state_id

    return g.fq_states

@click.command('init-db')
@click.argument('job_csv', type=str)
@with_appcontext
def init_db_command(job_csv):
    """Clear the existing data and create new tables."""
    engine = get_engine(echo=False)
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)

    session = get_session(echo=False)
    for state in ACC_STATES:
        session.add(AccessionState(name=state))

    for state in FQ_STATES:
        session.add(FastQState(name=state))

    states = acc_states()
    import csv
    with open(job_csv) as f:
        for line in csv.reader(f):
            acc = Accession(acc_state_id=states['new'],
                            acc=line[0],
                            sra_url=line[1],
                            split_cmd=line[2],
                            merge_cmd=line[3])
            session.add(acc)

    session.commit()
    click.echo('Initialized the database.')

def init_app(app):
    app.cli.add_command(init_db_command)

