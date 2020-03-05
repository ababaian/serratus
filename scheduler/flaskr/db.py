import click
from flask import current_app, g
from flask.cli import with_appcontext

from sqlalchemy import create_engine
from sqlalchemy import Column, Integer, String, ForeignKey, Boolean, Enum
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()

ACC_STATES = ('new', 'splitting', 'split_done', 'merge_wait', 'merging',
              'merge_done', 'split_err', 'merge_err')

FQ_STATES = ('new', 'aligning', 'done', 'fail')

class Accession(Base):
    __tablename__ = 'acc'

    acc_id = Column(Integer, primary_key=True)
    acc_state = Column(Enum(*ACC_STATES))

    acc = Column(String)
    sra_url = Column(String)
    split_cmd = Column(String)
    merge_cmd = Column(String)

    def to_dict(self):
        return {
            'acc_id': self.acc_id,
            'acc': self.acc,
            'sra_url': self.sra_url,
            'split_cmd': self.split_cmd,
            'merge_cmd': self.merge_cmd,
        }

class FastQ(Base):
    __tablename__ = 'fastq'

    fastq_id = Column(Integer, primary_key=True)
    fastq_state = Column(Enum(*FQ_STATES))
    acc_id = Column(Integer, ForeignKey('acc.acc_id'))
    n = Column(Integer)
    align_cmd = Column(String)
    paired = Column(Boolean)

def get_engine(echo=False):
    path = 'sqlite:///' + current_app.config['DATABASE']
    print(path)
    return create_engine(path, echo=echo)

def get_session(**kwargs):
    return sessionmaker(bind=get_engine(**kwargs))()

@click.command('init-db')
@click.argument('job_csv', type=str)
@with_appcontext
def init_db_command(job_csv):
    """Clear the existing data and create new tables."""
    engine = get_engine(echo=False)
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)

    # Using sqlite3 instead of SQLAlchemy here because SQLAlchemy uses a
    # transaction, which is slowwwwwwww (2M lines = 300s vs 5s)
    import csv, sqlite3
    conn = sqlite3.connect(current_app.config['DATABASE'])
    cursor = conn.cursor()
    with open(job_csv) as f:
        for line in csv.reader(f):
            cursor.execute('INSERT INTO acc VALUES (null, ?, ?, ?, ?, ?)',
                           ['new'] + line)

    conn.commit()

    click.echo('Initialized the database.')

def init_app(app):
    app.cli.add_command(init_db_command)

