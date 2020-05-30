from flask import Blueprint, make_response
import prometheus_client
from prometheus_client import Gauge
from . import db

bp = Blueprint("metrics", __name__, url_prefix="/metrics")


@bp.route("")
def metrics():
    r = make_response((prometheus_client.generate_latest(), 200))
    r.mimetype = "text/plain"
    return r


def acc_count(state):
    session = db.get_session()
    return session.query(db.Accession).filter_by(state=state).count()
    session.close()


def block_count(state):
    session = db.get_session()
    return session.query(db.Block).filter_by(state=state).count()
    session.close()


def register():
    """Create gauges for the database and register their callbacks"""
    acc_state_gauge = Gauge(
        "scheduler_accesion_state", "Number of accessions in each state", ["state"]
    )

    block_state_gauge = Gauge(
        "scheduler_block_state", "Number of blocks in each state", ["state"]
    )

    for state in db.ACC_STATES:
        count = lambda state=state: acc_count(state)
        acc_state_gauge.labels(state=state).set_function(count)

    for state in db.BLOCK_STATES:
        count = lambda state=state: block_count(state)
        block_state_gauge.labels(state=state).set_function(count)


register()
