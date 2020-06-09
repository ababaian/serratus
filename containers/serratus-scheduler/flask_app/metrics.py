from flask import Blueprint, make_response
import prometheus_client
from prometheus_client import Gauge
from sqlalchemy import func

from . import db

bp = Blueprint("metrics", __name__, url_prefix="/metrics")

acc_state_gauge = Gauge(
    "scheduler_accesion_state", "Number of accessions in each state", ["state"]
)

block_state_gauge = Gauge(
    "scheduler_block_state", "Number of blocks in each state", ["state"]
)


def update_counts():
    """Update counters for prometheus"""
    session = db.get_session()

    acc_bystate = (
        session.query(db.Accession.state, func.count(db.Accession.state))
        .group_by(db.Accession.state)
        .all()
    )
    for state, count in acc_bystate:
        acc_state_gauge.labels(state).set(count)

    block_bystate = (
        session.query(db.Block.state, func.count(db.Block.state))
        .group_by(db.Block.state)
        .all()
    )
    for state, count in block_bystate:
        block_state_gauge.labels(state).set(count)


@bp.route("")
def metrics():
    update_counts()
    r = make_response((prometheus_client.generate_latest(), 200))
    r.mimetype = "text/plain"
    return r
