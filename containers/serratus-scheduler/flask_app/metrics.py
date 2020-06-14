from flask import Blueprint, make_response
import prometheus_client
from prometheus_client import multiprocess, CollectorRegistry

bp = Blueprint("metrics", __name__, url_prefix="/metrics")

REGISTRY = CollectorRegistry()
multiprocess.MultiProcessCollector(REGISTRY)


@bp.route("")
def metrics():
    r = make_response((prometheus_client.generate_latest(REGISTRY), 200))
    r.mimetype = "text/plain"
    return r
