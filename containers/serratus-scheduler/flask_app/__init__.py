"""Master job scheduler.  Listens for requests from workers and feeds
them information about work that needs doing."""
import os
import json
import subprocess
import time

from flask import Flask, jsonify, request, current_app, redirect, url_for
from prometheus_client import Summary

from . import db, jobs, metrics, cron
from .metrics import REGISTRY

SCHEDULER_REQUEST_SUMMARY = Summary(
    "scheduler_request_latency_seconds",
    "Flask Request Latency",
    ["method", "endpoint"],
    registry=REGISTRY,
)


def before_request():
    request.start_time = time.time()


def after_request(response):
    request_latency = time.time() - request.start_time
    SCHEDULER_REQUEST_SUMMARY.labels(request.method, request.url_rule).observe(
        request_latency
    )

    return response


def wait_postgres(app):
    cmd = ["pg_isready", "-d", app.config["DATABASE"]]
    while True:
        ret = subprocess.call(cmd)
        if ret == 0:
            break
        elif ret in (1, 2):
            time.sleep(5)
            continue
        else:
            raise IOError("Invalid database command: {}".format(cmd))


def create_app(test_config=None):
    """Main entrypoint.  Run this program using:

    FLASK_APP=flask-app

    flask-3 run
    OR
    gunicorn 'flask-app:create_app()'
    """
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        DATABASE="postgresql://postgres@localhost:5432/postgres",
        AWS_REGION="us-east-1",
    )

    if test_config is None:
        app.config.from_pyfile("config.py", silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    @app.route("/status")
    def status():
        return jsonify({"status": "up"})

    @app.route("/config", methods=["GET", "PUT"])
    def config():
        if request.method == "PUT":
            db.update_config(json.loads(request.data))
        return jsonify(dict(db.get_config()))

    @app.route("/")
    def root():
        return redirect(url_for("jobs.show_jobs"))

    app.register_blueprint(db.bp)
    app.register_blueprint(jobs.bp)
    app.register_blueprint(metrics.bp)
    wait_postgres(app)

    db.init_app(app)
    cron.register(app)

    app.before_request(before_request)
    app.after_request(after_request)

    return app
