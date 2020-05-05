"""Master job scheduler.  Listens for requests from workers and feeds
them information about work that needs doing."""
import os
import json

from flask import Flask, jsonify, request, current_app, redirect, url_for

from . import db, jobs, metrics, cron

def create_app(test_config=None):
    """Main entrypoint.  Run this program using:

    FLASK_APP=flask-app

    flask-3 run
    OR
    gunicorn 'flask-app:create_app()'
    """
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        DATABASE=os.path.join(app.instance_path, 'scheduler.sqlite'),
        AWS_REGION="us-east-1",
    )

    if test_config is None:
        app.config.from_pyfile('config.py', silent=True)
    else:
        app.config.from_mapping(test_config)

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    @app.route('/status')
    def status():
        return jsonify({'status': 'up'})

    @app.route('/config', methods=["GET", "PUT"])
    def config():
        if request.method == "PUT":
            db.update_config(json.loads(request.data))
        return jsonify(dict(db.get_config()))

    @app.route('/')
    def root():
        return redirect(url_for('jobs.show_jobs'))

    app.register_blueprint(db.bp)
    app.register_blueprint(jobs.bp)
    app.register_blueprint(metrics.bp)
    db.init_app(app)

    cron.register(app)

    return app
