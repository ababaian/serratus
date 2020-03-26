"""Master job scheduler.  Listens for requests from workers and feeds
them information about work that needs doing."""
import os

from flask import Flask, jsonify
from . import db, jobs


def create_app(test_config=None):
    """Main entrypoint.  Run this program using:

    FLASK_APP=flask-app

    flask-3 run
    OR
    gunicorn 'flask-app:create_app()'
    """
    app = Flask(__name__, instance_relative_config=True)
    app.config.from_mapping(
        DATABASE=os.path.join(app.instance_path, 'scheduler.sqlite')
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

    app.register_blueprint(db.bp)
    app.register_blueprint(jobs.bp)

    db.init_app(app)
    return app
