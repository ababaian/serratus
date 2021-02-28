from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from psycogreen.gevent import patch_psycopg
from gevent import monkey

monkey.patch_all()
patch_psycopg()

DATABASE = "postgresql://postgres@localhost:5433/postgres"

app = Flask(__name__)
app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False
db = SQLAlchemy(app)


@app.route("/cpu")
def cpu():
    return str(sum(range(100000000)))


@app.route("/io")
def index():
    with db.engine.connect() as con:
        con.execute("SELECT pg_sleep(5);")
    return "hello there!\n"
