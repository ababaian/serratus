import pytest
from flask_app import db

def test_get_engine_cache(app):
    """There should be one database engine per app"""
    with app.app_context():
        engine = db.get_engine()
        assert bool(engine)

    with app.app_context():
        assert engine is db.get_engine()

def test_get_session_cache(app):
    with app.app_context():
        session = db.get_session()
        assert session is db.get_session()

    with app.app_context():
        assert session is not db.get_session()

def test_get_session_commit(app):
    with app.app_context():
        session = db.get_session()
        session.add(db.Accession())
        # No commit
        assert session.query(db.Accession).count() == 1

    with app.app_context():
        session = db.get_session()
        assert session.query(db.Accession).count() == 0
        session.add(db.Accession())
        session.commit()

    with app.app_context():
        assert session.query(db.Accession).count() == 1

def test_get_session_error(app):
    with pytest.raises(ValueError):
        with app.app_context():
            session = db.get_session()
            session.add(db.Accession())
            assert session.query(db.Accession).count() == 1
            raise ValueError()
            # Should rollback

    with app.app_context():
        assert db.get_session().query(db.Accession).count() == 0

def test_init_db(app):
    "init_db should reset the database"
    with app.app_context():
        session = db.get_session()
        session.add(db.Accession())
        session.commit()
    with app.app_context():
        db.init_db()
        assert db.get_session().query(db.Accession).count() == 0
        assert db.get_session().query(db.Block).count() == 0

def test_init_db_command(runner, monkeypatch):
    "Command should call init_db()"
    class Recorder():
        called = False

    def fake_init_db():
        Recorder.called = True

    monkeypatch.setattr('flask_app.db.init_db', fake_init_db)
    result = runner.invoke(args=['init-db'])
    assert Recorder.called

