from flask_app import create_app
from json import loads

def test_config():
    assert not create_app().testing
    assert create_app({'TESTING': True}).testing

def test_status(client):
    response = client.get('/status')
    assert loads(response.data)['status'] == 'up'
