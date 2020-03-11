import pytest
from flask_app import jobs
from json import loads

@pytest.fixture
def csv_preload(client):
    csv = b'a,b,c\n1,2,3\n4,5,6'
    response = client.post('/jobs/add_sra_run_info', data=csv)
    assert response.status_code == 200

def test_add_sra_run_info(client):
    csv = b'a,b,c\n1,2,3\n4,5,6'
    response = client.post('/jobs/add_sra_run_info', data=csv)
    assert response.status_code == 200
    assert loads(response.data)['inserted_rows'] == 2
    assert loads(response.data)['total_rows'] == 2

    response = client.post('/jobs/add_sra_run_info', data=csv)
    assert response.status_code == 200
    assert loads(response.data)['inserted_rows'] == 2
    assert loads(response.data)['total_rows'] == 4

    response = client.post('/jobs/split/')
    assert response.status_code == 200
    sra_data = loads(response.data)['sra_run_info'] == {'a': '1',
                                                        'b': '2',
                                                        'c': '3',
                                                       }

    response = client.post('/jobs/split/')
    assert response.status_code == 200
    sra_data = loads(response.data)['sra_run_info'] == {'a': '4',
                                                        'b': '5',
                                                        'c': '6',
                                                       }

def test_start_split_job(client, csv_preload):
    response = client.post('/jobs/split/1')
    # Don't care what code we get yet.  Just not success.
    assert response.status_code != 200

    response = client.post('/jobs/split/')
    assert response.status_code == 200
    data = loads(response.data)
    assert data['acc_id'] == 1

    # Missing status
    response = client.post('/jobs/split/{}'.format(data['acc_id']))
    assert response.status_code != 200

    response = client.post('/jobs/split/{}?status=split_done&N_paired=3'.format(data['acc_id']))
    assert response.status_code == 200

    # Repeated
    response = client.post('/jobs/split/{}?status=split_done&N_paired=3'.format(data['acc_id']))
    assert response.status_code != 200

def test_job_reset_commit(client, csv_preload):
    """When setting job to 'new' or 'split_err', state should be committed"""
    response = client.post('/jobs/split/')
    data = loads(response.data)
    response = client.post('/jobs/split/{}?status=new'.format(data['acc_id']))
    assert response.status_code == 200
    assert loads(response.data) == {'result': 'success'}

    response = client.post('/jobs/split/')
    data = loads(response.data)
    response = client.post('/jobs/split/{}?status=split_err'.format(data['acc_id']))
    assert response.status_code == 200
    assert loads(response.data) == {'result': 'success'}




