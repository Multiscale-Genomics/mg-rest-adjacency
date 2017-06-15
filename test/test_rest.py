"""
.. Copyright 2017 EMBL-European Bioinformatics Institute
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
   
       http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import os
import pytest
import tempfile
import json

from context import app

@pytest.fixture
def client(request):
    db_fd, app.APP.config['DATABASE'] = tempfile.mkstemp()
    app.APP.config['TESTING'] = True
    client = app.APP.test_client()
    
    def teardown():
        os.close(db_fd)
        os.unlink(app.APP.config['DATABASE'])
    request.addfinalizer(teardown)

    return client

def test_endpoints(client):
    rv = client.get('/mug/api/adjacency')
    print rv.data
    assert b'_links' in rv.data

def test_ping(client):
    rv = client.get('/mug/api/adjacency/ping')
    print rv.data
    assert b'status' in rv.data

def test_details(client):
    rv = client.get('/mug/api/adjacency/details')
    details = json.loads(rv.data)
    print details
    assert 'usage' in details

def test_details_00(client):
    rv = client.get('/mug/api/adjacency/details?user_id=test&file_id=test')
    details = json.loads(rv.data)
    print details
    assert 'resolutions' in details

def test_details_01(client):
    rv = client.get('/mug/api/adjacency/details?user_id=test&file_id=test')
    details = json.loads(rv.data)
    print details
    assert len(details['resolutions']) == 3

def test_getInteractions(client):
    rv = client.get('/mug/api/adjacency/getInteractions')
    details = json.loads(rv.data)
    print details
    assert 'usage' in details

def test_getInteractions_00(client):
    rv = client.get('/mug/api/adjacency/getInteractions?user_id=test&file_id=test&chr=19&res=1000000&start=1000000&end=5000000')
    details = json.loads(rv.data)
    print details.keys()
    assert 'values' in details

def test_getInteractions_01(client):
    rv = client.get('/mug/api/adjacency/getInteractions?user_id=test&file_id=test&chr=19&res=1000000&start=1000000&end=5000000')
    details = json.loads(rv.data)
    print len(details['values'])
    assert len(details['values']) == 11408

def test_getValue(client):
    rv = client.get('/mug/api/adjacency/getValue')
    details = json.loads(rv.data)
    print details
    assert 'usage' in details

def test_getValue_00(client):
    rv = client.get('/mug/api/adjacency/getValue?user_id=test&file_id=test&res=1000000&pos_x=0&pos_y=2673')
    details = json.loads(rv.data)
    print details
    assert 'values' in details
