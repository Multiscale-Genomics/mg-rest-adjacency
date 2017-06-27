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

from __future__ import print_function

import os
import tempfile
import json
import pytest

from context import app

@pytest.fixture
def client(request):
    """
    Definges the client object to make requests against
    """
    db_fd, app.APP.config['DATABASE'] = tempfile.mkstemp()
    app.APP.config['TESTING'] = True
    client = app.APP.test_client()
    
    def teardown():
        """
        Close the client once testing has completed
        """
        os.close(db_fd)
        os.unlink(app.APP.config['DATABASE'])
    request.addfinalizer(teardown)

    return client

def test_endpoints(client):
    """
    Test that the root endpoint is returning the expected keys
    """
    rest_value = client.get('/mug/api/adjacency')
    print(rest_value.data)
    assert '_links' in rest_value.data

def test_ping(client):
    """
    Test that the ping function is returning the status key
    """
    rest_value = client.get('/mug/api/adjacency/ping')
    print(rest_value.data)
    assert 'status' in rest_value.data

def test_details(client):
    """
    Test the details endpoint to ensure it returns the usage element
    """
    rest_value = client.get('/mug/api/adjacency/details')
    details = json.loads(rest_value.data)
    print(details)
    assert 'usage' in details

def test_details_00(client):
    """
    Test that details endpoint includes a resolution value when a test user and
    file are specified
    """
    rest_value = client.get('/mug/api/adjacency/details?user_id=test&file_id=test')
    details = json.loads(rest_value.data)
    print(details)
    assert 'resolutions' in details

def test_details_01(client):
    """
    Test that details endpoint includes 3 resolution values when a test user and
    file are specified
    """
    rest_value = client.get('/mug/api/adjacency/details?user_id=test&file_id=test')
    details = json.loads(rest_value.data)
    print(details)
    assert len(details['resolutions']) == 3

def test_getinteractions(client):
    """
    Test the interactions endpoint to ensure it returns the usage element
    """
    rest_value = client.get('/mug/api/adjacency/getInteractions')
    details = json.loads(rest_value.data)
    print(details)
    assert 'usage' in details

def test_getinteractions_00(client):
    """
    Test that interactions returns a values block when test parameters are
    provided
    """
    rest_value = client.get('/mug/api/adjacency/getInteractions?user_id=test&file_id=test&chr=19&res=1000000&start=1000000&end=5000000')
    details = json.loads(rest_value.data)
    print(details.keys())
    assert 'values' in details

def test_getinteractions_01(client):
    """
    Test that interactions returns a known number of values when test parameters
    are provided
    """
    rest_value = client.get('/mug/api/adjacency/getInteractions?user_id=test&file_id=test&chr=19&res=1000000&start=1000000&end=5000000')
    details = json.loads(rest_value.data)
    print(len(details['values']))
    assert len(details['values']) == 11408

def test_getvalue(client):
    """
    Test the values endpoint to ensure it returns the usage element
    """
    rest_value = client.get('/mug/api/adjacency/getValue')
    details = json.loads(rest_value.data)
    print(details)
    assert 'usage' in details

def test_getvalue_00(client):
    """
    Test that values returns with a value element when test parameters are
    provided
    """
    rest_value = client.get('/mug/api/adjacency/getValue?user_id=test&file_id=test&res=1000000&pos_x=0&pos_y=2673')
    details = json.loads(rest_value.data)
    print(details)
    assert 'values' in details
