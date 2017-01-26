"""
Copyright 2017 EMBL-European Bioinformatics Institute

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
    db_fd, app.app.config['DATABASE'] = tempfile.mkstemp()
    app.app.config['TESTING'] = True
    client = app.app.test_client()
    
    def teardown():
        os.close(db_fd)
        os.unlink(app.app.config['DATABASE'])
    request.addfinalizer(teardown)

    return client

def test_taxon_id(client):
    rv = client.get('/rest/v0.0/getInteractions')
    assert b'9606' in rv.data

def test_accession_id(client):
    rv = client.get('/rest/v0.0/getInteractions/9606')
    assert b'GCA_000001405.14' in rv.data

def test_datasets(client):
    rv = client.get('/rest/v0.0/getInteractions/9606/GCA_000001405.14')
    assert b'rao2014' in rv.data

def test_resolutions(client):
    rv = client.get('/rest/v0.0/getInteractions/9606/GCA_000001405.14/rao2014')
    assert b'resolution' in rv.data
    assert b'10000' in rv.data

def test_chromosomes(client):
    rv = client.get('/rest/v0.0/getInteractions/9606/GCA_000001405.14/rao2014/10000')
    assert b'1' in rv.data
    assert b'X' in rv.data
    assert b'Y' in rv.data

def test_chromosome(client):
    rv = client.get('/rest/v0.0/getInteractions/9606/GCA_000001405.14/rao2014/10000/1')
    assert b'bins' in rv.data
    assert b'size' in rv.data

def test_range(client):
    rv = client.get('/rest/v0.0/getInteractions/9606/GCA_000001405.14/rao2014/10000/1/0/10000')
    assert b'startA' in rv.data
    assert b'startB' in rv.data

def test_value(client):
    rv = client.get('/rest/v0.0/getValue/9606/GCA_000001405.14/rao2014/10000/1/1')
    assert b'46' in rv.data
