"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.

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
import json
from httplib2 import Http
from flask import request, abort

def validate_token(access_token):
    """
    Verify that a MuG access token is valid
    """

    with open(os.path.dirname(os.path.abspath(__file__)) + '/auth_meta.json') as data_file:
        data = json.load(data_file)

    http_handler = Http()
    resp, user_data = http_handler.request(
        data['auth_server']['url'],
        headers={
            'Authorization': access_token
        }
    )

    if not resp['status'] == '200':
        return None

    try:
        data = json.loads(user_data)
    except TypeError:
        # Python 3 returns byt objects
        data = json.loads(user_data.decode())

    return {
        'user_id': data['mug_id']
    }

def authorized(func):
    """
    Wrapper for authorization based on tokens
    """

    def _wrap(*args, **kwargs):
        if 'Authorization' not in request.headers:
            print('No token provided')
            abort(401)
            return None

        print('Checking token ...')
        user_id = validate_token(request.headers['Authorization'])
        if user_id is None:
            print('Check FAILED')
            abort(401)
            return None

        return func(user_id=user_id, *args, **kwargs)

    return _wrap
