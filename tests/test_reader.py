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

#import os
#import tempfile
#import json
import pytest

from reader.hdf5_adjacency import adjacency

def test_range():
    """
    Test the range function
    """
    hdf5_handle = adjacency('test', '', 10000)
    x = hdf5_handle.get_range('chr1', 100000, 200000, limit_chr='chr2')
    print(x)
    assert 'results' in x
    assert len(x['results']) > 0
