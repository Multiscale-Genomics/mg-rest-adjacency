#!/usr/bin/env python

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
import h5py
import numpy as np

def create_matrix(size):
    """
    Function to generate a matrix of size n with random integers to populate the
    grid
    """
    rand_matrix = np.reshape(
        np.random.choice(
            [0, 1], size*size, p=[0.9, 0.1]
        ),
        (size, size)
    )
    return rand_matrix

def main():
    """
    Main function
    """
    resolutions = [10000, 100000, 1000000]
    chromosomes = [
        ['chr1', 32000000],
        ['chr2', 16000000],
        ['chr3', 8000000],
        ['chr4', 4000000],
        ['chr5', 2000000],
        ['chr6', 1000000],
        ['X', 10000000]
    ]

    d_size = sum([c[1] for c in chromosomes])

    # Create the HDF5 file
    filename = os.path.join(os.path.dirname(__file__), "../test/data/sample_adjacency.hdf5")
    hdf5_handle = h5py.File(filename, "w")

    for resolution in resolutions:
        local_size = d_size/resolution
        print(resolution, d_size, local_size)
        d_sample = np.zeros([local_size, local_size], dtype='int32')
        d_sample += create_matrix(local_size)

        dset = hdf5_handle.create_dataset(
            str(resolution),
            (local_size, local_size),
            dtype='int32',
            chunks=True,
            compression="gzip"
        )
        dset.attrs['chromosomes'] = chromosomes
        dset[0:local_size, 0:local_size] += d_sample

    hdf5_handle.close()


if __name__ == '__main__':
    main()
