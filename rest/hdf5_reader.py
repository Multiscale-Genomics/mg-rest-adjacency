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

from dmp import dmp
from scripts.GenerateSampleAdjacency import GenerateSampleAdjacency

class hdf5:
    """
    Class related to handling the functions for interacting directly with the
    HDF5 files. All required information should be passed to this class.
    """

    def __init__(self, user_id, file_id, resolution=None):
        """
        Initialise the module and generate sample data if required

        Parameters
        ----------
        user_id : str
            Identifier to uniquely locate the users files. Can be set to
            "common" if the files can be shared between users
        file_id : str
            Location of the file in the file system
        resolution : int (Optional)
            Level of resolution. This is optional, but only the functions
            get_resolutions() and set_resolutions() can be called. Once the
            resolution has been set then all functions are callable.
        """
        self.user_id = user_id
        self.file_id = file_id

        if user_id == 'test':
            #resource_package = __name__
            #resource_path = os.path.join(os.path.dirname(__file__), 'rao2014.hdf5')
            resource_path = os.path.join(
                os.path.dirname(__file__),
                '../test/data/sample_adjacency.hdf5'
            )
            if os.path.isfile(resource_path) is False:
                gsa = GenerateSampleAdjacency()
                gsa.main()

            self.hdf5_handle = h5py.File(resource_path, "r")
        else:
            cnf_loc = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'mongodb.cnf'
            )
            dm_handle = dmp(cnf_loc)
            file_obj = dm_handle.get_file_by_id(user_id, file_id)
            self.hdf5_handle = h5py.File(file_obj["file_path"], "r")

        self.resolutions = map(int, self.hdf5_handle.keys())

        if resolution is None:
            self.dset = self.hdf5_handle[str(self.resolutions[0])]
            self.resolution = self.resolutions[0]
        else:
            self.dset = self.hdf5_handle[str(resolution)]
            self.resolution = resolution

        self.chr_param = self._calculate_chr_param(self.resolutions, self.dset.attrs["chromosomes"])

    def close(self):
        """
        Close the HDF5 data file handle
        """
        self.hdf5_handle.close()

    def get_resolutions(self):
        """
        List resolutions that models have been generated for

        Returns
        -------
        list : str
            Available levels of resolution that can be set
        """

        return [res for res in self.hdf5_handle]


    def set_resolution(self, resolution):
        """
        Set, or change, the resolution level

        Parameters
        ----------
        resolution : int
            Level of resolution
        """

        self.resolution = resolution

        self.dset = self.hdf5_handle[str(self.resolution)]
        chromosomes = self.dset.attrs['chromosomes']
        self.chr_param = self._calculate_chr_param(self.resolution, chromosomes)

    def get_details(self):
        """
        Return a list of the available resolutions in a given HDF5 file
        """

        return {
            "chromosomes" : [list(c) for c in self.dset.attrs["chromosomes"]],
            "chr_param"   : self.chr_param,
            "resolutions" : self.resolutions}

    def get_resolution(self):
        """
        List the current level of rseolution

        Returns
        -------
        resolution : int
            Current level of resolution
        """

        return self.resolution

    def get_chromosomes(self):
        """
        List of chromosomes that have models at a given resolution

        Returns
        -------
        chromosomes : list
            List of chromosomes at the set resolution
        """

        if self.resolution is None:
            return {}

        return [list(c) for c in self.dset.attrs['chromosomes']]

    def get_chromosome_parameters(self):
        """
        Return a list of the available resolutions in a given HDF5 file

        Returns
        -------
        dict
           chromosomes : list
           chr_param : dict
           resolitions

        Example
        -------

        .. code-block:: python
           :linenos:

           from reader import adjacency
           r = adjacency('test', '', 10000)
           value = r.get_chromosome_parameters()

        """

        return self.chr_param

    def get_range(
            self, chr_id, start, end,
            limit_chr=None, limit_start=None, limit_end=None,
            value_url='/api/getValue', no_links=None):
        """
        Get the interactions that happen within a defined region on a specific
        chromosome. Returns inter and intra interactions with the defined
        region.

        Parameters
        ----------
        chr_id : str
           Chromosomal name
        start : int
           Start position within the chromosome
        end : int
           End position within the chromosome
        limit_chr : str (Optional)
           Limit the results to a particular chromosome
        limit_start : int (Optional)
           Limit the range start position on the limit_chr paramter
        limit_end : int (Optional)
           Limit the range end position on the limit_chr parameter
        value_url : str (Optional)
           Define a custom URL snippet for the location of the file if different
           from the defaul
        no_links : bool (Optional)
           Will return the URL links to the individual points within the
           adjacency matrix. In cases where this generates a large number of
           points it is possible to turn off generating these links. Set this
           value to 1.

        Returns
        -------
        dict
           log : list
              List of messages about the state for debugging
           results : list
              List of values for given positions within the adjacency matrix

        Example
        -------

        .. code-block:: python
           :linenos:

           from reader import adjacency
           r = adjacency('test', '', 10000)
           value = r.get_range(2000000, 1000000)
        """

        # Defines columns to get extracted from the array
        x_pos = int(np.floor(float(start)/float(self.resolution)))
        y_pos = int(np.ceil(float(end)/float(self.resolution)))

        # xy_offset for the chromosome in the super array
        xy_offset = self.chr_param[chr_id]["bins"][self.resolution][1]

        dset = self.hdf5_handle[str(self.resolution)]

        start2 = 0
        end2 = 0
        if limit_chr != None:
            if limit_start != None and limit_end != None:
                start2 = int(np.floor(float(limit_start)/float(self.resolution)))
                end2 = int(np.ceil(float(limit_end)/float(self.resolution)))
                xy2_offset = self.chr_param[limit_chr]["bins"][self.resolution][1]

                result = dset[(x_pos+xy_offset):(y_pos+xy_offset), (start2+xy2_offset):(end2+xy2_offset)]
            else:
                start2 = self.chr_param[limit_chr]["bins"][self.resolution][1]
                end2 = start2 + self.chr_param[limit_chr]["bins"][self.resolution][0]

                result = dset[(x_pos+xy_offset):(y_pos+xy_offset), start2:end2]
        else:
            result = dset[(x_pos+xy_offset):(y_pos+xy_offset), :]

        # Iterate through slice and extract results greater than zero
        results = []
        log_text = []

        r_index = np.transpose(np.nonzero(result))
        log_text.append(
            {
                "coord" : {
                    "x0" : (x_pos+xy_offset),
                    "x1" : (y_pos+xy_offset)
                },
                "r_index" : len(r_index),
                "param": {
                    "start" : start,
                    "x" : x_pos,
                    "end" : end,
                    "y" : y_pos,
                    "xy_offset" : xy_offset,
                    "resolution" : self.resolution,
                    "chr_id" : chr_id,
                    "startB" : start2,
                    "endB" : end2,
                    "limit_chr" : limit_chr
                },
                'chr_param': self.chr_param
            }
        )

        for i in r_index:
            x_start = ((i[0]+x_pos)*int(self.resolution))
            y_chr = self.get_chromosome_from_array_index(i[1]+start2)
            if limit_chr != None:
                y_start = (i[1]+start2)*int(self.resolution)
            else:
                y_start = (i[1]-self.chr_param[y_chr]["bins"][self.resolution][1])*int(self.resolution)

            entry = {
                "chrA" : chr_id,
                "startA" : x_start,
                "chrB" : y_chr,
                "startB" : y_start,
                "value" : int(result[i[0], i[1]]),
                "pos_x" : i[0]+x_pos+xy_offset,
                "pos_y" : i[1]
            }
            if no_links is None:
                entry['_links'] = {'self': value_url + "?user_id=" + str(self.user_id) + "&file_id=" + str(self.file_id) + "&res=" + str(self.resolution) + "&pos_x=" + str(i[0]+x_pos+xy_offset) + "&pos_y=" + str(i[1])}
            results.append(entry)

        return {"log": log_text, "results": results}


    def get_value(self, bin_i, bin_j):
        """
        Get a specific value for a given dataset, resolution

        Parameters
        ----------
        bin_i : int
            Array position in the first dimension
        bin_j : int
            Array position in the second dimension

        Returns
        -------
        value : int
            Value for a given cell in the adjacency array

        Example
        -------

        .. code-block:: python
           :linenos:

           from reader import adjacency
           r = adjacency('test', '', 10000)
           value = r.get_value(2000000, 1000000)

        """
        value = self.dset[int(bin_i), int(bin_j)]
        return value

    def _calculate_chr_param(self, bin_sizes, chromosomes):
        """
        Load the self.chr_param object with the required information about the
        start and stop positions of the chromosomes and bins.
        """

        chr_param = {}

        genome_len = 0
        bin_count = [0]*len(bin_sizes)
        chromosome_count = len(chromosomes)
        for i in range(chromosome_count):
            chromosome = chromosomes[i]

            genome_len += int(chromosome[1])

            # Calculate the number of bins for a chromosome and then join with
            # the offset values for the start in the array
            bin_s = [int(np.ceil(int(chromosome[1])/float(y))) for y in bin_sizes]
            bin_c = dict(
                zip(
                    bin_sizes,
                    [[bin_s[j], bin_count[j], bin_s[j]+bin_count[j]] for j in range(len(bin_count))]
                )
            )

            chr_param[str(chromosome[0])] = {
                'size': [int(chromosome[1]), genome_len],
                'bins': bin_c
            }

            # Calculate the new offset values.
            bin_count = [bin_count[i]+bin_s[i] for i in range(len(bin_count))]

        total_bin_count = dict(
            zip(
                bin_sizes,
                [[bin_s[i], bin_count[i]] for i in range(len(bin_count))]
            )
        )
        chr_param["meta"] = {"genomeSize": genome_len, "totalBinCount": total_bin_count}

        return chr_param


    def get_chromosome_from_array_index(self, index):
        """
        Identify the chromosome based on either the x or y coordinate in the
        array.

        Parameters
        ----------
        index : int
            Location within the array

        Returns
        -------
        chr_id : str
            Identity of the chromosome

        Example
        -------

        .. code-block:: python
           :linenos:

           from reader import adjacency
           r = adjacency('test', '', 10000)
           cid = r.get_chromosome_from_array_index(1234567890)

        """
        for chr_id in self.chr_param.keys():
            if chr_id == "meta":
                continue
            #print(self.chr_param[chr_id]["bins"], type(self.chr_param[chr_id]["bins"]))
            chr_end = self.chr_param[chr_id]["bins"][int(self.resolution)][2]
            chr_start = self.chr_param[chr_id]["bins"][int(self.resolution)][1]
            if index >= chr_start and index <= chr_end:
                return chr_id
