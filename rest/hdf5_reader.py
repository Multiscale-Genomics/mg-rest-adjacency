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

    def __init__(self, user_id, file_id):
        """
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

    def close(self):
        """
        """
        self.hdf5_handle.close()

    def get_details(self):
        """
        Return a list of the available resolutions in a given HDF5 file
        """

        resolutions = self.hdf5_handle.keys()
        dset = self.hdf5_handle[str(resolutions[0])]
        chr_param = self._calculate_chr_param(resolutions, dset.attrs["chromosomes"])

        return {
            "chromosomes" : [list(c) for c in dset.attrs["chromosomes"]],
            "chr_param"   : chr_param,
            "resolutions" : resolutions}


    def get_range(
            self, resolution, chr_id, start, end,
            limit_chr=None, limit_start=None, limit_end=None,
            value_url='/api/getValue', no_links=None):
        """
        Get the interactions that happen within a defined region on a specific
        chromosome. Returns inter and intra interactions with the defined region.
        """

        # Defines columns to get extracted from the array
        x = int(np.floor(float(start)/float(resolution)))
        y = int(np.ceil(float(end)/float(resolution)))

        # Get meta data from HDF5 file
        resolutions = map(int, self.hdf5_handle.keys())
        dset = self.hdf5_handle[str(resolutions[0])]
        chr_param = self._calculate_chr_param(resolutions, dset.attrs["chromosomes"])

        # xy_offset for the chromosome in the super array
        xy_offset = chr_param[chr_id]["bins"][resolution][1]

        dset = self.hdf5_handle[str(resolution)]

        startB = 0
        endB = 0
        if limit_chr != None:
            if limit_start != None and limit_end != None:
                startB = int(np.floor(float(limit_start)/float(resolution)))
                endB = int(np.ceil(float(limit_end)/float(resolution)))
                xyB_offset = chr_param[limit_chr]["bins"][resolution][1]

                result = dset[(x+xy_offset):(y+xy_offset),(startB+xyB_offset):(endB+xyB_offset)]
            else:
                startB = chr_param[limit_chr]["bins"][resolution][1]
                endB = startB + chr_param[limit_chr]["bins"][resolution][0]

                result = dset[(x+xy_offset):(y+xy_offset),startB:endB]
        else:
            result = dset[(x+xy_offset):(y+xy_offset),:]

        # Iterate through slice and extract results greater than zero
        results = []
        rShape = result.shape
        logText = []

        r_index = np.transpose(np.nonzero(result))
        logText.append(
            {
                "coord" : {
                    "x0" : (x+xy_offset),
                    "x1" : (y+xy_offset)
                },
                "r_index" : len(r_index),
                "param": {
                    "start" : start,
                    "x" : x,
                    "end" : end,
                    "y" : y,
                    "xy_offset" : xy_offset,
                    "resolution" : resolution,
                    "chr_id" : chr_id,
                    "startB" : startB,
                    "endB" : endB,
                    "limit_chr" : limit_chr
                },
                'chr_param': chr_param
            }
        )

        for i in r_index:
            x_start = ((i[0]+x)*int(resolution))
            y_chr = self.get_chromosome_from_array_index(chr_param, int(resolution), i[1]+startB)
            if limit_chr != None:
                y_start = (i[1]+startB)*int(resolution)
            else:
                y_start = (i[1]-chr_param[y_chr]["bins"][resolution][1])*int(resolution)

            r = {
                "chrA" : chr_id,
                "startA" : x_start,
                "chrB" : y_chr,
                "startB" : y_start,
                "value": int(result[i[0], i[1]])
            }
            if no_links is None:
                r['_links'] = {'self': value_url + "?user_id=" + str(self.user_id) + "&file_id=" + str(self.file_id) + "&res=" + str(resolution) + "&pos_x=" + str(i[0]+x+xy_offset) + "&pos_y=" + str(i[1])}
            results.append(r)

        return {"log": logText, "results": results}


    def get_value(self, resolution, bin_i, bin_j):
        """
        Get a specific value for a given dataset, resoltuoin
        """
        dset = self.hdf5_handle[str(resolution)]

        value = dset[int(bin_i), int(bin_j)]

        return value

    def _calculate_chr_param(self, binSizes, chromosomes):
        """
        Load the self.chr_param object with the required information about the
        start and stop positions of the chromosomes and bins.
        """

        chr_param = {}

        genomeLen = 0
        binCount = [0]*len(binSizes)
        chromosome_count = len(chromosomes)
        for i in range(chromosome_count):
            c = chromosomes[i]

            genomeLen += int(c[1])

            # Calculate the number of bins for a chromosome and then join with
            # the offset values for the start in the array
            binS = [int(np.ceil(int(c[1])/float(y))) for y in binSizes]
            binC = dict(
                zip(
                    binSizes,
                    [[binS[j], binCount[j], binS[j]+binCount[j]] for j in range(len(binCount))]
                )
            )

            chr_param[str(c[0])] = {'size': [int(c[1]), genomeLen], 'bins': binC}

            # Calculate the new offset values.
            binCount = [binCount[i]+binS[i] for i in range(len(binCount))]

        totalBinCount = dict(zip(binSizes, [[binS[i], binCount[i]]for i in range(len(binCount))]))
        chr_param["meta"] = {"genomeSize": genomeLen, "totalBinCount": totalBinCount}

        return chr_param


    def get_chromosome_from_array_index(self, chr_param, resolution, index):
        """
        Identify the chromosome based on either the x or y coordinate in the
        array.
        """
        for chr_id in chr_param.keys():
            if chr_id == "meta":
                continue
            #print chr_param[chr_id]["bins"], type(chr_param[chr_id]["bins"])
            chr_end = chr_param[chr_id]["bins"][int(resolution)][2]
            chr_start = chr_param[chr_id]["bins"][int(resolution)][1]
            if index >= chr_start and index <= chr_end:
                return chr_id
