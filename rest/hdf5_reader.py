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

import pkg_resources, os, json, h5py
import numpy as np

from dmp import dmp

class hdf5:
    """
    Class related to handling the functions for interacting directly with the
    HDF5 files. All required information should be passed to this class.
    """
    
    def get_details(self, user_id, file_id):
        """
        Return a list of the available resolutions in a given HDF5 file
        """
        
        if user_id == 'test':
            resource_package = __name__
            resource_path = os.path.join(os.path.dirname(__file__), 'rao2014.hdf5')
            f = h5py.File(resource_path, "r")
        else:
            cnf_loc=os.path.dirname(os.path.abspath(__file__)) + '/mongodb.cnf'
            da = dmp(cnf_loc, test_data)
            file_obj = da.get_file_by_id(user_id, file_id)
            f = h5py.File(file_obj["file_path"], "r")
        
        resolutions = f.keys()
        dset = f[str(resolutions[0])]
        chr_param = _calculate_chr_param(resolutions, dset.attrs["chromosomes"])
        return {
            "chromosomes" : dset.attrs["chromosomes"],
            "chr_param"   : chr_param,
            "resolutions" : resolutions}
    
    
    def get_range(self, user_id, file_id, resolution, chr_id, start, end, limit_region=None, limit_chr=None, value_url = '/rest/v0.0/getValue/9606'):
        """
        Get the interactions that happen within a defined region on a specific
        chromosome. Returns inter and intra interactions with the defined region.
        """
        
        # Defines columns to get extracted from the array
        x = int(np.floor(float(start)/float(resolution)))
        y = int(np.ceil(float(end)/float(resolution)))
        
        if user_id == 'test':
            resource_package = __name__
            resource_path = os.path.join(os.path.dirname(__file__), 'rao2014.hdf5')
            f = h5py.File(resource_path, "r")
        else:
            # Open the hdf5 file
            cnf_loc=os.path.dirname(os.path.abspath(__file__)) + '/mongodb.cnf'
            da = dmp(cnf_loc)
            file_obj = da.get_file_by_id(user_id, file_id)
            f = h5py.File(file_obj["file_path"], "r")
        
        # Get meta data from HDF5 file
        resolutions = map(int, f.keys())
        dset = f[str(resolutions[0])]
        chr_param = self._calculate_chr_param(resolutions, dset.attrs["chromosomes"])
        
        # xy_offset for the chromosome in the super array
        xy_offset = chr_param[chr_id]["bins"][resolution][1]
        
        dset = f[str(resolution)]
        
        if limit_chr != None:
            startB = chr_param[limit_chr]["bins"][resolution][1]
            endB = startB + chr_param[limit_chr]["bins"][resolution][0]
            
            result = dset[(x+xy_offset):(y+xy_offset),startB:endB]
        else:
            result = dset[(x+xy_offset):(y+xy_offset),:]
        f.close()
        
        # Iterate through slice and extract results greater than zero
        results = []
        rShape = result.shape
        logText = []
        
        r_index = np.transpose(np.nonzero(result))
        logText.append({"coord": {"x0": (x+xy_offset), "x1": (y+xy_offset)}, "r_index": len(r_index), "param": {"start": start, "x": x, "end": end, "y": y, "xy_offset": xy_offset, "resolution": resolution, "chr_id": chr_id}, 'chr_param': chr_param})
        
        for i in r_index:
            x_start = ((i[0]+x)*int(resolution))
            y_chr = self.get_chromosome_from_array_index(chr_param, int(resolution), i[1])
            y_start = (i[1]-chr_param[y_chr]["bins"][resolution][1])*int(resolution)
            r = {"chrA": chr_id, "startA": x_start, "chrB": y_chr, "startB": y_start, "value": int(result[i[0],i[1]]), '_links': {'self': value_url + "/getInteractions?user_id=" + str(user_id) + "&file_id=" + str(file_id) + "&res=" + str(resolution) + "&pos_x=" + str(i[0]+x+xy_offset) + "&pos_y=" + str(i[1])}}
            results.append(r)
        
        return {"log": logText, "results": results}
    
    
    def get_value(self, user_id, file_id, resolution, bin_i, bin_j):
        """
        Get a specific value for a given dataset, resoltuoin
        """
        if user_id == 'test':
            resource_package = __name__
            resource_path = os.path.join(os.path.dirname(__file__), 'rao2014.hdf5')
            f = h5py.File(resource_path, "r")
        else:
            cnf_loc=os.path.dirname(os.path.abspath(__file__)) + '/mongodb.cnf'
            da = dmp(cnf_loc, test_data)
            file_obj = da.get_file_by_id(user_id, file_id)
            f = h5py.File(file_obj["file_path"], "r")
        
        dset = f[str(resolution)]
        value = dset[int(bin_i), int(bin_j)]
        f.close()
        
        return value
    
    def _calculate_chr_param(self, binSizes, chromosomes):
        """
        Load the self.chr_param object with the required information about the
        start and stop positions of the chromosomes and bins.
        """
        
        chr_param = {}
        
        genomeLen = 0
        binCount = [0]*len(binSizes)
        for i in xrange(len(chromosomes)):
            c = chromosomes[i]
            
            genomeLen += int(c[1])
    
            # Calculate the number of bins for a chromosome and then join with
            # the offset values for the start in the array
            binS = [int(np.ceil(int(c[1])/float(y))) for y in binSizes]
            binC = dict(zip(binSizes, [[binS[j], binCount[j], binS[j]+binCount[j]] for j in range(len(binCount))]))
            
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
            chr_end = chr_param[chr_id]["bins"][int(resolution)][2]
            chr_start = chr_param[chr_id]["bins"][int(resolution)][1]
            if index >= chr_start and index <= chr_end:
                return chr_id
