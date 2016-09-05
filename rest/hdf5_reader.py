"""
Copyright 2016 EMBL-European Bioinformatics Institute

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

import pkg_resources, os
import h5py
import numpy as np
import json

class hdf5:
    """
    Class related to handling the functions for interacting directly with the
    HDF5 files. All required information should be passed to this class.
    """
    
    def get_resolutions(self, accession_id, dataset):
        """
        Return a list of the available resolutions in a given HDF5 file
        """
        resource_package = __name__
        resource_path = os.path.join(os.path.dirname(__file__), dataset + '.hdf5')
        f = h5py.File(resource_path, "r")
        return {"accession_id": accession_id, "dataset": dataset, "resolutions": f.keys()}
    
    
    def get_range(self, ds, dataset, resolution, accession_id, chr_id, start, end, limit_region=None, limit_chr=None, value_url = '/rest/v0.0/getValue/9606'):
        """
        Get the interactions that happen within a defined region on a specific
        chromosome. Returns inter and intra interactions with the defined region.
        """
        
        # Defines columns to get extracted from the array
        x = int(np.floor(float(start)/float(resolution)))
        y = int(np.ceil(float(end)/float(resolution)))
        
        # xy_offset for the chromosome in the super array
        xy_offset = ds.getOffset(accession_id, resolution, chr_id)
        
        # Open the hdf5 file
        resource_package = __name__
        resource_path = os.path.join(os.path.dirname(__file__), dataset + '.hdf5')
        f = h5py.File(resource_path, "r")
        dset = f[str(resolution)]
        if limit_chr != None:
            startB = ds.getOffset(accession_id, resolution, limit_chr)
            endB = startB + ds.getBinCount(accession_id, resolution, limit_chr)
            
            result = dset[(x+xy_offset):(y+xy_offset),startB:endB]
        else:
            result = dset[(x+xy_offset):(y+xy_offset),:]
        f.close()
        
        # Iterate through slice and extract results greater than zero
        results = []
        rShape = result.shape
        logText = []
        
        r_index = np.transpose(np.nonzero(result))
        logText.append({"coord": {"x0": (x+xy_offset), "x1": (y+xy_offset)}, "r_index": len(r_index), "param": {"start": start, "x": x, "end": end, "y": y, "xy_offset": xy_offset, "resolution": resolution, "accession": accession_id, "chr_id": chr_id}, 'chr_param': ds.getChr_param()})
        
        for i in r_index:
            x_start = ((i[0]+x)*int(resolution))
            y_chr = ds.get_chromosome_from_array_index(accession_id, int(resolution), i[1])
            y_start = (i[1]-ds.getOffset(accession_id, resolution, y_chr))*int(resolution)
            r = {"chrA": chr_id, "startA": x_start, "chrB": y_chr, "startB": y_start, "value": int(result[i[0],i[1]]), '_links': {'self': value_url + "/" + str(accession_id) + "/" + str(dataset) + "/" + str(resolution) + "/" + str(i[0]+x+xy_offset) + "/" + str(i[1])}}
            results.append(r)
        
        return {"log": logText, "results": results}
    
    def get_value(self, dataset, resolution, bin_i, bin_j):
        """
        Get a specific value for a given dataset, resoltuoin
        """
        resource_package = __name__
        resource_path = os.path.join(os.path.dirname(__file__), dataset + '.hdf5')
        f = h5py.File(resource_path, "r")
        dset = f[str(resolution)]
        value = dset[int(bin_i), int(bin_j)]
        f.close()
        return value
