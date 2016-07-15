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
    
    
    def get_range(self, ds, dataset, resolution, accession_id, chr_id, start, end, limit_region, limit_chr, value_url = '/rest/v0.0/getValue/9606'):
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
        logText.append({"coord": {"x0": (x+xy_offset), "x1": (y+xy_offset)}, "r_index": len(r_index), "param": {"start": start, "x": x, "end": end, "y": y, "xy_offset": xy_offset, "resolution": resolution, "accession": accession_id, "chr_id": chr_id, 'chr_param': ds.getChr_param()}})
        
        for i in r_index:
            x_start = (i[0]+x+xy_offset)*10000
            y_chr = ds.get_chromosome_from_array_index(accession_id, int(resolution), i[1])
            chrB = ds.getChromosome(accession_id, int(resolution), y_chr)
            y_start = (i[1]*int(resolution))-(chrB["size"][1]-chrB["size"][0])
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
