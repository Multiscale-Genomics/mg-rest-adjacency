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
        f = h5py.File(dataset + '.hdf5', "r")
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
        f = h5py.File(dataset + '.hdf5', "r")
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
        for i in xrange(rShape[0]):
            # Convert array location to region start position
            x_start = i*int(resolution)
            
            for j in xrange(rShape[1]):
                if result[i,j] > 0:
                  # Convert array location to region start position
                  y_start = j*int(resolution)
                  
                  # Identify the chromosome from the array location
                  y_chr = ds.get_chromosome_from_array_index(accession_id, int(resolution), j)
                  
                  r = {"chrA": chr_id, "startA": x_start, "chrB": y_chr, "startB": y_start, "value": int(result[i,j]), '_links': {'self': value_url + "/" + str(accession_id) + "/" + str(dataset) + "/" + str(resolution) + "/" + str(i) + "/" + str(j)}}
                  
                  if limit_region != None:
                      if limit_region.encode("utf-8") == "inter":
                          if chr_id != y_chr:
                              results.append(r)
                      elif limit_region.encode("utf-8") == "intra":
                          if chr_id == y_chr:
                              results.append(r)
                      else:
                          results.append(r)
                  else:
                      results.append(r)
        
        return {"log": logText, "results": results}
    
    def get_value(self, dataset, resolution, bin_i, bin_j):
        """
        Get a specific value for a given dataset, resoltuoin
        """
        f = h5py.File(dataset + '.hdf5', "r")
        dset = f[str(resolution)]
        value = dset[int(bin_i), int(bin_j)]
        f.close()
        return value
