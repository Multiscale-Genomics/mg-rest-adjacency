import h5py
import numpy as np
import json

class hdf5:
    
    def __init__(self):
        # Initialise the chr_param as an empty index when the service is started
        self.chr_param = {}
    
    """
    Load the chromosomes for a given genomes along with the potential bins and 
    offsets for each bin.
    """
    def load_chromosome_sizes(self, genome):
        # Set the base bin sizes and initial offsets
        binSizes = [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
        binCount = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        # Load the chromosome sizes for the whole genome
        chrSizes = open(genome + ".size", "r")
        genomeLen = 0
        for line in chrSizes:
            line = line.rstrip()
            line = line.split("\t")
            
            c = line[0]
            l = int(line[1])
            
            genomeLen += l
            
            # Calculate the number of bins for a chromosome and then join with
            # the offset values for teh start in the array
            binS = [int(np.ceil(l/float(y))) for y in binSizes]
            binC = dict(zip(binSizes, [[binS[i], binCount[i]]for i in range(len(binCount))]))
            if self.chr_param.has_key(genome):
                self.chr_param[genome][c] = {'size': [l, genomeLen], 'bins': binC}
            else:
                self.chr_param[genome] = {c: {'size':[l, genomeLen], 'bins': binC}}
            
            # Calculate the new offset values.
            binCount = [binCount[i]+binS[i] for i in range(len(binCount))]
        chrSizes.close()
        self.chr_param[genome]["meta"] = {"genomeSize": genomeLen}
    
    def get_chromosome_sizes(self, genome):
        # Check that the genome indexes are loaded
        if self.chr_param.has_key(genome) == False:
             self.load_chromosome_sizes(genome)
        return self.chr_param
    
    """
    Identify the chromosome based on either the x or y coordinate in the array.
    """
    def get_chromosome_from_array_index(self, genome, resolution, index):
        # Check that the genome indexes are loaded
        if self.chr_param.has_key(genome) == False:
            self.load_chromosome_sizes(genome)
        
        for chr_id in self.chr_param[genome]:
            if chr_id == "meta":
                continue
            chr_end = self.chr_param[genome][chr_id]["bins"][int(resolution)][1] + self.chr_param[genome][chr_id]["bins"][int(resolution)][0]
            chr_start = self.chr_param[genome][chr_id]["bins"][int(resolution)][1]
            if index > chr_start and index < chr_end:
                return chr_id
    
    def get_genomes(self):
        datasets = json.loads(open('datasets.json').read())
        out = []
        for i in datasets.keys():
            out.append({"genome_id": i})
        return out
        
    def get_datasets(self, genome_id):
        datasets = json.loads(open('datasets.json').read())
        out = []
        for i in datasets[genome_id]:
            out.append({"dataset": i, "link": "<a href='/" + i + "'>"})
        return out
        
    def get_resolutions(self, genome_id, dataset):
        f = h5py.File(dataset + '.hdf5', "r")
        return {"genome_id": genome_id, "dataset": dataset, "resolutions": f.keys()}
    
    """
    Get the interactions that happen within a defined region on a specific
    chromosome. Returns inter and intra interactions with the defined region.
    """
    def get_range(self, dataset, resolution, genome, chr_id, start, end, limit_region, limit_chr):
        # Check that the genome indexes are loaded
        if self.chr_param.has_key(genome) == False:
          self.load_chromosome_sizes(genome)
        
        # Defines columns to get extracted from the array
        x = int(np.floor(float(start)/float(resolution)))
        y = int(np.ceil(float(end)/float(resolution)))
        
        # xy_offset for the chromosome in the super array
        xy_offset = self.chr_param[genome][chr_id]["bins"][int(resolution)][1]
        
        # Open the hdf5 file
        f = h5py.File(dataset + '.hdf5', "r")
        dset = f[resolution]
        if limit_chr != None:
            startB = self.chr_param[genome][limit_chr]["bins"][int(resolution)][1]
            endB = startB + self.chr_param[genome][limit_chr]["bins"][int(resolution)][0]
            
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
                  y_chr = self.get_chromosome_from_array_index(genome, int(resolution), j)
                  
                  r = {"chrA": chr_id, "startA": x_start, "chrB": y_chr, "startB": y_start, "value": int(result[i,j])}
                  
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
        
