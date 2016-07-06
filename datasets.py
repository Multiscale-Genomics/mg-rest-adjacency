import json
import numpy as np

"""
Functions to handle the datasets JSON file and the calaulation of the bin cut
off points. The reading of the HDF5 files should be left to the hdf5_reader.py
functions.
"""

class datasets:
    
    def __init__(self):
        # Initialise the chr_param as an empty index when the service is started
        self.datasets = json.loads(open('datasets.json').read())
        self.chr_param = {}
        self.load_datasets()
    
    def load_datasets(self):
        # Set the base bin sizes and initial offsets
        binSizes = [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
        binCount = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        for taxon_id in self.datasets["taxon_id"]:
            for accession_id in self.datasets["taxon_id"][taxon_id]["accession"]:
                genomeLen = 0
                for i in xrange(len(self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"])):
                    c = self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"][i]
                    
                    genomeLen += c[1]
            
                    # Calculate the number of bins for a chromosome and then join with
                    # the offset values for the start in the array
                    binS = [int(np.ceil(c[1]/float(y))) for y in binSizes]
                    binC = dict(zip(binSizes, [[binS[j], binCount[j]] for j in range(len(binCount))]))
                    if self.chr_param.has_key(accession_id):
                        self.chr_param[accession_id][c[0]] = {'size': [c[1], genomeLen], 'bins': binC}
                    else:
                        self.chr_param[accession_id] = {c[0]: {'size':[c[1], genomeLen], 'bins': binC}}
                    
                    # Calculate the new offset values.
                    binCount = [binCount[i]+binS[i] for i in range(len(binCount))]
                self.chr_param[accession_id]["meta"] = {"genomeSize": genomeLen}
    
    def getTaxon(self):
        return self.datasets["taxon_id"].keys()
    
    def getAccessions(self, taxon_id):
        return self.datasets["taxon_id"][taxon_id]["accession"].keys()
    
    def getDatasets(self, taxon_id, accession_id):
        return self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["datasets"]
    
    def getChromosomes(self, taxon_id, accession_id, resolution):
        return self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"]
    
    def getChromosome(self, accession_id, resolution, chr_id):
        return {
            "bins": self.chr_param[accession_id][chr_id]["bins"][resolution][0],
            "size": self.chr_param[accession_id][chr_id]["size"]
        }
    
    """
    Identify the chromosome based on either the x or y coordinate in the array.
    """
    def get_chromosome_from_array_index(self, accession_id, resolution, index):
        # Check that the genome indexes are loaded
        if self.chr_param.has_key(accession_id) == False:
            self.load_chromosome_sizes(accession_id)
        
        for chr_id in self.chr_param[accession_id]:
            if chr_id == "meta":
                continue
            chr_end = self.chr_param[accession_id][chr_id]["bins"][int(resolution)][1] + self.chr_param[accession_id][chr_id]["bins"][int(resolution)][0]
            chr_start = self.chr_param[accession_id][chr_id]["bins"][int(resolution)][1]
            if index > chr_start and index < chr_end:
                return chr_id
    
    def getOffset(self, accession_id, resolution, chr_id):
        return self.chr_param[accession_id][chr_id]["bins"][int(resolution)][1]
    
    def getBinCount(self, accession_id, resolution, chr_id):
        return self.chr_param[accession_id][chr_id]["bins"][int(resolution)][0]
