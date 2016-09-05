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
import json
import numpy as np

class datasets:
    """
    Class for handling all functions related to the interactions with the
    datasets and returning information about those datasets.
    """
    
    def __init__(self):
        """
        Intialisation function for the datasets class
        """
        
        # Load the datasets.json at the beginning
        resource_package = __name__
        resource_path = os.path.join('datasets.json')
        dataset_file = pkg_resources.resource_string(resource_package, resource_path)
        self.datasets = json.loads(dataset_file)
        
        # Initialise the chr_param as an empty index when the service is started
        self.chr_param = {}
        
        # Load the chr_param based on the info in the datasets.json file
        self.load_datasets()
    
    def load_datasets(self):
        """
        Load the self.chr_param object with the required information about the
        start and stop positions of the chromosomes and bins.
        """
        
        # Set the base bin sizes and initial offsets
        binSizes = [1000, 5000, 10000, 25000, 50000, 100000, 250000, 500000, 1000000]
        
        for taxon_id in self.datasets["taxon_id"]:
            for accession_id in self.datasets["taxon_id"][taxon_id]["accession"]:
                genomeLen = 0
                binCount = [0, 0, 0, 0, 0, 0, 0, 0, 0]
                for i in xrange(len(self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"])):
                    c = self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"][i]
                    
                    genomeLen += c[1]
            
                    # Calculate the number of bins for a chromosome and then join with
                    # the offset values for the start in the array
                    binS = [int(np.ceil(c[1]/float(y))) for y in binSizes]
                    binC = dict(zip(binSizes, [[binS[j], binCount[j], binS[j]+binCount[j]] for j in range(len(binCount))]))
                    if self.chr_param.has_key(accession_id):
                        self.chr_param[accession_id][c[0]] = {'size': [c[1], genomeLen], 'bins': binC}
                    else:
                        self.chr_param[accession_id] = {c[0]: {'size':[c[1], genomeLen], 'bins': binC}}
                    
                    # Calculate the new offset values.
                    binCount = [binCount[i]+binS[i] for i in range(len(binCount))]
                
                totalBinCount = dict(zip(binSizes, [[binS[i], binCount[i]]for i in range(len(binCount))]))
                self.chr_param[accession_id]["meta"] = {"genomeSize": genomeLen, "totalBinCount": totalBinCount}
    
    def getTaxon(self):
        """
        Return a list of taxon IDs
        """
        return self.datasets["taxon_id"].keys()
    
    def getAccessions(self, taxon_id):
        """
        Return a list of accession IDs for a given taxon ID
        """
        return self.datasets["taxon_id"][taxon_id]["accession"].keys()
    
    def getDatasets(self, taxon_id, accession_id):
        """
        Return a list of datasets for a given taxon and accession
        """
        return self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["datasets"]
    
    def getChromosomes(self, taxon_id, accession_id, resolution):
        """
        Return a list of chromosomes for a given taxon and accession
        """
        return self.datasets["taxon_id"][taxon_id]["accession"][accession_id]["chromosomes"]
    
    def getChromosome(self, accession_id, resolution, chr_id):
        """
        Return a dict of chromosome properties for a given accession, resolution
        and chromosome.
        """
        return {
            "bins": self.chr_param[accession_id][chr_id]["bins"][resolution][0],
            "bin_offset": self.chr_param[accession_id][chr_id]["bins"][resolution][1],
            "size": self.chr_param[accession_id][chr_id]["size"]
        }
    
    
    def get_chromosome_from_array_index(self, accession_id, resolution, index):
        """
        Identify the chromosome based on either the x or y coordinate in the
        array.
        """
        for chr_id in self.chr_param[accession_id].keys():
            if chr_id == "meta":
                continue
            chr_end = self.chr_param[accession_id][chr_id]["bins"][int(resolution)][2]
            chr_start = self.chr_param[accession_id][chr_id]["bins"][int(resolution)][1]
            if index >= chr_start and index <= chr_end:
                return chr_id
    
    def getOffset(self, accession_id, resolution, chr_id):
        """
        Return an integer for the offset for a given chromosome at a given
        resolution.
        """
        return self.chr_param[accession_id][chr_id]["bins"][int(resolution)][1]
    
    def getBinCount(self, accession_id, resolution, chr_id):
        """
        Return a count of the number of bins for a given dataset and chromosome
        """
        return self.chr_param[accession_id][chr_id]["bins"][int(resolution)][0]
    
    def getTotalBinCount(self, accession_id, resolution):
        """
        Return a count of the total number of bins for a given dataset and resolution
        """
        return self.chr_param[accession_id]["meta"]["totalBinCount"][int(resolution)][1]
    
    def getChr_param(self):
        return self.chr_param
