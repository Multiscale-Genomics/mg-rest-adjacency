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

from flask import Flask, make_response, request
from flask_restful import Api, Resource

#from datasets import datasets
from dmp import dmp
from hdf5_reader import hdf5

app = Flask(__name__)
api = Api(app)

@api.representation('application/tsv')
def output_tsv(data, code, headers=None):
    """
    TSV representation for interactions
    """
    if request.endpoint == "values":
        outstr = ''
        for v in data["values"]:
            outstr += str(v["chrA"]) + "\t" + str(v["startA"]) + "\t" + str(v["chrB"]) + "\t" + str(v["startB"]) + "\t" + str(v["value"]) + "\n"
        resp = make_response(outstr, code)
        resp.headers.extend(headers or {})
        return resp


class GetDetails(Resource):
    """
    Class to handle the http requests for the size of the chromosome, the
    number of bins and available resolutions
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution, chr_id):
        ds = datasets()
        chr_param = ds.getChromosome(accession_id, resolution, chr_id)
        
        request_path = request.path
        rp = request_path.split("/")
        
        h5 = hdf5()
        x = h5.get_details(user_id, file_id)
        chr_param = x["chr_param"]
        
        return {
            '_links': {
                'self': request.base_url,
                'child': request.base_url + "/0/" + str(resolution),
                'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id) + '/' + str(dataset) + '/' + str(resolution) + "/" + str(chr_id)
            },
            'chromosomes': x["chromosomes"],
            "bins": chr_param["bins"],
            "size": chr_param["size"]
        }

class GetInteractions(Resource):
    """
    Class to handle the http requests for retrieving ranges of interactions from
    a given dataset
    """
    
    def get(self):
        user_id = request.args.get('user_id')
        file_id = request.args.get('file_id')
        chr_id = request.args.get('chr')
        start = int(request.args.get('start'))
        end = int(request.args.get('end'))
        resolution = int(request.args.get('res'))
        limit_region = request.args.get('limit_region')
        limit_chr = request.args.get('limit_chr')
        
        if chr_id not in chromosomes:
            return {
                "error" : "No chr parameter provided"
            }
        
        if resolution == None:
            return {
                "error" : "No res parameter provided"
            }
        
        if start == None or end == None:
            return {
                "error" : "No start and/or end parameter provided"
            }
        
        
        request_path = request.path
        rp = request_path.split("/")
        value_url = request.url_root + 'rest/' + str(rp[2]) + '/getValue/' + str(taxon_id)
        
        #ds = datasets()
        h5 = hdf5()
        x = h5.get_range(user_id, file_id, resolution, accession_id, chr_id, start, end, limit_region, limit_chr, value_url)
        #app.logger.warn(x["log"])
        
        return {
            '_links': {
                'self': request.url,
                'parent': request.url_root + '/rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id) + '/' + str(dataset) + '/' + str(resolution) + "/" + str(chr_id)
            },
            'dataset': dataset,
            'resolution': resolution,
            'genome': accession_id,
            'chr_id': chr_id,
            'start': start,
            'end': end,
            'limit_region': limit_region,
            'limit_chr': limit_chr,
            'interaction_count': len(x["results"]),
            #'page_id': 
            #'pages': 
            #'page_size': 
            'values': x["results"],
            'log': x["log"]
        }

class GetValue(Resource):
    """
    Class to handle the http requests for retrieving a single value from a given
    dataset
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution, bin_i, bin_j):
        h5 = hdf5()
        meta_data = h5.get_details(user_id, file_id)
        value = h5.get_value(user_id, file_id, resolution, bin_i, bin_j)
        
        chrA_id = h5.get_chromosome_from_array_index(meta_data["chr_param"], resolution, bin_i)
        chrB_id = h5.get_chromosome_from_array_index(meta_data["chr_param"], resolution, bin_j)
        
        request_path = request.path
        rp = request_path.split("/")
        
        return {
            '_links': {
              'self': request.url_root + 'rest/' + str(rp[2]) + '/getValue/' + str(taxon_id) + "/" + str(accession_id) + "/" + str(dataset) + "/" + str(resolution) + "/" + str(bin_i) + "/" + str(bin_j)
            },
            'genome': accession_id,
            'chrA': chrA_id,
            'chrB': chrB_id,
            'dataset': dataset,
            'resolution': resolution,
            'bin_i': bin_i,
            'bin_j': bin_j,
            'value': int(value)
        }

class ping(Resource):
    """
    Class to handle the http requests to ping a service
    """
    
    def get(self):
        import release
        res = {
            "status":  "ready",
            "version": release.__version__,
            "author":  release.__author__,
            "license": release.__license__,
            "name":    release.__rest_name__,
            "description": release.__description__,
            "end_points":  {}
        }
        return res

"""
Define the URIs and their matching methods
"""
#   List the available end points for this service
api.add_resource(GetEndPoints, "/api/adjacency", endpoint='adjacency_root')

#   Show the size of the chromosome, the number of bins and available resolutions
#   Parameters:
#    - file_id - (string)
#    - user_id - (string)
api.add_resource(GetDetails, "/api/adjacency/details", endpoint='meta_data')

#   List the interactions for a given region
#   Parameters:
#    - chr     - chromosome (string)
#    - res     - resolution (int)
#    - start   - (int)
#    - end     - (int)
#    - file_id - (string)
#    - user_id - (string)
#   Parameters (optional):
#    - limit_region - 
#    - limit_chr    - Chromosome (string)
api.add_resource(GetInteractions, "/api/adjacency/getInteractions", endpoint='values')

#   Get a specific edge value for an interaction
#   Parameters:
#    - res     - resolution (int)
#    - pos_x   - (int)
#    - pox_y   - (int)
#    - file_id - (string)
#    - user_id - (string)
api.add_resource(GetValue, "/api/adjacency/getValue", endpoint="value")

#   Service ping
api.add_resource(ping, "/api/adjacency/ping", endpoint='adjacency-ping')


"""
Initialise the server
"""
if __name__ == "__main__":
    app.run()
