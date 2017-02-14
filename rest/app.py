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

from flask import Flask, make_response, request
from flask_restful import Api, Resource

from hdf5_reader import hdf5

app = Flask(__name__)
#app.config['DEBUG'] = False

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


class GetEndPoints(Resource):
    """
    Class to handle the http requests for returning information about the end
    points
    """
    
    def get(self):
        return {
            '_links': {
                '_self': request.base_url,
                '_details': request.url_root + 'mug/api/adjacency/details',
                '_getInteractions': request.url_root + 'mug/api/adjacency/getInteractions',
                '_getValue': request.url_root + 'mug/api/adjacency/getValue',
                '_ping': request.url_root + 'mug/api/adjacency/ping',
                '_parent': request.url_root + 'mug/api'
            }
        }


class GetDetails(Resource):
    """
    Class to handle the http requests for the size of the chromosome, the
    number of bins and available resolutions
    """
    
    def usage(self, error_message, status_code, parameters = {}):
        usage = {
                    '_links' : {
                        '_self' : request.base_url,
                        '_parent': request.url_root + 'mug/api/adjacency'
                    },
                    'parameters' : {
                        'user_id' : ['User ID', 'str', 'REQUIRED'],
                        'file_id' : ['File ID', 'str', 'REQUIRED'],
                    }
                }
        message = {
                      'usage' : usage,
                      'status_code' : status_code
                  }

        if len(parameters) > 0:
            message['provided_parameters'] = parameters
        
        if error_message != None:
            message['error'] = error_message

        return message
    
    def get(self):
        user_id = request.args.get('user_id')
        file_id = request.args.get('file_id')
        
        params = [user_id, file_id]

        # Display the parameters available
        if sum([x is None for x in params]) == len(params):
            return self.usage(None, 200)
        
        # ERROR - one of the required parameters is NoneType
        if sum([x is not None for x in params]) != len(params):
            return self.usage('MissingParameters', 400, {'user_id' : user_id, 'file_id' : file_id}), 400
        
        request_path = request.path
        rp = request_path.split("/")
        
        h5 = hdf5()
        x = h5.get_details(user_id, file_id)
        chr_param = x["chr_param"]
        
        
        return {
            '_links': {
                '_self': request.base_url,
                '_parent': request.url_root + 'mug/api/adjacency'
            },
            'chromosomes' : [{'chromosome' : c[0], 'length' : c[1]} for c in x["chromosomes"]],
            'resolutions' : x['resolutions']
        }


class GetInteractions(Resource):
    """
    Class to handle the http requests for retrieving ranges of interactions from
    a given dataset
    """
    
    def usage(self, error_message, status_code, parameters = {}):
        usage = {
                    '_links' : {
                        '_self' : request.base_url,
                        '_parent': request.url_root + 'mug/api/adjacency'
                    },
                    'parameters' : {
                        'user_id' : ['User ID', 'str', 'REQUIRED'],
                        'file_id' : ['File ID', 'str', 'REQUIRED'],
                        'chr'     : ['Chromosome ID', 'str', 'REQUIRED'],
                        'start'   : ['Chromosome start position', 'int', 'REQUIRED'],
                        'end'     : ['Chromosome end position', 'int', 'REQUIRED'],
                        'res'     : ['Resolution', 'int', 'REQUIRED'],
                        'limit_chr' : ['Limit interactions to interacting with a specific chromosome', 'str', 'OPTIONAL']
                    }
                }
        message = {
                      'usage' : usage,
                      'status_code' : status_code
                  }
        
        if len(parameters) > 0:
            message['provided_parameters'] = parameters
        
        if error_message != None:
            message['error'] = error_message

        return message
    
    def get(self):
        user_id = request.args.get('user_id')
        file_id = request.args.get('file_id')
        chr_id = request.args.get('chr')
        start = request.args.get('start')
        end = request.args.get('end')
        resolution = request.args.get('res')
        limit_chr = request.args.get('limit_chr')
        limit_start = request.args.get('limit_start')
        limit_end = request.args.get('limit_end')
        no_links   = request.args.get('no_links')
        
        params = [user_id, file_id, chr_id, start, end, resolution]

        # Display the parameters available
        if sum([x is None for x in params]) == len(params):
            return self.usage(None, 200)
        
        # ERROR - one of the required parameters is NoneType
        if sum([x is not None for x in params]) != len(params):
            return self.usage('MissingParameters', 400, {'user_id' : user_id, 'file_id' : file_id, 'chr' : chr_id, 'start' : start, 'end' : end, 'res' : resolution, 'limit_chr' : limit_chr}), 400
        
        h5 = hdf5()
        details = h5.get_details(user_id, file_id)
        
        # ERROR - the requested resolution is not available
        if resolution not in details["resolutions"]:
            return self.usage('Resolution Not Available', 400, {'user_id' : user_id, 'file_id' : file_id, 'chr' : chr_id, 'start' : start, 'end' : end, 'res' : resolution, 'limit_chr' : limit_chr})
        
        try:
            start = int(start)
            end = int(end)
            resolution = int(resolution)
        except Exception as e:
            # ERROR - one of the parameters is not of integer type
            return self.usage('IncorrectParameterType', 400, {'user_id' : user_id, 'file_id' : file_id, 'chr' : chr_id, 'start' : start, 'end' : end, 'res' : resolution, 'limit_chr' : limit_chr})
        
        if limit_start is not None or limit_end is not None:
            if limit_chr is None:
                return self.usage('MissingParameters', 400, {'user_id' : user_id, 'file_id' : file_id, 'chr' : chr_id, 'start' : start, 'end' : end, 'res' : resolution, 'limit_chr' : limit_chr, 'limit_start' : limit_start, 'limit_end' : limit_end})
            try:
                limit_start = int(limit_start)
                limit_end = int(limit_end)
            except Exception as e:
                # ERROR - one of the parameters is not of integer type
                return self.usage('IncorrectParameterType', 400, {'user_id' : user_id, 'file_id' : file_id, 'chr' : chr_id, 'start' : start, 'end' : end, 'res' : resolution, 'limit_chr' : limit_chr, 'limit_start' : limit_start, 'limit_end' : limit_end})
        
        request_path = request.path
        rp = request_path.split("/")
        value_url = request.url_root + 'mug/api/adjacency/getValue'
        
        x = h5.get_range(user_id, file_id, resolution, chr_id, start, end, limit_chr, limit_start, limit_end, value_url, no_links)
        #app.logger.warn(x["log"])
        
        return {
            '_links': {
                '_self': request.url,
                '_parent': request.url_root + 'mug/api/adjacency'
            },
            'resolution': resolution,
            'chr': chr_id,
            'start': start,
            'end': end,
            'limit_region': limit_region,
            'limit_chr': limit_chr,
            'interaction_count': len(x["results"]),
            'values': x["results"],
            'log': x["log"]
        }

class GetValue(Resource):
    """
    Class to handle the http requests for retrieving a single value from a given
    dataset
    """
    
    def usage(self, error_message, status_code, parameters = {}):
        usage = {
                    '_links' : {
                        '_self' : request.base_url,
                        '_parent': request.url_root + 'mug/api/adjacency'
                    },
                    'parameters' : {
                        'user_id' : ['User ID', 'str', 'REQUIRED'],
                        'file_id' : ['File ID', 'str', 'REQUIRED'],
                        'res'     : ['Resolution', 'int', 'REQUIRED'],
                        'bin_i'   : ['Position i', 'int', 'REQUIRED'],
                        'bin_j'   : ['Position j', 'int', 'REQUIRED'],
                    }
                }
        message = {
                      'usage' : usage,
                      'status_code' : status_code
                  }
        
        if len(parameters) > 0:
            message['provided_parameters'] = parameters
        
        if error_message != None:
            message['error'] = error_message

        return message
    
    def get(self):
        user_id = request.args.get('user_id')
        file_id = request.args.get('file_id')
        resolution = request.args.get('res')
        bin_i = request.args.get('pos_x')
        bin_j = request.args.get('pos_y')
        
        params = [user_id, file_id, resolution, bin_i, bin_j]

        # Display the parameters available
        if sum([x is None for x in params]) == len(params):
            return self.usage(None, 200)
        
        # ERROR - one of the required parameters is NoneType
        if sum([x is not None for x in params]) != len(params):
            return self.usage('MissingParameters', 400, {'user_id' : user_id, 'file_id' : file_id, 'resolution' : resolution, 'pos_x' : bin_i, 'pos_y' : bin_j}), 400

        try:
            bin_i = int(bin_i)
            bin_j = int(bin_j)
            resolution = int(resolution)
        except Exception as e:
            # ERROR - one of the parameters is not of integer type
            return self.usage('IncorrectParameterType', 400, {'user_id' : user_id, 'file_id' : file_id, 'resolution' : resolution, 'pos_x' : bin_i, 'pos_y' : bin_j}), 400
        
        h5 = hdf5()
        meta_data = h5.get_details(user_id, file_id)
        value = h5.get_value(user_id, file_id, resolution, bin_i, bin_j)
        
        chrA_id = h5.get_chromosome_from_array_index(meta_data["chr_param"], resolution, bin_i)
        chrB_id = h5.get_chromosome_from_array_index(meta_data["chr_param"], resolution, bin_j)
        
        request_path = request.path
        rp = request_path.split("/")
        
        return {
            '_links': {
              '_self': request.url_root + 'mug/api/adjacency/getValue?user_id=' + str(user_id) + "&file_id=" + str(file_id) + "&res=" + str(resolution) + "&pos_x=" + str(bin_i) + "&pos_y=" + str(bin_j)
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
        from . import release
        res = {
            "status":  "ready",
            "version": release.__version__,
            "author":  release.__author__,
            "license": release.__license__,
            "name":    release.__rest_name__,
            "description": release.__description__,
            "_links" : {
                '_self' : request.url_root + 'mug/api/adjacency/ping',
                '_parent' : request.url_root + 'mug/api/adjacency'
            }
        }
        return res

"""
Define the URIs and their matching methods
"""
#   List the available end points for this service
api.add_resource(GetEndPoints, "/mug/api/adjacency", endpoint='adjacency_root')

#   Show the size of the chromosome, the number of bins and available resolutions
#   Parameters:
#    - file_id - (string)
#    - user_id - (string)
api.add_resource(GetDetails, "/mug/api/adjacency/details", endpoint='meta_data')

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
api.add_resource(GetInteractions, "/mug/api/adjacency/getInteractions", endpoint='values')

#   Get a specific edge value for an interaction
#   Parameters:
#    - res     - resolution (int)
#    - pos_x   - (int)
#    - pox_y   - (int)
#    - file_id - (string)
#    - user_id - (string)
api.add_resource(GetValue, "/mug/api/adjacency/getValue", endpoint="value")

#   Service ping
api.add_resource(ping, "/mug/api/adjacency/ping", endpoint='adjacency-ping')


"""
Initialise the server
"""
if __name__ == "__main__":
    app.run()
