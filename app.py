#import logging
#from logging.handlers import RotatingFileHandler

"""
TODO:
* Work on including children and parent links
   * Need to find out how to get the current URL then munge it to get the parent then conctenate it for the relevant children
* Work out the best way to document the REST API
   * Can this be made language agnostic?
   * Needs to be simple and flexible
* Generate README.md and prepare for loading to GitHub
"""

from flask import Flask, jsonify, request
from flask_restful import Api, Resource

from hdf5_reader import hdf5

app = Flask(__name__)
api = Api(app)

class GetValue(Resource):
    def get(self, genome_id, chr_id, dataset, resolution, bin_id):
        return jsonify(
            {
                'genome': genome_id,
                'chromosome': chr_id,
                'dataset': dataset,
                'resolution': resolution,
                'bin_id': bin_id
            }
        )

class GetDatasets(Resource):
    def get(self, genome_id):
        h5 = hdf5()
        chr_param = h5.get_datasets(genome_id)
        return jsonify(chr_param)

class GetResolutions(Resource):
    def get(self, genome_id, dataset):
        h5 = hdf5()
        chr_param = h5.get_resolutions(genome_id, dataset)
        return jsonify(chr_param)

class GetChromosomes(Resource):
    def get(self, genome_id, dataset, resolution):
        h5 = hdf5()
        chr_param = h5.get_chromosome_sizes(genome_id)
        return jsonify(chr_param[genome_id].keys())

class GetChromosome(Resource):
    def get(self, genome_id, dataset, resolution, chr_id):
        h5 = hdf5()
        chr_param = h5.get_chromosome_sizes(genome_id)
        return jsonify(chr_param[genome_id][chr_id])

class GetInteractions(Resource):
    def get(self, genome, dataset, resolution, chr_id, start, end):
        limit_region = request.args.get('limit_region')
        limit_chr = request.args.get('limit_chr')
        
        h5 = hdf5()
        x = h5.get_range(dataset, resolution, genome, chr_id, start, end, limit_region, limit_chr)
        #app.logger.warn(x["log"])
        return jsonify(
            {
                'dataset': dataset,
                'resolution': resolution,
                'genome': genome,
                'chr_id': chr_id,
                'start': start,
                'end': end,
                'limit_region': limit_region,
                'limit_chr': limit_chr,
                'result': x["results"]
            }
        )

# /rest/v0.0/getChromosomes/<string:genome_id>
api.add_resource(GetChromosomes, "/rest/v0.0/getChromosomes/<string:genome_id>", endpoint="chromosomes")

# /rest/v0.0/getInteractions
#   List the available datasets with links
api.add_resource(GetDatasets, "/rest/v0.0/getInteractions/<string:genome_id>", endpoint='datasets')

# /rest/v0.0/getInteractions
#   List the resolutions available with links
api.add_resource(GetResolutions, "/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>", endpoint='resolutions')

# /rest/v0.0/getInteractions
#   List the Chromosomes and their sizes and number of bins for the given resolution with links
api.add_resource(GetChromosomes, "/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>", endpoint='bins')

# /rest/v0.0/getInteractions
#   Show the size of the chromosome, the number of bins and a link to a minimal set
api.add_resource(GetChromosome, "/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>", endpoint='sizes')

# /rest/v0.0/getInteractions
api.add_resource(GetInteractions, "/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:start>/<int:end>", endpoint='values')

# /rest/v0.0/getValue/<string:genome_id>/<string:dataset>/<string:chr_id>/<int:res_lvl>/<int:bin_id>
#   Get a specific edge value for an interaction
api.add_resource(GetValue, "/rest/v0.0/getValue/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:bin_id>", endpoint="value")

if __name__ == "__main__":
    #handler = RotatingFileHandler('foo.log', maxBytes=10000, backupCount=1)
    #handler.setLevel(logging.DEBUG)
    #app.logger.addHandler(handler)
    app.run()
