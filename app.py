from flask import Flask, jsonify, request
from flask_restful import Api, Resource

from datasets import datasets
from hdf5_reader import hdf5

app = Flask(__name__)
api = Api(app)


class GetTaxons(Resource):
    """
    Class to handle the http requests for retrieving a list of taxon IDs
    """
    
    def get(self):
        ds = datasets()
        taxons = ds.getTaxon()
        request_path = request.path
        rp = request_path.split("/")
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': [request.base_url + "/" + str(d) for d in taxons]
                },
                'taxons': taxons
            }
        )

class GetAccessions(Resource):
    """
    Class to handle the http requests for retrieving a list of accessions for a 
    given taxonic ID
    """
    
    def get(self, taxon_id):
        ds = datasets()
        accessions = ds.getAccessions(taxon_id)
        request_path = request.path
        rp = request_path.split("/")
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': [request.base_url + "/" + str(d) for d in accessions],
                    'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id)
                },
                'accessions': accessions
            }
        )

class GetDatasets(Resource):
    """
    Class to handle the http requests for retrieving a list of datasets for a
    given accession
    """
    
    def get(self, taxon_id, accession_id):
        ds = datasets()
        dataset = ds.getDatasets(taxon_id, accession_id)
        request_path = request.path
        rp = request_path.split("/")
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': [request.base_url + "/" + str(d) for d in dataset],
                    'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id)
                },
                'datasets': dataset
            }
        )

class GetResolutions(Resource):
    """
    Class to handle the http requests for retrieving the available resolutions
    that have been loaded in a dataset
    """
    
    def get(self, taxon_id, accession_id, dataset):
        h5 = hdf5()
        chr_param = h5.get_resolutions(accession_id, dataset)
        request_path = request.path
        rp = request_path.split("/")
        children = []
        for c in chr_param["resolutions"]:
            children.append(request.base_url + "/" + c)
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': children,
                    'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id) + '/' + str(dataset)
                },
                'dataset': dataset,
                'accession_id': accession_id,
                'resolutions': chr_param['resolutions'],
            }
        )

class GetChromosomes(Resource):
    """
    Class to handle the http requests for retrieving the list of chromosomes for
    a given accession
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution):
        ds = datasets()
        chr_param = ds.getChromosomes(taxon_id, accession_id)
        request_path = request.path
        rp = request_path.split("/")
        children = []
        chromosomes = []
        for c in chr_param:
          children.append(request.base_url + "/" + c[0])
          chromosomes.append(c[0])
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': children,
                    'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id) + '/' + str(dataset) + '/' + str(resolution)
                },
                'chromosomes': chromosomes
            }
        )

class GetChromosome(Resource):
    """
    Class to handle the http requests for retrieving chromosome information for
    a given chromosome, include teh number of bins at a given resoltuoin and the
    size of the chromosome
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution, chr_id):
        ds = datasets()
        chr_param = ds.getChromosome(accession_id, resolution, chr_id)
        
        request_path = request.path
        rp = request_path.split("/")
        return jsonify(
            {
                '_links': {
                    'self': request.base_url,
                    'child': request.base_url + "/0/" + str(resolution),
                    'parent': request.url_root + 'rest/' + str(rp[2]) + '/' + str(rp[3]) + '/' + str(taxon_id) + '/' + str(accession_id) + '/' + str(dataset) + '/' + str(resolution) + "/" + str(chr_id)
                },
                'chromosome': chr_id,
                "bins": chr_param["bins"],
                "size": chr_param["size"]
            }
        )

class GetInteractions(Resource):
    """
    Class to handle the http requests for retrieving ranges of interactions from
    a given dataset
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution, chr_id, start, end):
        limit_region = request.args.get('limit_region')
        limit_chr = request.args.get('limit_chr')
        
        request_path = request.path
        rp = request_path.split("/")
        value_url = request.url_root + 'rest/' + str(rp[2]) + '/getValue/' + str(taxon_id)
        
        ds = datasets()
        h5 = hdf5()
        x = h5.get_range(ds, dataset, resolution, accession_id, chr_id, start, end, limit_region, limit_chr, value_url)
        #app.logger.warn(x["log"])
        
        return jsonify(
            {
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
                'values': x["results"]
            }
        )

class GetValue(Resource):
    """
    Class to handle the http requests for retrieving a single value from a given
    dataset
    """
    
    def get(self, taxon_id, accession_id, dataset, resolution, bin_i, bin_j):
        h5 = hdf5()
        value = h5.get_value(dataset, resolution, bin_i, bin_j)
        
        ds = datasets()
        chrA_id = ds.get_chromosome_from_array_index(accession_id, resolution, bin_i)
        chrB_id = ds.get_chromosome_from_array_index(accession_id, resolution, bin_j)
        
        request_path = request.path
        rp = request_path.split("/")
        
        return jsonify(
            {
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
        )

"""
Define the URIs and their matching methods
"""
#   List the available species for which there are datasets available
api.add_resource(GetTaxons, "/rest/v0.0/getInteractions", endpoint='taxons')

#   List the available assemblies for a given species with links
api.add_resource(GetAccessions, "/rest/v0.0/getInteractions/<string:taxon_id>", endpoint='accessions')

#   List the available datasets for a given genome with links
api.add_resource(GetDatasets, "/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>", endpoint='datasets')

#   List the resolutions available with links
api.add_resource(GetResolutions, "/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>", endpoint='resolutions')

#   List the Chromosomes and their sizes and number of bins for the given resolution with links
api.add_resource(GetChromosomes, "/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>", endpoint='bins')

#   Show the size of the chromosome, the number of bins and a link to a minimal set
api.add_resource(GetChromosome, "/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chr_id>", endpoint='sizes')

#   List the interactions for a given region
api.add_resource(GetInteractions, "/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:start>/<int:end>", endpoint='values')

#   Get a specific edge value for an interaction
api.add_resource(GetValue, "/rest/v0.0/getValue/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:bin_i>/<int:bin_j>", endpoint="value")


"""
Initialise the server
"""
if __name__ == "__main__":
    app.run()
