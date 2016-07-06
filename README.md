# mg-rest-hdf5

Microservice RESTful API for the querying of HDF5 files that have been generated using the code from the mg-storage-hdf5 scripts

# Requirements
- Python 2.7+
- Python Modules
  - h5py
  - NumPy
  - Flask
  - Flask-Restful
  - json
- virtualenv
- pip

# Installation
## Initialise
```
virtualenv env
source env/bin/activate

cd rest
```

## Require data files
The dataset files for each of the genomes:

- datasets.json
```
{
    "taxon_id": {
        "9606": {
            "accession": {
                "GCA_000001405.14": {
                    "datasets": ["rao2014"],
                    "chromosomes": [
                        ["1", 249250621],
                        ["2", 243199373],
                        ["3", 198022430],
                        ["4", 191154276],
                        ["5", 180915260],
                        ...
                    ]
                }
            }
        }
    }
}
```

Matching HDF5 files listed in the datasets.json file
```
<author><year>.hdf5
```

## Starting the Service
```
virtualenv env
source env/bin/activate

cd rest
python app.py
```

Place this in the boot scripts to get intialised as a service.

# RESTful API
## List taxon IDs
```
wget http://<host>/rest/v0.0/getInteractions/
```

## List accessions
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>
```

## List datasets
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>
```

## List available resolutions
List the avaiable resolutions the are loaded for a given dataset
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>
```

## List chromosomes
List the chromosomes at the given resolution
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>
```

## Size of the chromosome
Show the size of the chromosome, the number of bins and a link to a minimal set
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chr_id>
```

## Get interactions from chromosome range
```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:start>/<int:end>
```
### Optional arguments:
- limit_region
  - Permitted values are:
    - "all" (default): Any interactions
    - "intra": Interactions only within the same chromosome
    - "inter": Interactions only with other chromosomes
- limit_chr
  - Interactions with a specified chromosome

```
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:start>/<int:end>?limit_region=intra
wget http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chrA_id>/<int:start>/<int:end>?limit_chr=<string:chrB_id>
```

### Interactions in TSV format
By modifying the header to request `application/tsv` the following request:

`wget -S -q --header "Accept: application/tsv" http://<host>/rest/v0.0/getInteractions/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:chrA_id>/<int:start>/<int:end> -O test.out.tsv``

```

will return the interactions where columns represent:
1. Chromosome 1
2. Starting position for chromosome 1
3. Chromosome 2
4. Starting position for chromosome 2
5. Value

## Get individual value
```
wget http://<host>/rest/v0.0/getValue/<string:taxon_id>/<string:accession_id>/<string:dataset>/<int:resolution>/<string:bin_i>/<int:bin_j>
```


