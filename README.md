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
"<genome_id>": ["<author><year>"]
}
```

Matching genome size file for each genome in the datasets.json file:
- <genome_id>.size
```
1       249250621
2       243199373
3       198022430
4       191154276
5       180915260
...
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
## List datasets
```
wget http://<host>/rest/v0.0/getInteractions/<string:genome_id>
```

## List available resolutions
List the avaiable resolutions the are loaded for a given dataset
```
wget http://<host>/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>
```

## List chromosomes
List the chromosomes at the given resolution
```
wet http://<host>/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>
```

## Size of the chromosome
Show the size of the chromosome, the number of bins and a link to a minimal set
```
wget http://<host>/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>
```

## Get interactions from chromosome range
```
http://<host>/rest/v0.0/getInteractions/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:start>/<int:end>
```

## Get individual value
```
http://<host>/rest/v0.0/getValue/<string:genome_id>/<string:dataset>/<int:resolution>/<string:chr_id>/<int:bin_id>
```


