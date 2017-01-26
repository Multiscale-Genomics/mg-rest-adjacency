# mg-rest-adjacency

Microservice RESTful API for the querying of Adjacency data stored in HDF5 files that have been generated using the code from the mg-storage-hdf5 / mg-process-fastq scripts

# Requirements
- Python 2.7+
- pyenv
- pyenv virtualenv
- pip
- Python Modules
  - h5py
  - NumPy
  - Flask
  - Flask-Restful
  - json
  - pytest
  - Waitress

# Installation
Cloneing from GitHub:
```
git clone https://github.com/Multiscale-Genomics/mg-rest-adjacency.git
```
To get this to be picked up by pip if part of a webserver then:
```
pip install --editable .
```
This should install the required packages listed in the `setup.py` script.


Installation via pip:
```
pip install git+https://github.com/Multiscale-Genomics/mg-rest-adjacency.git
```

# Setting up a server
```
git clone https://github.com/Multiscale-Genomics/mg-rest-adjacency.git

cd mg-rest-adjacency
pyenv virtualenv 2.7.12 mg-rest-adjacency
pyenv activate mg-rest-service
pip install git+https://github.com/Multiscale-Genomics/mg-dm-api.git
pip install -e .
pip deactivate
```
Starting the service:
```
nohup ${PATH_2_PYENV}/versions/2.7.12/envs/mg-rest-adjacency/bin/waitress-serve --listen=127.0.0.1:5002 rest.app:app &
```

# RESTful API
## List end points
```
wget http://127.0.0.1:5002/api/adjacency/getInteractions
```

## List details about a file
```
wget http://127.0.0.1:5002/api/adjacency/getInteractions/getDetails?user_id=<string:user_id>&file_if=<string:file_id>
```

## Get interactions from chromosome range
```
wget http://127.0.0.1:5002/api/adjacency/getInteractions?user_id=<string:user_id>&file_id=<string:file_id>&res=<int:resolution>&chr=<string:chr_id>&start=<int:start>&end=<int:end>
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
wget http://127.0.0.1:5002/api/adjacency/getInteractions?user_id=<string:user_id>&file_id=<string:file_id>&res=<int:resolution>&chr=<string:chr_id>&start=<int:start>&end=<int:end>&limit_region=intra
wget http://127.0.0.1:5002/api/adjacency/getInteractions?user_id=<string:user_id>&file_id=<string:file_id>&res=<int:resolution>&chr=<string:chr_id>&start=<int:start>&end=<int:end>&limit_chr=<string:chrB_id>
```

### Interactions in TSV format
By modifying the header to request `application/tsv` the following request:

```
wget -S -q --header "Accept: application/tsv" http://127.0.0.1:5002/api/adjacency/getInteractions?user_id=<string:user_id>&file_id=<string:file_id>&res=<int:resolution>&chr=<string:chr_id>&start=<int:start>&end=<int:end> -O test.out.tsv``

```

will return the interactions where columns represent:

1. Chromosome 1
2. Starting position for chromosome 1
3. Chromosome 2
4. Starting position for chromosome 2
5. Value

## Get individual value
```
wget http://127.0.0.1:5002/api/adjacency/getValue?user_id=<string:user_id>&file_id=<string:file_id>&res=<int:resolution>&pox_x=<int:pos_x>&pos_y=<int:pos_y>
```

# Testing
Test scripts are located in the `test/` directory. Run `py.test` to from this directory to ensure that the API is working correctly.

The scripts require a valid hdf5 file generated using the scripts from mg-storage-hdf5 and a matching datasets.json file located in the `rest/` directory

