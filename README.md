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

# Testing
Test scripts are located in the `test/` directory. Run `pytest` to from this directory to ensure that the API is working correctly.

The scripts require a valid hdf5 file generated using the scripts from mg-storage-hdf5 and a matching datasets.json file located in the `rest/` directory

