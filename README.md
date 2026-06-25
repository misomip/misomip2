# misomip2
A python repo to process model outputs to standard MISOMIP2 format and to analyse MISOMIP2 multi-model outputs. This is a repo that is meant to help modellers prepare their model outputs to the standard format required in MISOMIP2 (and, later, to analyse the outputs). It is not plug and play and requires some work for individual modellers to adapt the scripts to their model outputs. Please carefully check the results of these scripts as only basic checks have been done so far. Modellers are also free to use other scripts to match the requirements listed in [De Rydt et al. (2024)](https://doi.org/10.5194/gmd-17-7105-2024). 

Exampls of use are provided below for outputs from NEMO, MITgcm and ROMS.

The scripts are written for python3 and largely based on [xarray](http://xarray.pydata.org) and [scipy.interpolate](https://docs.scipy.org/doc/scipy/reference/interpolate.html).

This repo includes:
* **preproc** : contains functions used to pre-process model outputs, i.e. to interpolate them to the MISOMIP2 grids and write files with standard attributes ([More details here](https://github.com/nicojourdain/misomip2/tree/master/preproc/README.md)).
* **examples** : contains scripts that can be used on provided test cases (e.g. to check that your python environment works) and that can be adapted to your specific output files.
* **analysis** : (TO BE COMPLETED) will contain scripts used to analyse multiple model outputs that were previously interpolated to the common grids ([More details here](https://github.com/nicojourdain/misomip2/tree/master/analysis/README.md)).

### Contributors
* Nicolas C. Jourdain (IGE, CNRS-UGA, Grenoble, France)
* Jan De Rydt (University Northumbria Newcastle, UK)
* Yoshihiro Nakayama (Hokkaido University, Japan)
* Ole Richter (Alfred Wegener Institute, Germany)

All MISOMIP participants who want to contribute to this repo are invited to fork this github repo and make a pull request as [documented here](https://opensource.com/article/19/7/create-pull-request-github).

### Install
If you don't have a python environment yet, you can install anaconda following [this page](https://docs.anaconda.com/anaconda/install/).

If these modules are not installed (check with ```conda list```), install them:
```bash
conda install numpy xarray scipy netcdf4
conda install dask pyproj
conda install -c conda-forge gsw 
# or conda install -c conda-forge/label/gcc7 gsw 
# or conda install -c conda-forge/label/cf201901 gsw 
# or conda install -c conda-forge/label/cf202003 gsw
```

Then, clone the misomip2 repo :

```bash
git clone https://github.com/misomip/misomip2.git
# or git clone git@github.com:misomip/misomip2.git
```

You can update the cloned directory anytime with ```git pull``` executed in that directory. 

<br/><br/>
----------
----------

# Examples of ocean interpolation to the MISOMIP2 grids

We provide a few test cases so that users can check that this package works well in their python environment before adapting it to their model outputs. These test cases are also used by the developers of this package to check that it is working with various types of model. 

For simplicity, we describe the example directly used in the misomip2 repositories, but the interpolation script can be run from anywhere.

```bash
cd misomip2/examples
vi interpolate_to_common_grid_oce.py # or using any other text editor than vi
```
The first thing to change in this file is the directory in which the misomip2 package has been cloned, i.e. change this path:
```python
sys.path.append("/Users/jourdain/MY_PACKAGES")
```
The second thing is to choose the test case you want to run (among NEMO, eORCA025, MITGCM and ROMS).
We recommend starting with 'NEMO_test' which is the smallest in size and fastest to run.
To make it work, you need to download the input files for the tes cases, which are all provided on [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4709851.svg)](https://doi.org/10.5281/zenodo.4709851). You can get the files directly using wget as:
```bash
cd test_cases/oce
wget https://zenodo.org/record/4709851/files/NEMO_test.zip
unzip NEMO_test.zip && rm NEMO_test.zip
```

Other test cases can be downloaded as follows (note the file sizes):
```bash
# MITGCM_test (242M):
wget https://zenodo.org/record/4709851/files/MITGCM_test.zip
unzip MITGCM_test.zip && rm MITGCM_test.zip
# eORCA025_test (7.8G):
wget https://zenodo.org/record/4709851/files/eORCA025_test.zip
unzip eORCA025_test.zip && rm eORCA025_test.zip
# ROMS_test (257M):
wget https://zenodo.org/record/4709851/files/ROMS_test.zip
unzip ROMS_test.zip && rm ROMS_test.zip
```

Then, execute the script as follows using the yaml card of the selected test case:
```bash
python interpolate_to_common_grid_oce.py --config test_cases/oce/NEMO_config.yml
```
The default value is ```test_cases/oce/NEMO_config.yml```

This should create the following files (here for NEMO and section 1 and mooring 2):
* Oce3d\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc
* OceSec1\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc
* OceMoor2\_NEMO3.6-IGE-CNRS-UGA\_a\_Ocean-A1\_201001-201002.nc

On a laptop (16Gb, 2GHz), the test cases took the following durations:
* NEMO\_test : 48s (smaller than MISOMIP2 domain; all variables calculated). 
* MITGCM_test : 3min 30s (smaller than MISOMIP2 domain; no sea-ice or surface fluxes).
* eORCA025_test : 2min 30s (global simulation; most variables calculated).
* ROMS_test : 6min 20s (circum-Antarctic simulation (sigma coordinates); no sea-ice or surface fluxes).

### Adapt to your own ocean configuration

You can copy interpolate_to_common_grid_oce.py and adapt it to your model. You need to adapt at least the yml card to fit your simulation, and section 3a (Files and variables). The yml card define what regions you want to process (*processing*), the simulation description (*simulation*) and the output metadata (*metadata*). In section 3a, you may need to create or modify a load\_oce\_mod\_xxxx.py function similar to the one existing for NEMO, MITgcm and ROMS if your model is not covered yet. If your model is already covered but variable names differ from what is understood by load\_oce\_mod\_xxxx.py, just add options for these variables in load\_oce\_mod\_xxxx.py and [make a pull request](https://opensource.com/article/19/7/create-pull-request-github) to upload it onto the official misomip2 repository (so that you keep this in case of updates).

### Performance of ocean interpolation

The scripts were written in best effort by non-experts in a way that should work for various grids, and suggestions that may boost perfromances are welcome. We nonetheless found that the computing times were acceptable for several example of simulations, the idea being that the higher the original model resolutoin, the shorter the period to be computed in a single call of interpolate\_to\_common\_grid\_oce.py. Here are some example of computing time:

* Processing the "NEMO-AMUXL025" data (231x190x75 points, reduced to 203x107x75 when "loading nemo") on a 2.6 GHz computer required ~10 minutes and ~2.5Gb of memory per year of simulation (tests up to 10 years in a single call of interpolate\_to\_common\_grid\_oce.py show that both walltime and memory increase linearly). 

* Processing the "NEMO-AMUXL12" data (687x567x75 points, reduced to 623x286x75 when "loading nemo") on a 2.6 GHz computer required 1h22 minutes and a bit less than 16Gb of memory per year of simulation.

<br/><br/>
----------
----------
 
# Ice test cases

**To be completed**

