# Multi-model Analysis 
Contains scripts to download and check MISOMIP2 files on the common grids and formats, and to plot multi-model diagnostics.

To use these scripts in another directory, you may need to put something like this in your scripts:
```
import sys
sys.path.append("/Users/whoever/whereever") # or update and export PYTHONPATH
import misomip2.analysis as ma
```

### Scripts to check your ocean files before uploading onto Zenodo

To check that all files are present with correct naming:
```
vi check_ocean_files_name_and_existence.py # or other editor ; adapt institute, model name, directory 
python check_ocean_files_name_and_existence.py
```

To check that every file has the right dimensions, variables and attributes:
```
vi check_ocean_dims_vars_attrs.py # or any similar editor ; adapt institute, model name, directory
python check_ocean_dims_vars_attrs.py
```

### Download and check all the Zenodo repositories for a given list of MISOMIP2 ocean experiments

As a prerequisite, install zenodo-get:
```
pip install zenodo-get
```

Then, update the list of existing Zenodo repositories in this script and execute it:
```
vi get_and_check_ocean_all.py # or any similar editor
python get_and_check_ocean_all.py
```

### Download the ocean datasets in MIPkit-A and MIPkit-W

This will download the latest version of both datasets:
```
python get_all_latest_MIPkit_ocean_data.py
```

### Making ocean plots for all models

The following scripts put plots in the ./figure directory.

To compare available model simulations to mean observed mooring data:
```
python plot_all_moorings_time_avg.py
ls -al figures/Moorings_*_mean_profiles_MISOMIP2.*
```

To compare available model simulations to all observed vertical sections:
```
python plot_thetao_all_sections.py
ls -al figures/Section*_MISOMIP2.*
```

To download and interpolate masks of individual ice shelf drainage basins:
```
python get_and_interpolate_ice_shelf_mask.py
ls -al DATA/Mask_Ice_Shelves_*.nc 
```


**TO BE COMPLETED**
