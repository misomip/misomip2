from check_utils import check_files
from check_utils import check_dims_vars_attrs
from download_utils import download_from_zenodo

# list the ID of all Zenodo submissions (number in the url, https://zenodo.org/records/xxxxxxxx):
zenodo_ID_list = [
'21511729', # IGE-CNRS-UGA NEMO4.0
'21626519'  # UCLA-UMD     MITgcm
]

# directory in which downloaded netcdf files are stored:
dir_nc = './DATA'

# MISOMIP2 experiments: 'OceanA-hind', 'OceanW-hind', ...
experiment_list = ['OceanA-hind', 'OceanW-hind']

for zID in zenodo_ID_list:

   for exp in experiment_list:

      # Download files:
      inst, mod = download_from_zenodo(zenodo_ID=zID,experiment=exp,output_dir=dir_nc)

      if ( ( exp[0:6] == 'OceanA' ) | ( exp[0:9] == 'IceOceanA' ) ):
         reg='A'
      elif ( ( exp[0:6] == 'OceanW' ) | ( exp[0:9] == 'IceOceanW' ) ):
         reg='W'
      else:
         reg='X'
 
      # Check all file name and existence:
      status = check_files(institute=inst,model=mod,dir_model=dir_nc,region=reg)

      # Check dimensions, variable names and attributes
      if ( status == 0 ):
         status = check_dims_vars_attrs(institute=inst,model=mod,dir_model=dir_nc,region=reg)
