from check_utils import check_files
from check_utils import check_dims_vars_attrs
from download_utils import download_from_zenodo
from pathlib import Path

# list the ID of all Zenodo submissions (number in the url, https://zenodo.org/records/xxxxxxxx):
zenodo_ID_list = [
# hincast:
#'21511729', # hind IGE-CNRS-UGA      NEMO4.0
#'21626519', # hind UCLA-UMD          MITgcm
#'21728621', # hind UTAS              ROMS
#'23339559', # hind ANU               MOM6
'23318739', # hind DartmouthCollege  MITgcm66j
'24167365', # hind BAS               MITgcm_67d (Weddell)
'24167347'  # hind BAS               MITgcm_67s (Amundsen)
# warm scenario:
            # warm UCLA-UMD          MITgcm
]

#'19568656', # hind AWI          FESOM2 (Amundsen)
#'21452420'  # hind AWI          FESOM2 (Weddell)

# directory in which downloaded netcdf files are stored:
dir_nc = './DATA'

# MISOMIP2 experiments: 'OceanA-hind', 'OceanW-hind', ...
experiment_list = ['OceanA-hind', 'OceanW-hind', 'OceanA-warm', 'OceanW-warm']

for zID in zenodo_ID_list:

   if zID == '19568656':
      experiment_list2 = ['OceanA-hind']
   elif zID == '21452420':
      experiment_list2 = ['OceanW-hind']
   else:
      experiment_list2 = experiment_list

   for exp in experiment_list2:

      # Download files:
      inst, mod = download_from_zenodo(zenodo_ID=zID,experiment=exp,output_dir=dir_nc)

      # Correct BAS MITgcm naming error
      if mod[0:9] == 'MITgcm-67': 
         print('*******************************************')
         print('*  WARNING: renaming BAS MITgcm files...  *')
         print('*******************************************')
         folder = Path(dir_nc)
         for file_path in folder.glob("*BAS_MITgcm_67*"):
            if file_path.is_file():
               new_stem = file_path.stem.replace("MITgcm_", "MITgcm-")
               new_name = f"{new_stem}{file_path.suffix}"
               file_path.rename(folder / new_name)

      if ( not inst == 'None' ):

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
