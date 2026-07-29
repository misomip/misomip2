import zenodo_get as zg
from pathlib import Path
import re

#===================================================================================
def download_from_zenodo(zenodo_ID,experiment,output_dir='.'):
    """
    Download a full Zenodo repository for a given model and experiment

    Args:
       * zenodo_ID [string]: number in the Zenodo url, e.g., '21511729' 
       * experiment [string]: MISOMIP2 experiment and region, e.g. 'OceanA-hind', 'OceanW-hind'
       * output_dir [string]: directory in which downloaded netcdf files are stored

    Returns: 
       Download files, makes some prints and returns institute, model

    Exemple:
       institute, model = download_from_zenodo(zenodo_ID='21511729',experiment='OceanA-hind',output_dir='./DATA') 

    """

    print('===================================================')
 
    # create output directory if it does not exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # get file list:
    filelist=[]
    zg.download(record_or_doi=zenodo_ID,output_dir='.',file_glob="*_"+experiment+"_*.nc",md5=True)
    with open("md5sums.txt", "r", encoding="ascii") as file:
      for line in file:
          filelist.append(line.strip().split('  ', 1)[1])
    Nfiles=len(filelist)
    print(Nfiles,' files for experiment ',experiment,' in Zenodo record ',zenodo_ID)
    Path("md5sums.txt").unlink() # delete file

    # get corresponding institute and model names:
    parts = re.split('_', filelist[0])
    ins = parts[1]
    mod = parts[2]
    if not parts[3] == 'a' :
       print('ERROR: NEED TO ADAPT THE SCRIPTS FROM MORE THAN ONE REALISATION PER MODEL AND INSTITUTE')
       exit()
    if not parts[4] in ['OceanA-hind', 'OceanW-hind', 'OceanA-Pgeom', 'OceanW-Pgeom', 'OceanA-Fgeom', 'OceanW-Fgeom', 'OceanA-warm', 'OceanW-warm', \
                        'IceOceanA-hind', 'IceOceanW-hind', 'IceOceanA-ctrl', 'IceOceanW-ctrl', 'IceOceanA-warm', 'IceOceanW-warm']:
       print('ERROR: EXPERIMENT DOES NOT EXIST IN MISOMIP2')
       exit()

    print(mod,' ',ins,'   Zenodo record: ',zenodo_ID)
    print('   ',Nfiles,' files for experiment ',experiment)

    # download all files (tries several times as can failed with poor connections):
    for ktry in range(10):
       dodnl = False
       for file in filelist:
          file_path = Path(output_dir+'/'+file) 
          if not file_path.is_file():
             dodnl = True
             print('   Missing file: ',file,'  >>>>>> starting or completing file download...')
       if dodnl:
          # Only download files that are not already there (use start_fresh=True to dowload from scratch):
          zg.download(record_or_doi=zenodo_ID,output_dir=output_dir,file_glob="*_"+experiment+"_*.nc",start_fresh=False,continue_on_error=True)
       else:
          #nico: I did not find a simple way to check that all downloads are complete, so next step to have more chance that this is the case...
          print('   All files are there, running a last download to complete potential unseccessful dowloads')
          zg.download(record_or_doi=zenodo_ID,output_dir=output_dir,file_glob="*_"+experiment+"_*.nc",start_fresh=False,continue_on_error=True)
          break

    return ins, mod


#===================================================================================
def download_MIPkit_from_zenodo(region,output_dir='.'):
    """
    Download a full Zenodo repository for a given MIPkit dataset

    Args:
      region [string]: 'A' (Amundsen) or 'W' (Weddell)
      output_dir [string]: directory in which we download the MIPkit data

    Returns:
      Download files and returns 0 if several files in the dataset

    Example:
      status = download_MIPkit_from_zenodo(region='A',output_dir='MIPkit/MIPkit-A')

    """

    if region == 'A':
       doi = '10.5281/zenodo.10062355' # latest version
    elif region == 'W':
       doi = '10.5281/zenodo.8316180' # latest version

    # create output directory if it does not exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # get file list:
    filelist=[]
    zg.download(record_or_doi=doi,output_dir='.',file_glob="Oce*.nc",md5=True)
    with open("md5sums.txt", "r", encoding="ascii") as file:
      for line in file:
          filelist.append(line.strip().split('  ', 1)[1])
    Nfiles=len(filelist)
    print(Nfiles,' files for MIPkit-'+region)
    Path("md5sums.txt").unlink() # delete file

    # download all files (tries several times as can failed with poor connections):
    for ktry in range(10):
       dodnl = False
       for file in filelist:
          file_path = Path(output_dir+'/'+file)
          if not file_path.is_file():
             dodnl = True
             print('   Missing file: ',file,'  >>>>>> starting or completing file download...')
       if dodnl:
          # Only download files that are not already there (use start_fresh=True to dowload from scratch):
          zg.download(record_or_doi=doi,output_dir=output_dir,file_glob="Oce*.nc",start_fresh=False,continue_on_error=True)
       else:
          #nico: I did not find a simple way to check that all downloads are complete, so next step to have more chance that this is the case...
          print('   All files are there, running a last download to complete potential unseccessful dowloads')
          zg.download(record_or_doi=doi,output_dir=output_dir,file_glob="Oce*.nc",start_fresh=False,continue_on_error=True)
          break

    if Nfiles > 0:
      return 0
    else:
      return -1
