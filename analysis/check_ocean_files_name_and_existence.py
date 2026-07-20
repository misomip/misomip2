import numpy as np
import xarray as xr
import glob

#==================================================
# User's specification for checking all outputs

institute='IGE-CNRS-UGA'
model='NEMO4.0'

dir_model='DATA/NEMO4.0-IGE-CNRS-UGA_a'

#==================================================
# 1- check file name and existence:

nfile_miss = 0

for region in ['A', 'W']:

   Nmoor = 8
   for kmoor in range(Nmoor):
      IDmoor = str(kmoor+1)
      file_mod=glob.glob(dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_Ocean'+region+'-hind_*.nc')
      if not file_mod:
         print('Missing file for region',region,' and Mooring ',IDmoor)
         print('  Expected: '+dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_Ocean'+region+'-hind_<period>.nc')
         nfile_miss = nfile_miss + 1

   if region == 'A':
      Nsec = 2
   else:
      Nsec = 4
   for ksec in range(Nsec):
      IDsec = str(ksec+1)
      file_mod=glob.glob(dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_Ocean'+region+'-hind_*.nc')
      if not file_mod:
         print('Missing file for region',region,' and Section ',IDsec)
         print('  Expected: '+dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_Ocean'+region+'-hind_<period>.nc')
         nfile_miss = nfile_miss + 1

   file_mod=glob.glob(dir_model+'/Oce3d_'+institute+'_'+model+'_Ocean'+region+'-hind_*.nc')
   if not file_mod:
      print('Missing 3d ocean file for region',region)
      print('  Expected: '+dir_model+'/Oce3d_'+institute+'_'+model+'_Ocean'+region+'-hind_<period>.nc')
      nfile_miss = nfile_miss + 1

print(' ')
print(nfile_miss,' files are missing')
