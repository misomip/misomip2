import numpy as np
import xarray as xr
import glob
import sys

#==================================================
# User's specification for checking all outputs

institute='IGE-CNRS-UGA'
model='NEMO3.6'

dir_model='./' #NEMO4.0-IGE-CNRS-UGA_a'

#==================================================

for region in ['A', 'W']:

   if ( region == 'A' ):
      Nlev_Moor = 116
   elif ( region == 'W' ):
      Nlev_Moor = 116

   #-------------------------------
   # Moorings :
   Nmoor = 8
   for kmoor in range(Nmoor):
      print('-----------------------------------------------')
      print('-----------------------------------------------')
      IDmoor = str(kmoor+1)
      file_mod=glob.glob(dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')[0]
      print(file_mod)
      ds = xr.open_dataset(file_mod)
      if not ( ( "lev" in ds.dims ) & ( "time" in ds.dims )):
         print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'time']")
      elif not ( ds.lev.size == Nlev_Moor ):
         print("   Wrong size for dimension 'lev', should be "+str(Nlev_Moor))
      for vv in [ "thetao", "so", "levof" ]:
         if not ( vv in ds.data_vars ):
            print('   Missing variable: '+vv)
         else:
            for aa in [ "units", "long_name", "cell_methods" ]:
               if not ( aa in ds[vv].attrs ):
                  print("   Variable '",vv,"' ->  Missing attribute : ",aa)

   #-------------------------------
   # Sections :
   if region == 'A':
      Nsec = 2
   else:
      Nsec = 4
   for ksec in range(Nsec):
      print('-----------------------------------------------')
      print('-----------------------------------------------')
      IDsec = str(ksec+1)
      file_mod=glob.glob(dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')[0]
      print(file_mod)
      ds = xr.open_dataset(file_mod)
      if not ( ( "x" in ds.dims ) & ( "lev" in ds.dims ) & ( "time" in ds.dims )):
          print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'x', 'time']")
      for vv in [ "thetao", "so", "levof" ]:
         if not ( vv in ds.data_vars ):
            print('   Missing variable: '+vv)
         else:
            for aa in [ "units", "long_name", "cell_methods" ]:
               if not ( aa in ds[vv].attrs ):
                  print("   Variable '",vv,"' ->  Missing attribute : ",aa)

   #-------------------------------
   # 3d Oce :
   if region == 'A':
      mlon=211
      mlat=501
   else:
      mlon=271
      mlat=251
   print('-----------------------------------------------')
   print('-----------------------------------------------')
   file_mod=glob.glob(dir_model+'/Oce3d_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')[0]
   #file_mod=glob.glob(dir_model+'/Oce3d_'+model+'-'+institute+'_?_Ocean'+region+'-hind_*.nc')[0]
   print(file_mod)
   ds = xr.open_dataset(file_mod)
   if not ( ( "lon" in ds.dims ) & ( "lat" in ds.dims ) & ( "lev" in ds.dims ) & ( "time" in ds.dims )):
      print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'lat', 'lon', 'time']")
   else:
      if not ( ds.lev.size == 12 ):
         print("   Wrong size for dimension 'lev', should be 12")
      if not ( ds.lon.size == mlon ):
         print("   Wrong size for dimension 'lon', should be ",mlon)
      if not ( ds.lat.size == mlat ):
         print("   Wrong size for dimension 'lat', should be ",mlat)
   for vv2d in [ "sftflf", "sftof", "deptho", "depflf" ]:
      if not ( vv2d in ds.data_vars ):
         print('   Missing variable: '+vv2d)
      else:
         if not ( np.size(ds[vv2d].shape) == 2 ):
            print('   Wrong number of dimensions for variable ',vv2d,'  -> should have 2 dimensions')
         for aa in [ "units", "long_name", "cell_methods" ]:
            if not ( aa in ds[vv2d].attrs ):
               print("   Variable '",vv2d,"' ->  Missing attribute : ",aa)
   for vv3d in [ "levof", "tob", "sob", "tauuo", "tauvo", "msftbarot", "zos", "wfoat", "flandice", "fsitherm", \
                 "wfocorr", "hfs", "libmassbffl", "dydrflf", "thdrflf", "hadrflf", "siconc", "sivol", "siu", "siv" ]:
      if not ( vv3d in ds.data_vars ): 
         print('   Missing variable: '+vv3d)
      else:
         if not ( np.size(ds[vv3d].shape) == 3 ):
            print('   Wrong number of dimensions for variable ',vv3d,'  -> should have 3 dimensions') 
         for aa in [ "units", "long_name", "cell_methods" ]:
            if not ( aa in ds[vv3d].attrs ):
               print("   Variable '",vv3d,"' ->  Missing attribute : ",aa)
   for vv4d in [ "thetao", "so", "uo", "vo" ]:
      if not ( vv4d in ds.data_vars ): 
         print('   Missing variable: '+vv4d)
      else:
         if not ( np.size(ds[vv4d].shape) == 4 ):
            print('   Wrong number of dimensions for variable ',vv4d,'  -> should have 4 dimensions')
         for aa in [ "units", "long_name", "cell_methods" ]:
            if not ( aa in ds[vv4d].attrs ):
               print("   Variable '",vv4d,"' ->  Missing attribute : ",aa)

