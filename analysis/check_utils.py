import numpy as np
import xarray as xr
import glob
import sys

#==================================================
def check_files(institute,model,dir_model,region):
    """
    Basic checks for MISOMIP2 ocean output file names.

    Args: 
       * institute [string]: e.g. 'IGE-CNRS-UGA'
       * model [string]: e.g. 'NEMO3.6'
       * dir_model [string]: directory in which are stored the netcdf files, e.g. './'
       * region [string]: 'A' (Amundsen) or 'W' (Weddell)

    Returns: 
       Nb of missing files (0 if all files there)

    Example:
       status = check_files(institute='IGE-CNRS-UGA',model='NEMO3.6','./DATA',region='A')

    """

    nfile_miss = 0
 
    Nmoor = 8
    for kmoor in range(Nmoor):
       IDmoor = str(kmoor+1)
       file_mod=glob.glob(dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_a_Ocean'+region+'-hind_*.nc')
       if not file_mod:
          print('Missing file for region',region,' and Mooring ',IDmoor)
          print('  Expected: '+dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_a_Ocean'+region+'-hind_<period>.nc')
          nfile_miss = nfile_miss + 1
 
    if region == 'A':
       Nsec = 2
    else:
       Nsec = 4
    for ksec in range(Nsec):
       IDsec = str(ksec+1)
       file_mod=glob.glob(dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_a_Ocean'+region+'-hind_*.nc')
       if not file_mod:
          print('Missing file for region',region,' and Section ',IDsec)
          print('  Expected: '+dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_a_Ocean'+region+'-hind_<period>.nc')
          nfile_miss = nfile_miss + 1
 
    file_mod=glob.glob(dir_model+'/Oce3d_'+institute+'_'+model+'_a_Ocean'+region+'-hind_*.nc')
    if not file_mod:
       print('Missing 3d ocean file for region',region)
       print('  Expected: '+dir_model+'/Oce3d_'+institute+'_'+model+'_a_Ocean'+region+'-hind_<period>.nc')
       nfile_miss = nfile_miss + 1
 
    print(' ')
    print(nfile_miss,' files are missing for region',region)
 
    return nfile_miss


#==================================================
def check_dims_vars_attrs(institute,model,dir_model,region):
    """
    Basic checks for MISOMIP2 ocean outputs: dimensions and variable attributes.

    Args: 
       * institute [string]: e.g. 'IGE-CNRS-UGA'
       * model [string]: e.g. 'NEMO3.6'
       * dir_model [string]: directory in which are stored the netcdf files, e.g. '.'
       * region [string]: 'A' (Amundsen) or 'W' (Weddell)

    Returns: 
       Many prints and the number of errors

    Example:
       status = check_dims_vars_attrs(institute='IGE-CNRS-UGA',model='NEMO3.6',dir_model='.',region='W')

    """

    nerr = 0
 
    #-------------------------------
    # Moorings :
    Nmoor = 8
    Nlev_Moor1 = 1151 # mooring output every 1m (as in the article)  
    Nlev_Moor2 = 116  # mooring output every 10m (more convenient, accepted format)
    for kmoor in range(Nmoor):
       IDmoor = str(kmoor+1)
       file_mod=glob.glob(dir_model+'/OceMoor'+IDmoor+'_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')
       if not file_mod:
          print('WARNING: NO FILE FOR MOORING ',IDmoor)
       else:
          file_mod=file_mod[0]
          print(file_mod)
          ds = xr.open_dataset(file_mod)
          if not ( ( "lev" in ds.dims ) & ( "time" in ds.dims )):
             print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'time'] (in any order)")
             nerr = nerr + 1
          elif not ( ( ds.lev.size == Nlev_Moor1 ) | ( ds.lev.size == Nlev_Moor2 ) ):
             print("   Wrong size for dimension 'lev', should be either "+str(Nlev_Moor1)+" or "+str(Nlev_Moor2))
             nerr = nerr + 1
          for vv in [ "thetao", "so", "levof" ]:
             if not ( vv in ds.data_vars ):
                print('   Missing variable: '+vv)
                nerr = nerr + 1
             else:
                for aa in [ "units", "long_name", "cell_methods" ]:
                   if not ( aa in ds[vv].attrs ):
                      print("   Variable '",vv,"' ->  Missing attribute : ",aa)
                      nerr = nerr + 1
 
    #-------------------------------
    # Sections :
    if region == 'A':
       Nsec = 2
    else:
       Nsec = 4
    for ksec in range(Nsec):
       IDsec = str(ksec+1)
       file_mod=glob.glob(dir_model+'/OceSec'+IDsec+'_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')[0]
       print(file_mod)
       ds = xr.open_dataset(file_mod)
       if not ( ( "x" in ds.dims ) & ( "lev" in ds.dims ) & ( "time" in ds.dims )):
          print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'x', 'time'] (in any order)")
          nerr = nerr + 1
       for vv in [ "thetao", "so", "levof" ]:
          if not ( vv in ds.data_vars ):
             print('   Missing variable: '+vv)
             nerr = nerr + 1
          else:
             for aa in [ "units", "long_name", "cell_methods" ]:
                if not ( aa in ds[vv].attrs ):
                   print("   Variable '",vv,"' ->  Missing attribute : ",aa)
                   nerr = nerr + 1
 
    #-------------------------------
    # 3d Oce :
    if region == 'A':
       mlon=501
       mlat=211
    else:
       mlon=271
       mlat=251
    file_mod=glob.glob(dir_model+'/Oce3d_'+institute+'_'+model+'_?_Ocean'+region+'-hind_*.nc')[0]
    print(file_mod)
    ds = xr.open_dataset(file_mod)
    if not ( ( "lon" in ds.dims ) & ( "lat" in ds.dims ) & ( "lev" in ds.dims ) & ( "time" in ds.dims )):
       print('   Wrong dimensions: ',list(ds.dims)," should be ['lev', 'lat', 'lon', 'time'] (in any order)")
       nerr = nerr + 1
    else:
       if not ( ds.lev.size == 12 ):
          print("   Wrong size for dimension 'lev', should be 12")
          nerr = nerr + 1
       if not ( ds.lon.size == mlon ):
          print("   Wrong size for dimension 'lon', should be ",mlon)
          nerr = nerr + 1
       if not ( ds.lat.size == mlat ):
          print("   Wrong size for dimension 'lat', should be ",mlat)
          nerr = nerr + 1
    for vv2d in [ "sftflf", "sftof", "deptho", "depflf" ]:
       if not ( vv2d in ds.data_vars ):
          print('   Missing variable: '+vv2d)
          nerr = nerr + 1
       else:
          if not ( np.size(ds[vv2d].shape) == 2 ):
             print('   Wrong number of dimensions for variable ',vv2d,'  -> should have 2 dimensions')
             nerr = nerr + 1
          for aa in [ "units", "long_name", "cell_methods" ]:
             if not ( aa in ds[vv2d].attrs ):
                print("   Variable '",vv2d,"' ->  Missing attribute : ",aa)
                nerr = nerr + 1
    for vv3d in [ "levof", "tob", "sob", "tauuo", "tauvo", "msftbarot", "zos", "wfoat", "flandice", "fsitherm", \
                  "wfocorr", "hfs", "libmassbffl", "dydrflf", "thdrflf", "hadrflf", "siconc", "sivol", "siu", "siv" ]:
       if not ( vv3d in ds.data_vars ): 
          print('   Missing variable: '+vv3d)
          nerr = nerr + 1
       else:
          if not ( np.size(ds[vv3d].shape) == 3 ):
             print('   Wrong number of dimensions for variable ',vv3d,'  -> should have 3 dimensions') 
             nerr = nerr + 1
          for aa in [ "units", "long_name", "cell_methods" ]:
             if not ( aa in ds[vv3d].attrs ):
                print("   Variable '",vv3d,"' ->  Missing attribute : ",aa)
                nerr = nerr + 1
    for vv4d in [ "thetao", "so", "uo", "vo" ]:
       if not ( vv4d in ds.data_vars ): 
          print('   Missing variable: '+vv4d)
          nerr = nerr + 1
       else:
          if not ( np.size(ds[vv4d].shape) == 4 ):
             print('   Wrong number of dimensions for variable ',vv4d,'  -> should have 4 dimensions')
             nerr = nerr + 1
          for aa in [ "units", "long_name", "cell_methods" ]:
             if not ( aa in ds[vv4d].attrs ):
                print("   Variable '",vv4d,"' ->  Missing attribute : ",aa)
                nerr = nerr + 1

    if ( nerr == 0 ):
       print('[All good !]')
 
    return nerr
