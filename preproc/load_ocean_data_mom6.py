# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
# 2026-10 : Modified for MOM6 using some code from ROMS version
#====================================================================================

import numpy as np
import xarray as xr
import gsw
from .def_grids import grid_bounds_oce
from .def_attrs import rename_dimensions
from datetime import datetime
import os, psutil

#====================================================================================
def load_oce_mod_mom6(file_mesh_mask='ocean_static.nc',\
                      file_bathy='ocean_static.nc',\
                      files_gridT='dummy',\
                      files_gridU='dummy',\
                      files_gridV='dummy',\
                      files_SBC='dummy',\
                      files_ice='dummy',\
                      files_is='dummy',\
                      files_BSF='dummy',\
                      rho0=1035.0, teos10=False, region='Weddell', parallel=False, chunks = {'xh':1080,'yh':360,'xq':1080, 'yq':360,'yT':360,'xT':1080} ):
   """ Read MOM6 outputs and define an xarray dataset containing 
       all variables required in MISOMIP2. It is adapted to
       C-grid and standard 
       mesh/mask variables.
       WARNING: Regrids to flat z levels, but assumes grid is static and takes vertical coordinate
       levels from first timestep of provided model output

       Input:

         file_mesh_mask: mesh mask file [default='mesh_mask.nc']

         file_bathy: bathymetry file [optional, look for data in file_mesh_mask if not provided]

         files_gridT: file or list of files containing the temperature and salinity variables
      
         files_gridU: file or list of files containing the x-velocity and related variables

         files_gridV: file or list of files containing the y-velocity and related variables

         files_ice: file or list of files containing the sea-ice variables 

         files_SBC: file or list of files containing the surface fluxes variables 

         files_is: file or list of files containing the ice shelf thickness

         files_BSF: file or list of files containing umo_2d (vertically integrated zonal mass transport) for calculating
         barotropic streamfunction calculated at U-points [look for data in files_gridU if not provided]

         rho0 corresponds to reference density RHO_0 in MOM6 i.e. reference volumic mass [kg m-3]

         teos10=False -> assumes the outputs are in potential temperature & practical salinity (EOS80)
               =True  -> assumes the outputs are in CT and AS and convert to PT and PS

       region = 'Amundsen' or 'Weddell' (default)

       parallel: If True, the open and preprocess steps of this function will be performed in parallel

       Output:
          xarray dataset of coordinates ("time", "z", "sxy") (sxy= one-dimensionalized horizontal space)

       Example:
        #Assumes output stored in output000, output001 etc folders with access-om3 naming conventions
        data_dir='/path/to/dir'
        startoutput = 120
        endoutput = 121
        
        f_mesh   = [data_dir+'/output'+np.int64(startoutput).astype('str').zfill(3)+'/access-om3.mom6.static.nc']
        f_bathy  = [data_dir+'/output'+np.int64(startoutput).astype('str').zfill(3)+'/access-om3.mom6.static.nc'] 
        fs_gridT = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.3d.temp.z.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]+ \
                  [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.tob.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                  [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.3d.e.z.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                   [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.3d.salt.z.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]+ \
                  [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.sob.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]
        fs_gridU = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.tauuo.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                  [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.3d.uo.z.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]
        fs_gridV = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.tauvo.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                  [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.3d.vo.z.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]
        fs_ice   = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.cice.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] 
        fs_is = [data_dir+'/output'+np.int64(startoutput).astype('str').zfill(3)+'/access-om3.mom6.2d.h_shelf.1mon.mean.*.nc' ]
        fs_SBC   = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.zos.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.mass_flux.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.exch_vel_t.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.thermal_driving.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.haline_driving.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.hfds.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.fsitherm.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.wfo.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.ficeberg.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] + \
                [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.salt_flux_added.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)] 
        fs_BSF = [data_dir+'/output'+month.astype('str').zfill(3)+'/access-om3.mom6.2d.umo_2d.1mon.mean.*.nc' for month in np.arange(startoutput,endoutput+1)]
                
        ds = mp.load_oce_mod_mom6( file_mesh_mask=[glob.glob(a)[0] for a in f_mesh], 
                                    file_bathy=[glob.glob(a)[0] for a in f_bathy], 
                                    files_gridT=[glob.glob(a)[0] for a in fs_gridT],
                                    files_gridU=[glob.glob(a)[0] for a in fs_gridU], 
                                    files_gridV=[glob.glob(a)[0] for a in fs_gridV], 
                                    files_SBC=[glob.glob(a)[0] for a in fs_SBC],
                                    files_ice=[glob.glob(a)[0] for a in fs_ice],
                                    files_is = [glob.glob(a)[0] for a in fs_is],
                                    files_BSF=[glob.glob(a)[0] for a in fs_BSF],
                                    rho0=1035.0, teos10=False,parallel=True )
    
   """

   startTime = datetime.now()

   if ( file_bathy == 'None' ):
     file_bathy=file_mesh_mask
   if ( files_BSF == 'None' ):
     files_BSF=files_gridU

   ncM = xr.open_mfdataset(file_mesh_mask,decode_coords=False,decode_timedelta=True,parallel=parallel,chunks=chunks)
   ncM=ncM.isel(time=0)
   ncB = xr.open_mfdataset(file_bathy,decode_coords=False,decode_timedelta=True,parallel=parallel,chunks=chunks)
   ncT = xr.open_mfdataset(files_gridT,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks); ncT=rename_dimensions(ncT)
   ncU = xr.open_mfdataset(files_gridU,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks).drop('xh').drop('h'); ncU=rename_dimensions(ncU)
   ncV = xr.open_mfdataset(files_gridV,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks).drop('yh').drop('h'); ncV=rename_dimensions(ncV)
   ncS = xr.open_mfdataset(files_SBC,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks).drop(['xT','yT']); ncS=rename_dimensions(ncS)
   ncI = xr.open_mfdataset(files_ice,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks); ncI=rename_dimensions(ncI)
   ncIS = xr.open_mfdataset(files_is,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks); ncIS=rename_dimensions(ncIS)
   ncP = xr.open_mfdataset(files_BSF,decode_coords=False,decode_timedelta=True, parallel=parallel,chunks=chunks); ncP=rename_dimensions(ncP)
    
   # fill all nans with a non-nan value (useful to keep points in triangular interpolation, even if they will be masked)
   # from ROMS example

   def fill_float_nan(ds, fill_value=9.969209968386869e36):
    ds = ds.copy()
    for v in ds.data_vars:
        if np.issubdtype(ds[v].dtype, np.floating):
            ds[v] = ds[v].fillna(fill_value)
    return ds

   ncM, ncB, ncT, ncS, ncU, ncV, ncI, ncIS, ncP = [fill_float_nan(ds) for ds in [ncM, ncB, ncT, ncS, ncU, ncV, ncI, ncIS, ncP]]

   mtime = ncT.time.shape[0]

   # Define useful masks on U, V, T grids 
   maskT = rename_dimensions(ncM.wet); maskT = maskT.where(maskT<1e18).fillna(0) # 3d mask (=1 if ocean, =0 elsewhere)
   maskU = rename_dimensions(ncM.wet_u); maskU = maskU.where(maskU<1e18).fillna(0) # 3d mask
   maskV = rename_dimensions(ncM.wet_v); maskV = maskV.where(maskV<1e18).fillna(0) # 3d mask

   [mz, my, mx] = ncT.temp.isel(time = 0).shape

   # Domain mask (ones with a halo of nans), used not to interpolate beyond the original domain:
   halonan = np.ones((my,mx))
   halonan[0,:] = np.nan ; halonan[-1,:] = np.nan
   halonan[:,0] = np.nan ; halonan[:,-1] = np.nan
   DOMMSKT = xr.DataArray( halonan, dims=['y', 'x'] )
   DOMMSKU = xr.DataArray( halonan, dims=['y', 'x'] )
   DOMMSKV = xr.DataArray( halonan, dims=['y', 'x'] )

   # longitude & latitude on U, V, T grids
   lonT = rename_dimensions(ncM.geolon.where(ncM.geolon<1e18)).bfill('y').ffill('y').bfill('x').ffill('x') ; latT = rename_dimensions(ncM.geolat.where(ncM.geolat<1e18)).bfill('x').ffill('x').bfill('y').ffill('y')
   lonU = rename_dimensions(ncM.geolon_u.where(ncM.geolon_u<1e18)).bfill('y').ffill('y').bfill('x').ffill('x') ; latU = rename_dimensions(ncM.geolat_u.where(ncM.geolat_u<1e18)).bfill('x').ffill('x').bfill('y').ffill('y')
   lonV = rename_dimensions(ncM.geolon_v.where(ncM.geolon_v<1e18)).bfill('y').ffill('y').bfill('x').ffill('x') ; latV = rename_dimensions(ncM.geolat_v.where(ncM.geolat_v<1e18)).bfill('x').ffill('x').bfill('y').ffill('y')

   # save original domain boundaries:
   domain_minlat = latT.min().values
   domain_maxlat = latT.max().values
   domain_minlon = lonT.min().values
   domain_maxlon = lonT.max().values

   dxT = rename_dimensions(ncM.dxt) # grid mesh width along x-axis in meters for sea ice grid (T)
   dyT = rename_dimensions(ncM.dyt) # grid mesh width along y-axis in meters for sea ice grid (T)

   dlatTdx = 0.500000000 * ( latT.shift(x=-1) - latT.shift(x=1) )
   dlatUdx = 0.500000000 * ( latU.shift(x=-1) - latU.shift(x=1) )
   dlatVdx = 0.500000000 * ( latV.shift(x=-1) - latV.shift(x=1) )

   dlonTdx = 0.500000000 * ( lonT.shift(x=-1) - lonT.shift(x=1) )
   dlonUdx = 0.500000000 * ( lonU.shift(x=-1) - lonU.shift(x=1) )
   dlonVdx = 0.500000000 * ( lonV.shift(x=-1) - lonV.shift(x=1) )

   # depth of U, V, T grids (neglecting the effects of partial steps in the interpolation)
   depTUV=ncT.z

   # local T, U, V grid rotation angle compared to the (zonal,meridional) direction [rad]
   thetaT = np.arctan2( dlatTdx, dlonTdx*np.cos(latT*np.pi/180.) )
   thetaU = np.arctan2( dlatUdx, dlonUdx*np.cos(latU*np.pi/180.) )
   thetaV = np.arctan2( dlatVdx, dlonVdx*np.cos(latV*np.pi/180.) )
   print('    Minimum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.min().values*180./np.pi)
   print('    Maximum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.max().values*180./np.pi)

   # Ocean fraction at each level:
   # this will be overwritten
   LEVOFT = maskT*100.0
   LEVOFU = maskU*100.0
   LEVOFV = maskV*100.0

   # Ocean fraction at surface:
   # this will be overwritten
   SFTOF = LEVOFT

   # 2d ice-shelf fraction (assumed 0 or 100%):
   SFTFLF = rename_dimensions(ncIS.h_shelf.isel(time = 0))
   SFTFLF = (SFTFLF.where( (SFTFLF > 0) & (SFTFLF < 1e18) )*0+100).fillna(0)

   # Bathymetry (including under ice shelves) [m]
   if ( "deptho" in ncM.data_vars ):
     DEPTHO = rename_dimensions(ncM.deptho)
   else:
     print('    WARNING :   No data found for DEPTHO  -->  filled with NaNs')
     DEPTHO = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['y', 'x'] )
  
   # Depth of ice shelf draft [m]:
   if ( "e" in ncT.data_vars ):
     DEPFLF = -1 * ncT.e.isel(time = 0).isel(z_i = 0)
     # if in open ocean, replace surface value with e = 0
     DEPFLF = DEPFLF.where(SFTFLF >0).fillna(0)
   else:
     print('    WARNING :   No data found for DEPFLF  -->  filled with NaNs')
     DEPFLF = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['y', 'x'] )
 
   # ocean temperature [degC]
   if ( "temp" in ncT.data_vars ):
     TT = ncT.temp
   elif ( "thetao" in ncT.data_vars ):
     TT = ncT.thetao
   else:
     print('    WARNING :   No data found for TT  -->  filled with NaNs')
     TT = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # ocean salinity [1.e-3]
   if ( "salt" in ncT.data_vars ):
     SS = ncT.salt
   elif ( "so" in ncT.data_vars ):
     SS = ncT.so
   else:
     print('    WARNING :   No data found for SS  -->  filled with NaNs')
     SS = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # sea bottom ocean temperature [degC]
   if ( "sbt" in ncT.data_vars ):
     TTB = ncT.sbt
   elif ( "sosbt" in ncT.data_vars ):
     TTB = ncT.sosbt
   elif ( "tob" in ncT.data_vars ):
     TTB = ncT.tob
   else:
     print('    WARNING :   No data found for TTB  -->  filled with NaNs')
     TTB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea bottom ocean salinity [1.e-3]
   if ( "sbs" in ncT.data_vars ):
     SSB = ncT.sbs
   elif ( "sosbs" in ncT.data_vars ):
     SSB = ncT.sosbs
   elif ( "sob" in ncT.data_vars ):
     SSB = ncT.sob
   else:
     print('    WARNING :   No data found for SSB  -->  filled with NaNs')
     SSB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # ocean potential temperature and practical salinity :
   if (teos10):
     TOB = xr.apply_ufunc(gsw.pt_from_CT, SSB, TTB, dask = 'allowed')
     SOB = xr.apply_ufunc(gsw.SP_from_SA, SSB, DEPTHO, lonT, latT, dask = 'allowed')
     THETAO = xr.apply_ufunc(gsw.pt_from_CT, SS, TT, dask = 'allowed')
     SO = xr.apply_ufunc(gsw.SP_from_SA, SS, DEPTHO, lonT, latT, dask = 'allowed')
   else: 
     TOB = TTB
     SOB = SSB
     THETAO = TT
     SO = SS
   
   # ocean x-ward velocity [m s-1]
   if ( "uoce" in ncU.data_vars ):
     UX = ncU.uoce
   elif ( "vozocrtx" in ncU.data_vars ):
     UX = ncU.vozocrtx
   elif ( "uo" in ncU.data_vars ):
     UX = ncU.uo
   else:
     print('    WARNING :   No data found for UX  -->  filled with NaNs')
     UX = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # ocean y-ward velocity [m s-1]
   if ( "voce" in ncV.data_vars ):
     VY = ncV.voce
   elif ( "vomecrty" in ncV.data_vars ):
     VY = ncV.vomecrty
   elif ( "vo" in ncV.data_vars ):
     VY = ncV.vo
   else:
     print('    WARNING :   No data found for VY  -->  filled with NaNs')
     VY = xr.DataArray( np.zeros((mtime,mz,my,mx))*np.nan, dims=['time', 'z', 'y', 'x'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "tauuo" in ncU.data_vars ):
     TAUX = ncU.tauuo
   elif ( "sozotaux" in ncU.data_vars ):
     TAUX = ncU.sozotaux
   else:
     print('    WARNING :   No data found for TAUX  -->  filled with NaNs')
     TAUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "tauvo" in ncV.data_vars ):
     TAUY = ncV.tauvo
   elif ( "sometauy" in ncV.data_vars ):
     TAUY = ncV.sometauy
   else:
     print('    WARNING :   No data found for TAUY  -->  filled with NaNs')
     TAUY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Sea surface height [m]
   if ( "zos" in ncS.data_vars ):
     ZOS = ncS.zos
   elif ( "sossh" in ncS.data_vars ):
     ZOS = ncS.sossh
   elif ( "ssh" in ncS.data_vars ):
     ZOS = ncS.ssh
   elif ( "sossheig" in ncS.data_vars ):
     ZOS = ncS.sossheig
   elif ( "zos" in ncT.data_vars ):
     ZOS = ncT.zos
   elif ( "sossh" in ncT.data_vars ):
     ZOS = ncT.sossh
   elif ( "ssh" in ncT.data_vars ):
     ZOS = ncT.ssh
   elif ( "sossheig" in ncT.data_vars ):
     ZOS = ncT.sossheig
   else:
     print('    WARNING :   No data found for ZOS  -->  filled with NaNs')
     ZOS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # mass barotropic streamfunction
   # see Griffies et al. (2016, section H26): d(psi)/dy=-U (U: x-ward mass transport), d(psi)/dx=V (V: yward mass transport)
   if ( "umo_2d" in ncP.data_vars ):
     MSFTBAROT = ncP.umo_2d.where(ncP.umo_2d<1e18).cumsum(dim='y')
   else:
     print('    WARNING :   No data found for MSFTBAROT  -->  filled with NaNs')
     MSFTBAROT = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Water Mass Flux Into Seawater From Land Ice [kg m−2 s−1, positive downward (ie into sea water)]
   lisf = ("mass_flux" in ncS.data_vars) 
   licb = ("ficeberg" in ncS.data_vars)
   lrnf = ("lrunoff" in ncS.data_vars)
   FLANDICE = xr.DataArray( np.zeros((mtime,my,mx)), dims=['time', 'y', 'x'] )
   if ( lisf or licb or lrnf ) :
     if ( lisf ):
       if ("mass_flux" in ncS.data_vars) & ("areacello" in ncM.data_vars):
         FLANDICE = FLANDICE + (ncS.mass_flux/rename_dimensions(ncM.areacello)*maskT) #convert units
     if ( licb ): 
       FLANDICE = FLANDICE + ncS.ficeberg
     if ( lrnf ):
       FLANDICE = FLANDICE + ncS.lrunoff
   else:
     print('    WARNING :   No data found for FLANDICE  -->  filled with NaNs')
     FLANDICE = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Water Mass Flux Into Seawater Due to Sea Ice Thermodynamics [kg m-2 s-1, positive for actual melting]
   if ( "fsitherm" in ncS.data_vars ):
     FSITHERM = ncS.fsitherm
   else:
     print('    WARNING :   No data found for FSITHERM  -->  filled with NaNs')
     FSITHERM = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # ice shelf melt [kg m-2 s-1, positive for actual melting] :
   if ( "mass_flux" in ncS.data_vars )  & ("areacello" in ncM.data_vars):
     LIBMASSBFFL = (ncS.mass_flux/rename_dimensions(ncM.areacello)*maskT) #convert units
   else:
     print('    WARNING :   No data found for LIBMASSBFFL  -->  filled with NaNs')
     LIBMASSBFFL = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # ice shelf dynamical driving (heat exchange velocity) [m s-1]:
   if ( "exch_vel_t" in ncS.data_vars ):
     DYDRFLF = ncS.exch_vel_t
   else:
     print('    WARNING :   No data found for DYDRFLF  -->  filled with NaNs')
     DYDRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] ) 

   # ice shelf thermal driving [degC]:
   if ( "thermal_driving" in ncS.data_vars ):
     THDRFLF = ncS.thermal_driving
   else:
     print('    WARNING :   No data found for THDRFLF  -->  filled with NaNs')
     THDRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # ice shelf haline driving [0.001]:
   if ( "haline_driving" in ncS.data_vars ):
     HADRFLF = ncS.haline_driving
   elif ( "halined_cav" in ncS.data_vars ):
     HADRFLF = ncS.halined_cav
   else:
     print('    WARNING :   No data found for HADRFLF  -->  filled with NaNs')
     HADRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice concentration [0-100]
   if ( "aice_m" in ncI.data_vars ):
     SICONC = ncI.aice_m*100.0
     #SICONC = SICONC.where( (~np.isnan(SICONC)) & (~np.isinf(SICONC)), 0.e0 )
   else:
     print('    WARNING :   No data found for SICONC  -->  filled with NaNs')
     SICONC = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )   

   # sea-ice volume per area [m]
   if ( "hi_m" in ncI.data_vars ):
     SIVOL = ncI.hi_m
   elif ( "sivol" in ncI.data_vars ):
     SIVOL = ncI.sivol
   else:
     print('    WARNING :   No data found for SIVOL  -->  filled with NaNs')
     SIVOL = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice x-ward velocity [m/s]
   if ( "uvel_m" in ncI.data_vars ):
     SIUX = ncI.uvel_m
   elif ("siu" in ncI.data_vars ):
     SIUX = ncI.siu
   else:
     print('    WARNING :   No data found for SIUX  -->  filled with NaNs')
     SIUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # sea-ice y-ward velocity [m/s]
   if ( "vvel_m" in ncI.data_vars ):
     SIVY = ncI.vvel_m
   elif ("siv" in ncI.data_vars ):
     SIVY = ncI.siv
   else:
     print('    WARNING :   No data found for SIUY  -->  filled with NaNs')
     SIVY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Total heat flux received by the ocean surface (including ice-shelf/ocean interface) [W m-2] 
   # see Griffies et al. (2016, section K4-K5) NB: here, including correction if any unlike Griffies (to avoid 2 variables)
   if ( ("hfds" in ncS.data_vars) ):
     HFS = ncS.hfds
   else:
     print('    WARNING :   No data found for HFS  -->  filled with NaNs')
     HFS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )

   # Water flux entering the ocean due to surface correction (SSS restoring)
   # (= wfocorr in Griffies 2016 section K2) [kg m-2 s-1]
   # - if no vprec, use salt restoring, which is converted into equiv. FW
   if ("vprec" in ncS.data_vars):
     WFOCORR = ncS.vprec
   elif ("salt_flux_added" in ncS.data_vars) and ("salt" in ncT.data_vars):
       # convert salt flux due to salinity restoring into equivalent FW flux (don't forget sign)
     WFOCORR = -1*ncS.salt_flux_added * 1000 / ncT.salt.isel(z = 0)
   else:
     print('    WARNING :   No data found for WFOCORR  -->  filled with NaNs')
     WFOCORR = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )
    
   # Water flux entering the ocean due to rainfall, snowfall, condensation - evap, 
   # river runoff [kg m-2 s-1]  
   # (= pr+prs+evs+friver in Griffies 2016, section K2)
   if ( "wfo" in ncS.data_vars ) & ("fsitherm" in ncS.data_vars) &  ("mass_flux" in ncS.data_vars) & ("areacello" in ncM.data_vars) & ("ficeberg" in ncS.data_vars):
       #("vprec" in ncS.data_vars) &
     WFOAT = ncS.wfo - ncS.fsitherm*maskT - (ncS.mass_flux/rename_dimensions(ncM.areacello)*maskT) - ncS.ficeberg*maskT
   else:
     print('    WARNING :   No data found for WFOAT  -->  filled with NaNs')
     WFOAT = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'y', 'x'] )
  
   #----------
   # Reduce the size of ocean dataset
 
   [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region=region)
   wdeg = 1 #1.5 * np.max([dxT.where(latT <= latmax).max(),dyT.where(latT <= latmax).max()]) / 6.37e6 * 180. / np.pi / np.cos(0.5*(latmin+latmax)*np.pi/180.)
   lonmin=lonmin - wdeg # take a bit more for interpolation
   lonmax=lonmax + wdeg
   latmin=latmin - wdeg
   latmax=latmax + wdeg

   condT2d = ( (latT >= latmin) & (latT <= latmax) & (lonT >= lonmin) & (lonT <= lonmax) )

   # Parallelised, but assumes lat/lon grid and no crossing over edge of domain
   
   imax = (condT2d.where(condT2d>0)*0+condT2d.x).max('y').argmax('x').values
   imin = (condT2d.where(condT2d>0)*0+condT2d.x).max('y').argmin('x').values
   jmax = (condT2d.where(condT2d>0)*0+condT2d.y).max('x').argmax('y').values
   jmin = (condT2d.where(condT2d>0)*0+condT2d.y).max('x').argmin('y').values


   #otherwise, use commented code

   # for ii in np.arange(latT.shape[1]):
   #    if ( np.sum(condT2d.isel(x=ii).values) == 0 ):
   #      imin=ii
   #    else:
   #      imin=ii
   #      break
   # for ii in np.arange(latT.shape[1]-1,0,-1):
   #    if ( np.sum(condT2d.isel(x=ii).values) == 0 ):
   #      imax=ii
   #    else:
   #      imax=ii
   #      break
   # for jj in np.arange(latT.shape[0]):
   #    if ( np.sum(condT2d.isel(y=jj).values) == 0 ):
   #      jmin=jj
   #    else:
   #      jmin=jj
   #      break
   # for jj in np.arange(latT.shape[0]-1,0,-1):
   #    if ( np.sum(condT2d.isel(y=jj).values) == 0 ):
   #      jmax=jj
   #    else:
   #      jmax=jj
   #      break

   print('Reducing domain to useful part assuming regular lat-lon grid, i.e. : ',[imin,imax,jmin,jmax])

   #----------
   # Interpolate 3d fileds from sigma to Z
   # modified from ROMS version
   #----------

   print('    Interpolating to Z levels...')

   #newdepth = np.append( np.arange(0.,70.,10.), np.arange(80.,1520.,20.) )
   newdepth = np.append(np.array([0.]), np.append( np.arange(10.,70.,10.), np.append(np.arange(80.,1520.,20.),np.arange(1600.,5000.,100.) )))
   # small array for testing
   #newdepth = np.append(np.array([0,50]),np.arange(100,1510,100))

   mz = np.size(newdepth)

   # compute middle-of-layer-depth from layer interface heights e (assumes e not varying with time!)
   ztmp = 0.5*(ncT.e.isel(time = 0).roll(z_i = -1) + ncT.e.isel(time = 0))
   ztmp = ztmp[:-1,...]
   # rename vertical coordinate to match temperature/other variables
   ztmp['z_i'] = THETAO.z.values; ztmp = ztmp.rename({'z_i':'z'})
   # if in open ocean, replace surface value with e = 0
   ztmp[0,...] = ztmp.isel(z = 0).where(SFTFLF > 0).fillna(0) 
   # fix sign and take smaller region
   ztmp=ztmp*(-1); ztmp = ztmp.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   
   # sea fraction in Z space:
   DepTUV = xr.DataArray( newdepth, dims=['z'], coords=dict( z=(['z'], newdepth) ) )
   DEPTHOT = DEPTHO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   DEPFLFT = DEPFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
    #arithmetic mean
   DEPTHOU = xr.DataArray( 0.5 * (DEPTHOT.values+DEPTHOT.shift(x=1).values), dims=['y', 'x'] )
   DEPFLFU = xr.DataArray( 0.5 * (DEPFLFT.values+DEPFLFT.shift(x=1).values), dims=['y', 'x'] )
   DEPTHOV = xr.DataArray( 0.5 * (DEPTHOT.values+DEPTHOT.shift(y=1).values), dims=['y', 'x'] )
   DEPFLFV = xr.DataArray( 0.5 * (DEPFLFT.values+DEPFLFT.shift(y=1).values), dims=['y', 'x'] )
    
   LEVOFT = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'y', 'x'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFT = LEVOFT.where( ( (DepTUV==0.) & (DEPFLFT==0.) ) | ( (DEPTHOT>=DepTUV) & (DEPFLFT<DepTUV) ), 0.e0 ) \
            * rename_dimensions(ncM.wet).isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   LEVOFT = LEVOFT.where(LEVOFT<1.e18).fillna(0)
   LEVOFU = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'y', 'x'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFU = LEVOFU.where( ( (DepTUV==0.) & (DEPFLFU==0.) ) | ( (DEPTHOU>=DepTUV) & (DEPFLFU<DepTUV) ), 0.e0 ) \
            * rename_dimensions(ncM.wet_u).isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   LEVOFU = LEVOFU.where(LEVOFU<1.e18).fillna(0)
   LEVOFV = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'y', 'x'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFV = LEVOFV.where( ( (DepTUV==0.) & (DEPFLFV==0.) ) | ( (DEPTHOV>=DepTUV) & (DEPFLFV<DepTUV) ), 0.e0 ) \
            * rename_dimensions(ncM.wet_v).isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   LEVOFV = LEVOFV.where(LEVOFV<1.e18).fillna(0)

   ###
   #  Calculate new variable
   ###
   # Ocean fraction at surface:
   SFTOF = LEVOFT.isel(z=0)*maskT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   ZMODT = xr.DataArray( ztmp, dims=['z', 'y', 'x'] )
   ZMODU = xr.DataArray( ztmp, dims=['z', 'y', 'x'] )
   ZMODV = xr.DataArray( ztmp, dims=['z', 'y', 'x'] )
   SO_red = SO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   THETAO_red = THETAO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   UX_red = UX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))
   VY_red = VY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1))

   ZMODU['x'] = UX_red.x
   ZMODV['y'] = VY_red.y
    
   SO_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   THETAO_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   UX_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   VY_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
    
   nz = np.size(newdepth)
   print(newdepth)
   newdepth = np.append(np.array([-1.0e-10]), np.append( np.arange(10.,70.,10.), np.append(np.arange(80.,1520.,20.),np.arange(1600.,5000.,100.) ))) #otherwise surface layer doesn't work, for some reason
   #newdepth = np.append(np.array([-1.0e-10,50]),np.arange(100,1510,100)) #otherwise surface layer doesn't work, for some reason
   print(newdepth)
   for kk in np.arange(np.size(newdepth)):
     print(f"\nDepth {kk+1}/{nz}")
     rss = psutil.Process(os.getpid()).memory_info().rss / 1024**3
     print(f"RSS = {rss:.1f} GB", flush=True)
     print('...')

     tmpaT = ZMODT - newdepth[kk]
     tmpaT = tmpaT.where( tmpaT >= 0., 1.e36 )
     msksupT = tmpaT.where( tmpaT==tmpaT.min('z'), 0.e0 )
     msksupT = msksupT.where( msksupT==0.e0, 1.e0 )
     dzsupT = tmpaT.where( tmpaT==tmpaT.min('z'), 0.e0 ).sum('z')
     tmpbT = ZMODT - newdepth[kk]
     tmpbT = tmpbT.where( tmpbT < 0., -1.e36 )
     mskinfT = tmpbT.where( tmpbT==tmpbT.max('z'), 0.e0 )
     mskinfT = mskinfT.where( mskinfT==0.e0, 1.e0 )
     dzinfT = -tmpbT.where( tmpbT==tmpbT.max('z'), 0.e0 ).sum('z')

     tmpaU = ZMODU - newdepth[kk]
     tmpaU = tmpaU.where( tmpaU >= 0., 1.e36 )
     msksupU = tmpaU.where( tmpaU==tmpaU.min('z'), 0.e0 )
     msksupU = msksupU.where( msksupU==0.e0, 1.e0 )
     dzsupU = tmpaU.where( tmpaU==tmpaU.min('z'), 0.e0 ).sum('z')
     tmpbU = ZMODU - newdepth[kk]
     tmpbU = tmpbU.where( tmpbU < 0., -1.e36 )
     mskinfU = tmpbU.where( tmpbU==tmpbU.max('z'), 0.e0 )
     mskinfU = mskinfU.where( mskinfU==0.e0, 1.e0 )
     dzinfU = -tmpbU.where( tmpbU==tmpbU.max('z'), 0.e0 ).sum('z')
     
     tmpaV = ZMODV - newdepth[kk]
     tmpaV = tmpaV.where( tmpaV >= 0., 1.e36 )
     msksupV = tmpaV.where( tmpaV==tmpaV.min('z'), 0.e0 )
     msksupV = msksupV.where( msksupV==0.e0, 1.e0 )
     dzsupV = tmpaV.where( tmpaV==tmpaV.min('z'), 0.e0 ).sum('z')
     tmpbV = ZMODV - newdepth[kk]
     tmpbV = tmpbV.where( tmpbV < 0., -1.e36 )
     mskinfV = tmpbV.where( tmpbV==tmpbV.max('z'), 0.e0 )
     mskinfV = mskinfV.where( mskinfV==0.e0, 1.e0 )
     dzinfV = -tmpbV.where( tmpbV==tmpbV.max('z'), 0.e0 ).sum('z')

     #getting weird values in land if I don't do this:
     dzsupT = dzsupT.where(dzinfT<1e18).fillna(0)
     dzsupU = dzsupU.where(dzinfU<1e18).fillna(0)
     dzsupV = dzsupV.where(dzinfV<1e18).fillna(0)

     dzinfT = dzinfT.where(dzsupT<1e18).fillna(0)
     dzinfU = dzinfU.where(dzsupU<1e18).fillna(0)
     dzinfV = dzinfV.where(dzsupV<1e18).fillna(0)

     # time varying 3d fields:
     for ll in np.arange(mtime):
        print(f"\rTime step {ll+1}/{mtime} ({100*(ll+1)/mtime:5.1f}%)", end="", flush=True)

        Saaa = msksupT*SO_red.isel(time=ll)
        Sbbb = mskinfT*SO_red.isel(time=ll)
        Sccc = ( Saaa.sum('z')*dzinfT + Sbbb.sum('z')*dzsupT ) / (dzinfT+dzsupT)
        SO_z[ll,kk,:,:] = Sccc.values

        Taaa = msksupT*THETAO_red.isel(time=ll)
        Tbbb = mskinfT*THETAO_red.isel(time=ll)
        Tccc = ( Taaa.sum('z')*dzinfT + Tbbb.sum('z')*dzsupT ) / (dzinfT+dzsupT)
        THETAO_z[ll,kk,:,:] = Tccc.values

        Uaaa = msksupU*UX_red.isel(time=ll)
        Ubbb = mskinfU*UX_red.isel(time=ll)
        Uccc = ( Uaaa.sum('z')*dzinfU + Ubbb.sum('z')*dzsupU ) / (dzinfU+dzsupU)
        UX_z[ll,kk,:,:] = Uccc.values

        Vaaa = msksupV*VY_red.isel(time=ll)
        Vbbb = mskinfV*VY_red.isel(time=ll)
        Vccc = ( Vaaa.sum('z')*dzinfV + Vbbb.sum('z')*dzsupV ) / (dzinfV+dzsupV)
        VY_z[ll,kk,:,:] = Vccc.values
     print()


   #----------
   # Create new xarray dataset including all useful variables:
   # reshaping (x,y) as 1-dimensional (sxy)

   nxy=(jmax-jmin+1)*(imax-imin+1)

   newdepth[0]=0.0 # so that 1st level is taken at the surface without extrapolation

   time_conv=ncT.time.dtype
   if ( time_conv == 'datetime64[ns]' ):
     time_val = ncT.time.values # standard calendar 
   else:
     time_val = ncT.indexes['time'].to_datetimeindex().values # to enable dealing with non-standard calendar (e.g. noleap)

   ds = xr.Dataset(
      {
       "SO":        (["time", "z", "sxy"], np.reshape( SO_z, (mtime,mz,nxy)) ),
       "THETAO":    (["time", "z", "sxy"], np.reshape( THETAO_z, (mtime,mz,nxy)) ),
       "UX":        (["time", "z", "sxy"], np.reshape( UX_z, (mtime,mz,nxy)) ),
       "VY":        (["time", "z", "sxy"], np.reshape( VY_z, (mtime,mz,nxy)) ),
       "TAUX":      (["time", "sxy"], np.reshape( TAUX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TAUY":      (["time", "sxy"], np.reshape( TAUY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "ZOS":       (["time", "sxy"], np.reshape( ZOS.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TOB":       (["time", "sxy"], np.reshape( TOB.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SOB":       (["time", "sxy"], np.reshape( SOB.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "FLANDICE":  (["time", "sxy"], np.reshape( FLANDICE.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "FSITHERM":  (["time", "sxy"], np.reshape( FSITHERM.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "LIBMASSBFFL": (["time", "sxy"], np.reshape( LIBMASSBFFL.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "DYDRFLF":   (["time", "sxy"], np.reshape( DYDRFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "THDRFLF":   (["time", "sxy"], np.reshape( THDRFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HADRFLF":   (["time", "sxy"], np.reshape( HADRFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "MSFTBAROT": (["time", "sxy"], np.reshape( MSFTBAROT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HFS":       (["time", "sxy"], np.reshape( HFS.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOAT":     (["time", "sxy"], np.reshape( WFOAT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOCORR":   (["time", "sxy"], np.reshape( WFOCORR.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SICONC":    (["time", "sxy"], np.reshape( SICONC.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVOL":     (["time", "sxy"], np.reshape( SIVOL.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIUX":      (["time", "sxy"], np.reshape( SIUX.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVY":      (["time", "sxy"], np.reshape( SIVY.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "LEVOFT":    (["z", "sxy"], np.reshape( LEVOFT.values, (mz,nxy)) ),
       "LEVOFU":    (["z", "sxy"], np.reshape( LEVOFU.values, (mz,nxy)) ),
       "LEVOFV":    (["z", "sxy"], np.reshape( LEVOFV.values, (mz,nxy)) ),
       "SFTOF":     (["sxy"], np.reshape( SFTOF.values, nxy) ),
       "SFTFLF":    (["sxy"], np.reshape( SFTFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DEPFLF":    (["sxy"], np.reshape( DEPFLF.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DEPTHO":    (["sxy"], np.reshape( DEPTHO.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKT":   (["sxy"], np.reshape( DOMMSKT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKU":   (["sxy"], np.reshape( DOMMSKU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKV":   (["sxy"], np.reshape( DOMMSKV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonT":      (["sxy"], np.reshape( lonT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonU":      (["sxy"], np.reshape( lonU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "lonV":      (["sxy"], np.reshape( lonV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latT":      (["sxy"], np.reshape( latT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latU":      (["sxy"], np.reshape( latU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "latV":      (["sxy"], np.reshape( latV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "dxT":       (["sxy"], np.reshape( dxT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "dyT":       (["sxy"], np.reshape( dyT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaT":    (["sxy"], np.reshape( thetaT.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaU":    (["sxy"], np.reshape( thetaU.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "thetaV":    (["sxy"], np.reshape( thetaV.isel(x=slice(imin,imax+1),y=slice(jmin,jmax+1)).values, nxy) ),
       "depTUV":    (['z'], newdepth)
      },
      coords={
      "time": time_val,
      "z": newdepth
      },
      attrs={
      "original_minlat": domain_minlat,
      "original_maxlat": domain_maxlat,
      "original_minlon": domain_minlon,
      "original_maxlon": domain_maxlon
      },
   )

   print('    Load duration: ',datetime.now() - startTime)

   return ds

 
