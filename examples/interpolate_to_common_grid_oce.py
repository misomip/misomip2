###################################################################################
# Script to interpolate native ocean model outputs to the standard MISOMIP2 grids.
#
# Notes:
#     * interpolations are done on 1-dimensional vectors, so the script should be
#       easily adaptable to non-structured grids.
#     * needs some modificationds to work with sigma or HAL coordinates.
#     
# History: 
#     * 2021-03: initial code + tests for NEMO & MITgcm outputs [N. Jourdain, IGE-CNRS]
#
###################################################################################
import numpy as np
import xarray as xr
import sys
sys.path.append("/Users/jourdain/MY_SCRIPTS")
#sys.path.append("/Users/nakayama/Documents/GitHub/")
import misomip2.preproc as mp
import os
from datetime import datetime

startTime = datetime.now()

np.seterr(divide='ignore', invalid='ignore') # to avoid warning due to divide by zero

#--------------------------------------------------------------------------
# 0- General information:

# Choose a test case:
#test_case='NEMO_test'
#test_case='MITGCM_test'
test_case='ROMS_test'
#test_case='eORCA025_test'

reg='Amundsen' # 'Amundsen' or 'Weddell'

data_dir='models/oce/'+test_case

missval=9.969209968386869e36 # missing value in created netcdf files

epsfr=1. # cut sea, sea-ice or ice-shelf fractions below epsfr (in %)
         # (useful to avoid extreme values that can occur with very small fractions)

filmis=False  # keep False (may become a deprecitated option)
              # if True, will fill nans (missing values in input file or values that can't be triangulated) 
              # through nearest-neighbor interpolation. This can be useful for the edges of non-structured grids.

#--------------------------------------------------------------------------
# 1- Files and variables

# loading an xarray dataset containing all useful variables with (x,y) reshaped 
# as 1 dimension in case of original structured grid :

if ( (test_case[0:4] == 'NEMO') | (test_case[0:5] == 'eORCA') ):

   if ( test_case[0:4] == 'NEMO' ):
     model='NEMO3.6'        # model name, possibly including a version number
     original_sim_name='random Amundsen test case'
   else:
     model='NEMO4.0'        # model name, possibly including a version number
     original_sim_name='random Global test case'

   institute='IGE-CNRS-UGA' # name of the institute(s) that produced the simulation
                            #    (use "-" rather than "_" for multiple entities)
   abc='a'                  # single letter used to distinguish multiple set-up produced by a given institute
                            #    (e.g., difference or model parameters, resolution, initial states or boundary conditions)
   exp='Ocean-A1'           # MISOMIP2 experiment name (e.g., Ocean-A1, Ocean-W1, IceOcean-A1...)

   print('LOADING NEMO...')

   fs_gridT = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_T.nc' for month in np.arange(1,3)]
   fs_gridU = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_U.nc' for month in np.arange(1,3)]
   fs_gridV = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_V.nc' for month in np.arange(1,3)]
   fs_ice   = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_icemod.nc' for month in np.arange(1,3)]
   if ( test_case[0:4] == 'NEMO' ):
     f_mesh   = data_dir+'/mesh_mask_test.nc'
     f_bathy  = data_dir+'/bathy_meter_test.nc'
     fs_SBC   = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_SBC.nc' for month in np.arange(1,3)]
     # Barotropic Streamfunction calculated at U points using the cdfpsi function which is part of the cdftools
     #    (see https://github.com/meom-group/CDFTOOLS):
     fs_BSF   = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_psi.nc' for month in np.arange(1,3)]
   else:
     f_mesh   = data_dir+'/eORCA025_test_mesh_mask.nc'
     f_bathy  = 'None'
     fs_SBC   = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_flxT.nc' for month in np.arange(1,3)]
     fs_BSF   = 'None'

   oce = mp.load_oce_mod_nemo( file_mesh_mask=f_mesh, file_bathy=f_bathy, files_gridT=fs_gridT,\
                  files_gridU=fs_gridU, files_gridV=fs_gridV, files_SBC=fs_SBC, files_ice=fs_ice,\
                  files_BSF=fs_BSF, rho0=1026.0, teos10=True )

elif ( test_case[0:6] == 'MITGCM' ):

   model='MITGCM'           # model name, possibly including a version number
   institute='UNN'          # name of the institute(s) that produced the simulation
                            #    (use "-" rather than "_" for multiple entities)
   abc='a'                  # single letter used to distinguish multiple set-up produced by a given institute
                            #    (e.g., difference or model parameters, resolution, initial states or boundary conditions)
   exp='Ocean-A1'           # MISOMIP2 experiment name (e.g., Ocean-A1, Ocean-W1, IceOcean-A1...)
   original_sim_name='random high-res Amundsen test case'

   print('LOADING MITGCM...')
   fT = data_dir+'/'+test_case+'_THETA.nc'
   fS = data_dir+'/'+test_case+'_SALT.nc'
   fU = data_dir+'/'+test_case+'_UVEL.nc'
   fV = data_dir+'/'+test_case+'_VVEL.nc'

   oce = mp.load_oce_mod_mitgcm( files_T=fT, files_S=fS, files_U=fU, files_V=fV, rho0=1026.0, teos10=False )

elif ( test_case[0:4] == 'ROMS' ):

   model='ROMS'             # model name, possibly including a version number
   institute='UTAS'         # name of the institute(s) that produced the simulation
                            #    (use "-" rather than "_" for multiple entities)
   abc='a'                  # single letter used to distinguish multiple set-up produced by a given institute
                            #    (e.g., difference or model parameters, resolution, initial states or boundary conditions)
   exp='Ocean-A1'           # MISOMIP2 experiment name (e.g., Ocean-A1, Ocean-W1, IceOcean-A1...)
   original_sim_name='random cirum-Antarctic test case'

   print('LOADING ROMS...')
   f_grid = data_dir+'/ROMS_test_grid.nc'
   f_ALL  = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'.nc' for month in np.arange(1,4)]

   oce = mp.load_oce_mod_roms( files_M=f_grid, files_T=f_ALL, rho0=1026.0, teos10=False )

else:

   sys.exit('Unknown test_case ==> Write a function to load this new test case')

#--------------------------------------------------------------------------
# 2- Global attributes of output netcdf :

def put_global_attrs(ds,experiment='TBD', avg_hor_res_73S=0.0, original_sim_name='TBD', ocean_model='TBD', institute='TBD', \
                     original_min_lat=-90.0, original_max_lat=90.0, original_min_lon=-180.0, original_max_lon=180.0):
   """ Put global attributes to the ds xarray dataset

   """
   ds.attrs['project'] = 'MISOMIP2'
   ds.attrs['contact'] = 'TBD <tbd@misomip2.com>'     # Name(s) of the person(s) who produced the simulation <email>
   ds.attrs['institute'] = institute                  # Name of the institute(s) that produced the simulation 
                                                      #     (use ``-'' to separate multiple entities)
   ds.attrs['computing_facility'] = 'TBD'             # Computing center where the simulation was run
   ds.attrs['interpolation_method'] = 'linear triangular barycentric'  # e.g. 'linear triangular barycentric', 'bi-linear',
                                                                       #      'nearest-neighbor', 'conservative', ...       
   ds.attrs['ocean_model'] = ocean_model              # Model name and version
   ds.attrs['reference'] = 'TBD'                      # Main publication or website describing the simulation or a similar one
   ds.attrs['original_sim_name'] = original_sim_name  # Original simulation name (to keep track of what was used in MISOMIP2)
   ds.attrs['experiment'] = experiment                # MISOMIP2 experiment; e.g. 'Ocean-A1', 'Ocean-W1', 'Ocean-A2', ...
   ds.attrs['bathymetry'] = 'TBD'                     # Bathymetry dataset (specify exact version), e.g. `BedMachine-1.33', `Bedmap2', 
                                                      #    `RTopo-2.0.4', `Merge of Millan et al. (2017) and Bedmap2'
   ds.attrs['ice_shelf_draft'] = 'TBD'                # Dataset for the depth of ice-shelf base (similar to previous)
   ds.attrs['atmosphere'] = 'TBD'                     # Atmospheric forcing, with a reference, e.g., `ERA5 (Hersbach et al. 2020)',
                                                      #    `ERAint (Dee et al. 2011)', `JRA55do (Tsujino et al. 2018)',
                                                      #    `MARv3.9.3 (Donat-Magnin et al. 2020)'
   ds.attrs['iceberg'] = 'TBD'                        # Method used to account for melting icebergs, with a reference,
                                                      #     e.g. `None', `Lagrangian model (Martin and Adcroft 2010)',
                                                      #          `Prescribed freshwater (Merino et al. 2016)',
                                                      #          `Prescribed Freshwater and Heat (Merino et al. 2016)' 
   ds.attrs['sea_ice'] = 'TBD'                        # Method used to simulate or prescribe the ocean--sea-ice interaction,
                                                      #     with a reference, e.g. 'Dynamics-Thermodynamics Model (Hibler 1979)',
                                                      #                            'Thermodynamics Model (Bitz and Lipscomb 1999)',
                                                      #                            'Prescribed Freshwater and Heat'
   ds.attrs['ocean_lateral_bdy'] = 'TBD'              # Type of lateral boundary conditions, e.g. 'Simulation (Merino et al. 2018)',
                                                      #     'Reanalysis (Mazloff et al. 2016)', 'Observations (Locarnini et al. 2018)',
                                                      #     'Corrected simulation (explain method)', 'None' 
   ds.attrs['tides'] = 'TBD'                          # Method used to account for the effect of tides on ice-shelf melt, and dataset
                                                      #     e.g. `Barotropic tidal harmonics prescribed at lateral boundaries (CATS)',
                                                      #          `Forced by a tidal potential', `None',
                                                      #          `Parameterized through uniform tidal velocity (utide=0.1 m s-1)',
                                                      #          'Parameterized through non-uniform tidal velocity (FES2012)'
   ds.attrs['vertical_coordinate'] = 'TBD'            # e.g. `Geopotential (Z)', `Stretched Geopotential (Zstar)', `Pressure (P)',
                                                      #      `Stretched Pressure (P*)', `Isopycnal', `Terrain-Following (Sigma)',
                                                      #      `Arbitrary Lagrangian-Eulerian (ALE)'
   ds.attrs['is_melt_param'] = 'TBD'                  # Parameterization used to calculate ice-shelf basal melt rates,
                                                      #      e.g. `3-equation (constant gamma)', `3-equation (velocity-dependent gamma)',
                                                      #           `3-equation (stability and velocity-dependent gamma)'
   ds.attrs['eos'] = 'TBD'                            # Equation of state in original simulation, e.g. `TEOS10', `EOS80', `linear'
   ds.attrs['advection'] = 'TBD'                      # Brief description of the momentum and tracer advection schemes 
                                                      #      (centered, third-order with limiter, etc)
   ds.attrs['horizontal_mixing'] = 'TBD'              # Brief description of how ``horizontal'' mixing was performed 
                                                      #      (harmonic, biharmonic, etc.; within model levels, along geopotentials, 
                                                      #       along isopycnals, etc.; using the Gent-McWilliams parameterization; etc)
   ds.attrs['vertical_mixing'] = 'TBD'                # Brief description of how ``vertical'' mixing was performed (constant diffusivity,
                                                      #       k-profile parameterization, etc.; harmonic, biharmonic, etc)
   ds.attrs['convection'] = 'TBD'                     # Brief description of the procedure for handling convection, 
                                                      #       e.g. `Explicitly modeled', `Parameterized using enhanced vertical mixing'
   ds.attrs['avg_hor_res_73S'] = avg_hor_res_73S      # Average horizontal resolution (m) at 73degS in the MISOMIP2 domain (average of x and y resolution)
   ds.attrs['original_min_lat'] = original_min_lat    # Minimum latitude of the original domain in [-90:90]
   ds.attrs['original_max_lat'] = original_max_lat    # Minimum latitude of the original domain in [-90:90]
   ds.attrs['original_min_lon'] = original_min_lon    # Minimum longitude of the original domain in [-180:180]
   ds.attrs['original_max_lon'] = original_max_lon    # Minimum longitude of the original domain in [-180:180]

#--------------------------------------------------------------------------
# 3a- Interpolate to common 3d grid :


# Characteristics of MISOMIP 3d grid:
[lon_miso,lat_miso,dep_miso] = mp.generate_3d_grid_oce(region=reg)
mlon = np.size(lon_miso)
mlat = np.size(lat_miso)
mdep = np.size(dep_miso)
lon2d_miso, lat2d_miso = np.meshgrid( lon_miso, lat_miso )
lon_miso1d = np.reshape( lon2d_miso, mlon*mlat )
lat_miso1d = np.reshape( lat2d_miso, mlon*mlat )

# model coordinates:
lonT=oce.lonT ; lonU=oce.lonU ; lonV=oce.lonV
latT=oce.latT ; latU=oce.latU ; latV=oce.latV

# scale factor to interpolate in the (lon*cos(<lat>),lat) space:
aa = np.cos( np.mean(lat_miso)*np.pi/180. )

# Some quantities needed to define the global attributes:
res_73S = 0.5 * (  oce.dxT.where(  (latT < -72.5) & (latT >= -73.5) & (lonT >= lon_miso.min()) & (lonT <= lon_miso.max()) \
                                 & (oce.LEVOFT.isel(z=0)>99.)   ).mean().values \
                 + oce.dyT.where(  (latT < -72.5) & (latT >= -73.5) & (lonT >= lon_miso.min()) & (lonT <= lon_miso.max()) \
                                 & (oce.LEVOFT.isel(z=0)>99.) ).mean().values )
print('Average horizontal resolution at 73S :',res_73S)

# 3d sea fraction
LEVOFT=oce.LEVOFT
LEVOFU=oce.LEVOFU
LEVOFV=oce.LEVOFV
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
seafracU=LEVOFU.max('z',skipna=True)
seafracV=LEVOFV.max('z',skipna=True)
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point
openseafracU=LEVOFU.isel(z=0)
openseafracV=LEVOFV.isel(z=0)

mtime = oce.time.size
yeari = str(oce.time.isel(time=0).dt.year.values)
yearf = str(oce.time.isel(time=mtime-1).dt.year.values)
monthi = str(oce.time.isel(time=0).dt.month.values).zfill(2)
monthf = str(oce.time.isel(time=mtime-1).dt.month.values).zfill(2)
period = yeari+monthi+'-'+yearf+monthf
print('Data from '+monthi+'-'+yeari+' to '+monthf+'-'+yearf)

# mask showing the original domain on the MISOMIP grid (=nan where interpolation of any of T, U, V grid is nan, =1 elsewhere):
DOMMSK_miso =               mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DOMMSKT )
DOMMSK_miso = DOMMSK_miso + mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DOMMSKU )
DOMMSK_miso = DOMMSK_miso + mp.horizontal_interp( latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DOMMSKV )

# horizontal interpolation grid roation angles :
theT = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.thetaT )
theU = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.thetaU )
theV = mp.horizontal_interp( latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.thetaV )

# Ice shelf fraction on MISOMIP grid (in [0-100], =nan beyond model domain)
SFTFLI_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SFTFLI )
SFTFLI_miso[ np.isnan(DOMMSK_miso) ] = 0.e0 # will update to missval at the end of the calculation

# Depth of floating ice on MISOMIP grid (ice-shelf draft)
DEPFLI_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DEPFLI, weight=oce.SFTFLI, skipna=True, filnocvx=filmis, threshold=epsfr )
DEPFLI_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPFLI_miso) ] = 0.e0  # will be replaced with missval later on

# Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
DEPTHO_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DEPTHO, weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
DEPTHO_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPTHO_miso) ] = 0.e0 # will be replaced with missval later on

# Sea area fraction at each vertical level on the MISOMIP grid: 
tmp_LEVT = LEVOFT.interp(z=dep_miso) # extrapolate to avoid nan at surface but better to
tmp_LEVU = LEVOFU.interp(z=dep_miso) # bring first level to zero when loading model data.
tmp_LEVV = LEVOFV.interp(z=dep_miso)

LEVOF_miso = np.zeros((mdep,mlat,mlon))
for kk in np.arange(mdep):
  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_LEVT[kk,:] )
  tzz[ np.isnan(DOMMSK_miso) ] = np.nan # will update to missval at the end of the calculation
  LEVOF_miso[kk,:,:] = tzz

# vertical interpolation of T, S, U, V to common vertical grid :

tmpxS = LEVOFT*oce.SO
tmp_SS = tmpxS.interp(z=dep_miso) / tmp_LEVT

tmpxT = LEVOFT*oce.THETAO
tmp_TT = tmpxT.interp(z=dep_miso) / tmp_LEVT

tmpxU = LEVOFU*oce.UX
tmp_UX = tmpxU.interp(z=dep_miso) / tmp_LEVU

tmpxV = LEVOFV*oce.VY
tmp_VY = tmpxV.interp(z=dep_miso) / tmp_LEVV

#----- interpolation of time-varying fields:

SO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
THETAO_miso = np.zeros((mtime,mdep,mlat,mlon)) + missval
UO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
VO_miso     = np.zeros((mtime,mdep,mlat,mlon)) + missval
ZOS_miso    = np.zeros((mtime,mlat,mlon)) + missval
TOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
SOB_miso    = np.zeros((mtime,mlat,mlon)) + missval
FICESHELF_miso = np.zeros((mtime,mlat,mlon)) + missval
DYDRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
THDRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
HADRFLI_miso = np.zeros((mtime,mlat,mlon)) + missval
MSFTBAROT_miso = np.zeros((mtime,mlat,mlon)) + missval
HFDS_miso = np.zeros((mtime,mlat,mlon)) + missval
WFOATRLI_miso = np.zeros((mtime,mlat,mlon)) + missval
WFOSICOR_miso = np.zeros((mtime,mlat,mlon)) + missval
SICONC_miso = np.zeros((mtime,mlat,mlon)) + missval
SIVOL_miso = np.zeros((mtime,mlat,mlon)) + missval
SIU_miso = np.zeros((mtime,mlat,mlon)) + missval
SIV_miso = np.zeros((mtime,mlat,mlon)) + missval
TAUUO_miso = np.zeros((mtime,mlat,mlon)) + missval
TAUVO_miso = np.zeros((mtime,mlat,mlon)) + missval

domcond = (np.isnan(DOMMSK_miso))

for ll in np.arange(mtime):

  ## fields interpolated from sea cells ( C/T grid ):

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.ZOS.isel(time=ll), \
                              weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  ZOS_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TOB.isel(time=ll), \
                              weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  TOB_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SOB.isel(time=ll), \
                              weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  SOB_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.HFDS.isel(time=ll), \
                              weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  HFDS_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.WFOATRLI.isel(time=ll), \
                              weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  WFOATRLI_miso[ll,:,:] = tzz

  ## fields interpolated from sea cells ( U & V grids ):

  tzz = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.MSFTBAROT.isel(time=ll), \
                              weight=seafracU, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  MSFTBAROT_miso[ll,:,:] = tzz

  TAUX_notrot = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TAUX.isel(time=ll), \
                                      weight=seafracU, skipna=True, filnocvx=filmis, threshold=epsfr )
  TAUY_notrot = mp.horizontal_interp( latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TAUY.isel(time=ll), \
                                      weight=seafracV, skipna=True, filnocvx=filmis, threshold=epsfr )
  TAUUO_miso[ll,:,:] = TAUX_notrot * np.cos(theU) - TAUY_notrot * np.sin(theV) # rotated to zonal
  TAUVO_miso[ll,:,:] = TAUY_notrot * np.cos(theV) + TAUX_notrot * np.sin(theU) # rotated to meridional
  TAUUO_miso[ ( np.isnan(TAUUO_miso) ) | domcond ] = missval
  TAUVO_miso[ ( np.isnan(TAUVO_miso) ) | domcond ] = missval

  ## fileds interpolated from open-sea cells :

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.WFOSICOR.isel(time=ll), \
                              weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  WFOSICOR_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SICONC.isel(time=ll), \
                              weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = np.nan # will be updated to missval later on
  SICONC_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIVOL.isel(time=ll), \
                              weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  SIVOL_miso[ll,:,:] = tzz

  ## fileds interpolated from ice-shelf cells :

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.FICESHELF.isel(time=ll), \
                              weight=oce.SFTFLI, skipna=True, filnocvx=filmis, threshold=epsfr )
  tzz[ np.isnan(tzz) | domcond ] = missval
  FICESHELF_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DYDRFLI.isel(time=ll), \
                              weight=oce.SFTFLI, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  DYDRFLI_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.THDRFLI.isel(time=ll), \
                              weight=oce.SFTFLI, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  THDRFLI_miso[ll,:,:] = tzz

  tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.HADRFLI.isel(time=ll), \
                              weight=oce.SFTFLI, skipna=True, filnocvx=filmis, threshold=epsfr ) 
  tzz[ np.isnan(tzz) | domcond ] = missval
  HADRFLI_miso[ll,:,:] = tzz

  ## fileds interpolated from sea-ice cells :

  SIUX_notrot = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIUX.isel(time=ll), \
                                      weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=filmis, threshold=epsfr )
  SIVY_notrot = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIVY.isel(time=ll), \
                                      weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=filmis, threshold=epsfr )
  SIU_miso[ll,:,:] = SIUX_notrot * np.cos(theT) - SIVY_notrot * np.sin(theT) # rotated to zonal
  SIV_miso[ll,:,:] = SIVY_notrot * np.cos(theT) + SIUX_notrot * np.sin(theT) # rotated to meridional
  SIU_miso [ np.isnan(SIU_miso) | domcond ] = missval
  SIV_miso [ np.isnan(SIV_miso) | domcond ] = missval

  for kk in np.arange(mdep):

    # horizontal interpolations of time-varying 3d fields to common horizontal grid : 

    condkk = ( (LEVOF_miso[kk,:,:] < epsfr) | domcond )

    tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_SS.isel(time=ll,z=kk), \
                                weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    SO_miso[ll,kk,:,:] = tzz

    tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_TT.isel(time=ll,z=kk), \
                                weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    THETAO_miso[ll,kk,:,:] = tzz

    UX_notrot = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_UX.isel(time=ll,z=kk), \
                                      weight=tmp_LEVU.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
    VY_notrot = mp.horizontal_interp( latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_VY.isel(time=ll,z=kk), \
                                      weight=tmp_LEVV.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
    tzz = UX_notrot * np.cos(theU) - VY_notrot * np.sin(theV) # rotated to zonal
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    UO_miso[ll,kk,:,:] = tzz
    tzz = VY_notrot * np.cos(theV) + UX_notrot * np.sin(theU) # rotated to meridional
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    VO_miso[ll,kk,:,:] = tzz


LEVOF_miso[ np.isnan(LEVOF_miso) ] = missval
SFTFLI_miso[ np.isnan(SFTFLI_miso) | np.isnan(DOMMSK_miso) ] = missval 
SICONC_miso[ np.isnan(SICONC_miso) ] = missval
DEPFLI_miso[ (DEPFLI_miso==0.e0) | np.isnan(DOMMSK_miso) ] = missval
DEPTHO_miso[ np.isnan(DOMMSK_miso) | (DEPTHO_miso==0.) ] = missval

#--------------------------------------------------------------------------
# 3b- Create new xarray dataset and save to netcdf

dsmiso3d = xr.Dataset(
    {
    "so":        (["time", "depth", "latitude", "longitude"], np.float32(SO_miso)),
    "thetao":    (["time", "depth", "latitude", "longitude"], np.float32(THETAO_miso)),
    "uo":        (["time", "depth", "latitude", "longitude"], np.float32(UO_miso)),
    "vo":        (["time", "depth", "latitude", "longitude"], np.float32(VO_miso)),
    "tauuo":     (["time", "latitude", "longitude"], np.float32(TAUUO_miso)),
    "tauvo":     (["time", "latitude", "longitude"], np.float32(TAUVO_miso)),
    "zos":       (["time", "latitude", "longitude"], np.float32(ZOS_miso)),
    "tob":       (["time", "latitude", "longitude"], np.float32(TOB_miso)),
    "sob":       (["time", "latitude", "longitude"], np.float32(SOB_miso)),
    "ficeshelf": (["time", "latitude", "longitude"], np.float32(FICESHELF_miso)),
    "dydrfli":   (["time", "latitude", "longitude"], np.float32(DYDRFLI_miso)),
    "thdrfli":   (["time", "latitude", "longitude"], np.float32(THDRFLI_miso)),
    "hadrfli":   (["time", "latitude", "longitude"], np.float32(HADRFLI_miso)),
    "msftbarot": (["time", "latitude", "longitude"], np.float32(MSFTBAROT_miso)),
    "hfds":      (["time", "latitude", "longitude"], np.float32(HFDS_miso)),
    "wfoatrli":  (["time", "latitude", "longitude"], np.float32(WFOATRLI_miso)),
    "wfosicor":  (["time", "latitude", "longitude"], np.float32(WFOSICOR_miso)),
    "siconc":    (["time", "latitude", "longitude"], np.float32(SICONC_miso)),
    "sivol":     (["time", "latitude", "longitude"], np.float32(SIVOL_miso)),
    "siu":       (["time", "latitude", "longitude"], np.float32(SIU_miso)),
    "siv":       (["time", "latitude", "longitude"], np.float32(SIV_miso)),
    "levof":     (["depth", "latitude", "longitude"], np.float32(LEVOF_miso)),
    "sftfli":    (["latitude", "longitude"], np.float32(SFTFLI_miso)),
    "depfli":    (["latitude", "longitude"], np.float32(DEPFLI_miso)),
    "deptho":    (["latitude", "longitude"], np.float32(DEPTHO_miso)),
    },
    coords={
    "longitude":np.float32(lon_miso),
    "latitude": np.float32(lat_miso),
    "depth": np.float32(dep_miso),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dsmiso3d,miss=missval)
put_global_attrs(dsmiso3d, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name, \
                 ocean_model=model, institute=institute, \
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )

file_miso3d = 'Oce3d_'+model+'-'+institute+'_'+abc+'_'+exp+'_'+period+'.nc'

print('Creating ',file_miso3d)
dsmiso3d.to_netcdf(file_miso3d,unlimited_dims="time")

del dsmiso3d
del SO_miso, THETAO_miso, UO_miso, VO_miso, ZOS_miso, FICESHELF_miso, DYDRFLI_miso, THDRFLI_miso, HADRFLI_miso, MSFTBAROT_miso
del HFDS_miso, WFOATRLI_miso, WFOSICOR_miso, SICONC_miso, SIVOL_miso, SIU_miso, SIV_miso

print('  Execution time: ',datetime.now() - startTime)

#--------------------------------------------------------------------------
# 4a- Interpolate to common section :

# Characteristics of MISOMIP section grid:
[lon_sect1d,lat_sect1d,dep_sect] = mp.generate_section_grid_oce(region=reg)
mlonlatsec = np.size(lon_sect1d)
mdepsect = np.size(dep_sect)

# Reduce input data size
wdeg = 1.5 * res_73S / ( 6.37e6 * np.pi / 180. )
lonmin_sect = np.min( lon_sect1d ) - wdeg/aa
lonmax_sect = np.max( lon_sect1d ) + wdeg/aa
latmin_sect = np.min( lat_sect1d ) - wdeg
latmax_sect = np.max( lat_sect1d ) + wdeg

cond_secT=( (oce.latT>latmin_sect) & (oce.latT<latmax_sect) & (oce.lonT>lonmin_sect) & (oce.lonT<lonmax_sect) )

# model coordinates:
lonT=oce.lonT.where(cond_secT,drop=True)
latT=oce.latT.where(cond_secT,drop=True)

# 3d sea fraction
LEVOFT=oce.LEVOFT.where(cond_secT,drop=True)
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point

# mask showing the original domain (nan where interpolation of any of T, U, V grid is nan):
DOMMSK_sect = np.squeeze(mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, oce.DOMMSKT.where(cond_secT,drop=True) ))

domcond = ( np.isnan(DOMMSK_sect) )

# Depth of floating ice on MISOMIP grid (ice-shelf draft)
DEPFLI_sect = np.squeeze( mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, oce.DEPFLI.where(cond_secT,drop=True), \
                                                weight=oce.SFTFLI.where(cond_secT,drop=True), \
                                                skipna=True, filnocvx=filmis, threshold=epsfr ) )
DEPFLI_sect[ domcond | np.isnan(DEPFLI_sect) ] = 0.e0 # will be replaced with missval later on

# Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
DEPTHO_sect = np.squeeze( mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, oce.DEPTHO.where(cond_secT,drop=True), \
                                                weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr ))
DEPTHO_sect[ domcond | np.isnan(DEPTHO_sect) ] = 0.e0 # will be replaced with missval later on

# vertical then horizontal interpolation of constant 3d fields to common grid :
tmp_LEVT = LEVOFT.interp(z=dep_sect)

LEVOF_sect = np.zeros((mdepsect,mlonlatsec))
for kk in np.arange(mdepsect):
  tzz = np.squeeze( mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_LEVT[kk,:] ) )
  tzz[ domcond ] = missval
  LEVOF_sect[kk,:] = tzz

# vertical interpolation of T, S, U, V to common vertical grid :
tmpxS = LEVOFT*oce.SO.where(cond_secT,drop=True)
tmp_SS = tmpxS.interp(z=dep_sect) / tmp_LEVT
tmpxT = LEVOFT*oce.THETAO.where(cond_secT,drop=True)
tmp_TT = tmpxT.interp(z=dep_sect) / tmp_LEVT
SO_sect     = np.zeros((mtime,mdepsect,mlonlatsec)) + missval
THETAO_sect = np.zeros((mtime,mdepsect,mlonlatsec)) + missval

for ll in np.arange(mtime):
  for kk in np.arange(mdepsect):

    condkk = ( (LEVOF_sect[kk,:] < epsfr) | domcond )

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_SS.isel(time=ll,z=kk), \
                                           weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    SO_sect[ll,kk,:] = tzz

    tzz = np.squeeze(mp.horizontal_interp( latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_TT.isel(time=ll,z=kk), \
                                           weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr ) )
    tzz[ condkk | (np.isnan(tzz)) ] = missval
    THETAO_sect[ll,kk,:] = tzz

DEPFLI_sect[ DEPFLI_sect==0.e0 ] = missval
DEPTHO_sect[ np.isnan(DOMMSK_sect) | (DEPTHO_sect==0.) ] = missval

#--------------------------------------------------------------------------
# 4b- Create new xarray dataset and save to netcdf

dssect = xr.Dataset(
    {
    "so":        (["time", "depth", "x"], np.float32(SO_sect)),
    "thetao":    (["time", "depth", "x"], np.float32(THETAO_sect)),
    "levof":     (["depth", "x"], np.float32(LEVOF_sect)),
    "longitude": (["x"], np.float32(lon_sect1d)),
    "latitude": (["x"], np.float32(lat_sect1d)),
    },
    coords={
    "depth": np.float32(dep_sect),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dssect,miss=missval)
put_global_attrs(dssect, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name,\
                 ocean_model=model, institute=institute, \
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )

file_sect = 'OceSec_'+model+'-'+institute+'_'+abc+'_'+exp+'_'+period+'.nc'

print('Creating ',file_sect)
dssect.to_netcdf(file_sect,unlimited_dims="time")

del dssect
del SO_sect, THETAO_sect

print('  Execution time: ',datetime.now() - startTime)

#--------------------------------------------------------------------------
# 5a- Interpolate to common mooring location :
#
#     For this mooring, we take the closest grid point that is not under 
#     an ice shelf (observational mooring in front of a real ice shelf
#     can be located in the cavity of some models).

# Characteristics of MISOMIP mooring:
[lon_moor0d,lat_moor0d,dep_moor] = mp.generate_mooring_grid_oce(region=reg)
mdepmoor = np.size(dep_moor)

# Reduce input data size
#   NB: for some reason, a small number of points (low wdeg below) may lead to the following error here:
#       scipy.spatial.qhull.QhullError: QH6214 qhull input error: not enough points(1) to construct initial simplex (need 4) 
#       If it happens with your dataset, increase it further...
#wdeg = 1.5 * np.sqrt( res_73S / 6.37e6 * 180. / np.pi / aa )
wdeg = 3.0 * res_73S / ( 6.37e6 * np.pi / 180. )
lonmin_moor = np.min( lon_moor0d ) - wdeg/aa
lonmax_moor = np.max( lon_moor0d ) + wdeg/aa
latmin_moor = np.min( lat_moor0d ) - wdeg
latmax_moor = np.max( lat_moor0d ) + wdeg

cond_mooT=( (oce.latT>latmin_moor) & (oce.latT<latmax_moor) & (oce.lonT>lonmin_moor) & (oce.lonT<lonmax_moor) )

# model coordinates:
lonT=oce.lonT.where(cond_mooT,drop=True)
latT=oce.latT.where(cond_mooT,drop=True)

# 3d sea fraction
LEVOFT=oce.LEVOFT.where(cond_mooT,drop=True)
# 2d sea fraction (including under-ice-shelf seas)
seafracT=LEVOFT.max('z',skipna=True) # 100 if any ocean point, including under ice shelves, =0 or =nan if no ocean point
# 2d sea fraction (NOT including under-ice-shelf seas)
openseafracT=LEVOFT.isel(z=0) # 100 only if open ocean point

mtime = np.shape(oce.SO)[0]

# Coordinates of closest point:
indxy = ((lat_moor0d-latT.where(openseafracT >= 50.))**2+((lon_moor0d-lonT.where(openseafracT >= 50.))*aa)**2).argmin()

print('Closest point of open sea near mooring of (lon,lat) = ',lon_moor0d,lat_moor0d)
print('   is at model point of (lon,lat) = ',lonT.where(openseafracT >= 50.).isel(sxy=indxy).values, latT.where(openseafracT >= 50.).isel(sxy=indxy).values )
DOMMSK_moor = oce.DOMMSKT.where(cond_mooT,drop=True).where(openseafracT >= 50.).isel(sxy=indxy).values
DEPFLI_moor = oce.DEPFLI.where(cond_mooT,drop=True).where(openseafracT >= 50.).isel(sxy=indxy).values
DEPFLI_moor[ np.isnan(DOMMSK_moor) ] = missval
DEPTHO_moor = oce.DEPTHO.where(cond_mooT,drop=True).where(openseafracT >= 50.).isel(sxy=indxy).values
DEPTHO_moor[ np.isnan(DOMMSK_moor) ] = missval

tmp_LEVT = LEVOFT.interp(z=dep_moor)
LEVOF_moor = np.zeros((mdepmoor))
for kk in np.arange(mdepmoor):
  tzz = tmp_LEVT[kk,:].where(openseafracT >= 50.).isel(sxy=indxy).values
  tzz[ np.isnan(DOMMSK_moor) ] = missval
  LEVOF_moor[kk] = tzz

# vertical interpolation of T, S to common vertical grid :
tmpxS = LEVOFT*oce.SO.where(cond_mooT,drop=True)
tmp_SS = tmpxS.interp(z=dep_moor) / tmp_LEVT
tmpxT = LEVOFT*oce.THETAO.where(cond_mooT,drop=True)
tmp_TT = tmpxT.interp(z=dep_moor) / tmp_LEVT
SO_moor     = np.zeros((mtime,mdepmoor)) + missval
THETAO_moor = np.zeros((mtime,mdepmoor)) + missval

for ll in np.arange(mtime):
  for kk in np.arange(mdepmoor):

    tzz = tmp_SS[kk,:,ll].where(openseafracT >= 50.).isel(sxy=indxy).values
    tzz[ np.isnan(DOMMSK_moor) | (LEVOF_moor[kk] < epsfr) ] = missval
    SO_moor[ll,kk] = tzz

    tzz = tmp_TT[kk,:,ll].where(openseafracT >= 50.).isel(sxy=indxy).values
    tzz[ np.isnan(DOMMSK_moor) | (LEVOF_moor[kk] < epsfr) ] = missval
    THETAO_moor[ll,kk] = tzz

#--------------------------------------------------------------------------
# 5b- Create new xarray dataset and save to netcdf

dsmoor = xr.Dataset(
    {
    "so":        (["time", "depth"], np.float32(SO_moor)),
    "thetao":    (["time", "depth"], np.float32(THETAO_moor)),
    "levof":     (["depth"], np.float32(LEVOF_moor)),
    },
    coords={
    "depth": np.float32(dep_moor),
    "time": oce.time.values
    },
)

mp.add_standard_attributes_oce(dsmoor,miss=missval)
put_global_attrs(dsmoor, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name,\
                 ocean_model=model, institute=institute, \
                 original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                 original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon') )
dsmoor.attrs['mooring_longitude'] = np.float32(lon_moor0d)
dsmoor.attrs['mooring_latitude'] = np.float32(lat_moor0d)

file_moor = 'OceMoor_'+model+'-'+institute+'_'+abc+'_'+exp+'_'+period+'.nc'

print('Creating ',file_moor)
dsmoor.to_netcdf(file_moor,unlimited_dims="time")

print('  Execution time: ',datetime.now() - startTime)
