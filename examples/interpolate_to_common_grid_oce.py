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
sys.path.append("/ccc/work/cont003/gen6035/mathiotp/MISOMIP2/")
#sys.path.append("/Users/nakayama/Documents/GitHub/")
import misomip2.preproc as mp
import os
from datetime import datetime

# script input managment
import argparse
import yaml

from scipy.spatial import Delaunay
from scipy.interpolate import LinearNDInterpolator

startTime = datetime.now()

np.seterr(divide='ignore', invalid='ignore') # to avoid warning due to divide by zero

#--------------------------------------------------------------------------
# Function to define global attributes of output netcdf :
#--------------------------------------------------------------------------
def put_global_attrs(ds, experiment='TBD', avg_hor_res_73S=0.0, original_sim_name='TBD', ocean_model='TBD', institute='TBD', \
                     original_min_lat=-90.0, original_max_lat=90.0, original_min_lon=-180.0, original_max_lon=180.0, metadata={}):
   """ Put global attributes to the ds xarray dataset

   """
   
   ds.attrs = metadata.copy()

   ds.attrs['institute'] = institute
   ds.attrs['interpolation_method'] = 'linear triangular barycentric'  # e.g. 'linear triangular barycentric', 'bi-linear',
   ds.attrs['ocean_model']          = ocean_model
   ds.attrs['original_sim_name']    = original_sim_name
   ds.attrs['experiment']           = experiment
   ds.attrs['avg_hor_res_73S']      = avg_hor_res_73S      # Average horizontal resolution (m) at 73degS in the MISOMIP2 domain (average of x and y resolution)
   ds.attrs['original_min_lat']     = original_min_lat    # Minimum latitude of the original domain in [-90:90]
   ds.attrs['original_max_lat']     = original_max_lat    # Minimum latitude of the original domain in [-90:90]
   ds.attrs['original_min_lon']     = original_min_lon    # Minimum longitude of the original domain in [-180:180]
   ds.attrs['original_max_lon']     = original_max_lon    # Minimum longitude of the original domain in [-180:180]

#--------------------------------------------------------------------------
# Function to manage script arguments
#--------------------------------------------------------------------------
def load_arguments():
    parser = argparse.ArgumentParser(description="Process ISOMIP2 data.")
    parser.add_argument("--config"  , type=str, default='test_cases/oce/NEMO_test.yml', help="Path to YAML configuration file (optional).")

# may be useful for production of the full simulation set
    parser.add_argument("--yearb"   , type=int, default=None, help="start year")
    parser.add_argument("--yeare"   , type=int, default=None, help="end year")
    
    return parser.parse_args()

#--------------------------------------------------------------------------
# Function to load yaml file
#--------------------------------------------------------------------------
def load_config(yaml_file, section):
    """Load configuration from YAML file."""
    if yaml_file and os.path.exists(yaml_file):
        with open(yaml_file, "r") as f:
            yaml_config = yaml.safe_load(f)
            # Merge YAML config into defaults
            config = yaml_config[section]
    else:
        print(f'card {yaml_file} not found')
        sys.exit(42)

    return config

def generate_filenames(ddir: str, sim_name: str, ftype: str, yearb: int, yeare: int) -> list:
    return [f"{ddir}/{year}/{sim_name}_y{year}m{month:02d}.1m_{ftype}.nc" for year in range(yearb, yeare+1) for month in range(1, 13)]

def compute_tri_and_mask(data,x,y,nandrop=False):
    '''
    compute reusable triangulation and nanmask for interpolation
    '''
    mask = (~np.isnan(data)) & (~np.isinf(data))
    tri = Delaunay(np.column_stack((x.where( mask, drop=nandrop ), latT.where( mask, drop=nandrop ))))  # Reusable triangulation
    return tri, mask
 
#--------------------------------------------------------------------------
# 0- General information:
#--------------------------------------------------------------------------

# load arguments
args = load_arguments()

# update regions, sections and moorings values
tasks = load_config(args.config, 'processing')

# update simulation value
config = load_config(args.config, 'simulation')

# initialise variables
test_case        =config['test_case']
data_dir         =config['data_dir']
institute        =config['institute']
model            =config['model']
original_sim_name=config['original_sim_name']
abc              =config['abc']

# load output metadata
metadata = load_config(args.config, 'metadata')

missval=9.969209968386869e36 # missing value in created netcdf files

epsfr=1. # cut sea, sea-ice or ice-shelf fractions below epsfr (in %)
         # (useful to avoid extreme values that can occur with very small fractions)

filmis=False  # keep False (may become a deprecitated option)
              # if True, will fill nans (missing values in input file or values that can't be triangulated) 
              # through nearest-neighbor interpolation. This can be useful for the edges of non-structured grids.

# dimension names
cdlon='lon'
cdlat='lat'
cdlev='lev'

#--------------------------------------------------------------------------
# 3- Interpolate to common 3d grid :
#--------------------------------------------------------------------------

for reg in tasks.keys():

    task=tasks[reg]

    # replace region placeholder
    exp =config['exp'].format(regname=task['name'])

    print('')
    print(f'---------------------')
    print(f'process region: {reg}')
    print(f'---------------------')
    print('')

    #--------------------------------------------------------------------------
    # 3a- Files and variables
    # loading an xarray dataset containing all useful variables with (x,y) reshaped 
    # as 1 dimension in case of original structured grid :
    #--------------------------------------------------------------------------
    if ( (test_case[0:4] == 'NEMO') | ( (test_case[0:5] == 'eORCA') and (test_case != 'eORCA025_prod' )) ):
    
       print('LOADING NEMO...')
    
       fs_gridT = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_T.nc' for month in np.arange(1,3)]
       fs_gridU = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_U.nc' for month in np.arange(1,3)]
       fs_gridV = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_grid_V.nc' for month in np.arange(1,3)]
       fs_ice   = [data_dir+'/'+test_case+'_m0'+month.astype('str')+'_icemod.nc' for month in np.arange(1,3)]
    
       print(fs_gridT)
    
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
                      files_BSF=fs_BSF, rho0=1026.0, teos10=True, region=reg )
    
    elif ( test_case == 'eORCA025_prod' ):
    
       yearb=args.yearb
       yeare=args.yeare
    
       print(f'Loading monthly NEMO from {yearb} to {yeare} ...')
       f_mesh   = f'{data_dir}/{original_sim_name}_mesh_mask.nc'
       f_bathy  = 'None'
    
       fs_SBC   = generate_filenames(data_dir, original_sim_name, 'flxT'  , yearb, yeare)
       fs_BSF   = generate_filenames(data_dir, original_sim_name, 'BSF'   , yearb, yeare)
       fs_gridT = generate_filenames(data_dir, original_sim_name, 'gridT' , yearb, yeare)
       fs_gridU = generate_filenames(data_dir, original_sim_name, 'gridU' , yearb, yeare)
       fs_gridV = generate_filenames(data_dir, original_sim_name, 'gridV' , yearb, yeare)
       fs_ice   = generate_filenames(data_dir, original_sim_name, 'icemod', yearb, yeare)
    
       oce = mp.load_oce_mod_nemo( file_mesh_mask=f_mesh, file_bathy=f_bathy, files_gridT=fs_gridT,\
                      files_gridU=fs_gridU, files_gridV=fs_gridV, files_SBC=fs_SBC, files_ice=fs_ice,\
                      files_BSF=fs_BSF, rho0=1026.0, teos10=True, region=reg )
    
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

    print(' Data loaded, start processing')
    print('')
    
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
    
    # Sea Area Percentage at the Surface on MISOMIP grid (in [0-100], =nan beyond model domain)
    SFTOF_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SFTOF )
    SFTOF_miso[ np.isnan(DOMMSK_miso) ] = 0.e0 # will update to missval at the end of the calculation

    # Ice shelf fraction on MISOMIP grid (in [0-100], =nan beyond model domain)
    SFTFLF_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SFTFLF )
    SFTFLF_miso[ np.isnan(DOMMSK_miso) ] = 0.e0 # will update to missval at the end of the calculation
    
    # Depth of floating ice on MISOMIP grid (ice-shelf draft)
    DEPFLF_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DEPFLF, weight=oce.SFTFLF, skipna=True, filnocvx=filmis, threshold=epsfr )
    DEPFLF_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPFLF_miso) ] = 0.e0  # will be replaced with missval later on
    
    # Ocean depth on MISOMIP grid (=nan where land or grounded ice or beyond model domain)
    DEPTHO_miso = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DEPTHO, weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
    DEPTHO_miso[ np.isnan(DOMMSK_miso) | np.isnan(DEPTHO_miso) ] = 0.e0 # will be replaced with missval later on
    
    # Sea area fraction at each vertical level on the MISOMIP grid: 
    tmp_LEVT = LEVOFT.interp(z=dep_miso) # extrapolate to avoid nan at surface but better to
    tmp_LEVU = LEVOFU.interp(z=dep_miso) # bring first level to zero when loading model data.
    tmp_LEVV = LEVOFV.interp(z=dep_miso)
    
    LEVOF_miso = np.zeros((mdep,mlat,mlon))

    # compute Delaunay tri-angulation
    # Pre-compute ONCE (outside your loop)
    print('compute triangulation (skipna = False, missval=False ...')
    triT = Delaunay(np.column_stack((lonT*aa, latT)))  # Reusable triangulation
    triV = Delaunay(np.column_stack((lonU*aa, latU)))  # Reusable triangulation
    triU = Delaunay(np.column_stack((lonV*aa, latV)))  # Reusable triangulation

    for kk in np.arange(mdep):
        print(f'{kk+1}/{mdep}')
        nanmask = (~np.isnan(tmp_LEVT[kk,:])) & (~np.isinf(tmp_LEVT[kk,:])) 
        tzz = mp.horizontal_interp_tri( triT, nanmask, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_LEVT[kk,:] )
        tzz[ np.isnan(DOMMSK_miso) ] = np.nan # will update to missval at the end of the calculation
        LEVOF_miso[kk,:,:] = tzz

    
    #----- interpolation of time-varying fields:
    
    domcond = (np.isnan(DOMMSK_miso))

    #------ write static fields:
    LEVOF_miso[ np.isnan(LEVOF_miso) ] = missval
    SFTOF_miso[ np.isnan(SFTOF_miso) | np.isnan(DOMMSK_miso) ] = missval 
    SFTFLF_miso[ np.isnan(SFTFLF_miso) | np.isnan(DOMMSK_miso) ] = missval 
    DEPFLF_miso[ (DEPFLF_miso==0.e0) | np.isnan(DOMMSK_miso) ] = missval
    DEPTHO_miso[ np.isnan(DOMMSK_miso) | (DEPTHO_miso==0.) ] = missval
 
    dsmiso3d = xr.Dataset(
        {
        "levof":     ([cdlev, cdlat, cdlon], np.float32(LEVOF_miso)),
        "sftof":     ([cdlat, cdlon], np.float32(SFTOF_miso)),
        "sftflf":    ([cdlat, cdlon], np.float32(SFTFLF_miso)),
        "depflf":    ([cdlat, cdlon], np.float32(DEPFLF_miso)),
        "deptho":    ([cdlat, cdlon], np.float32(DEPTHO_miso)),
        },
        coords={
        cdlon:np.float32(lon_miso),
        cdlat: np.float32(lat_miso),
        cdlev: np.float32(dep_miso),
        "time": oce.time.values
        },
    )
    del SFTFLF_miso, DEPFLF_miso, DEPTHO_miso
   
    #---- 2D T fields

    dsmiso3d['zos']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['tob']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['sob']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['hfs']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['flandice'] = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    for ll in np.arange(mtime):
    
        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print('2D interpolation T fields')
        print(datetime.now() - startTime)

        #-----------------
        #tri for ZOS, TOB, SOB, HFS, WFOAT
        triT, nanmaskT = compute_tri_and_mask(oce.ZOS.isel(time=ll), lonT*aa, latT, nandrop=True)
   
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.ZOS.isel(time=ll), \
                                  weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['zos'][ll,:,:] = tzz
    
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TOB.isel(time=ll), \
                                  weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['tob'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SOB.isel(time=ll), \
                                    weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['sob'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.HFS.isel(time=ll), \
                                    weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['hfs'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.FLANDICE.isel(time=ll), \
                                    weight=seafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['flandice'][ll,:,:] = tzz

    #----- U & V grids 2D
    dsmiso3d['tauvo']    = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['tauuo']    = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['msftbarot']= (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    for ll in np.arange(mtime):

        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print('2D interpolation U & V fields')
        print(datetime.now() - startTime)
      
        ## fields interpolated from sea cells ( U & V grids ):
        tzz = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.MSFTBAROT.isel(time=ll), \
                                    weight=seafracU, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['msftbarot'][ll,:,:] = tzz
      
        TAUX_notrot = mp.horizontal_interp( latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TAUX.isel(time=ll), \
                                            weight=seafracU, skipna=True, filnocvx=filmis, threshold=epsfr )

        TAUY_notrot = mp.horizontal_interp( latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.TAUY.isel(time=ll), \
                                            weight=seafracV, skipna=True, filnocvx=filmis, threshold=epsfr )

        tzz = TAUX_notrot * np.cos(theU) - TAUY_notrot * np.sin(theV) # rotated to zonal
        tzz[ ( np.isnan(tzz) ) | domcond ] = missval
        dsmiso3d['tauuo'][ll,:,:] = tzz

        tzz = TAUY_notrot * np.cos(theV) + TAUX_notrot * np.sin(theU) # rotated to meridional
        tzz[ ( np.isnan(tzz) ) | domcond ] = missval
        dsmiso3d['tauvo'][ll,:,:] = tzz

    #------- sea ice fields 2D
    dsmiso3d['wfocorr']  = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['wfoat']    = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['siconc']   = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['sivol']    = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['fsitherm'] = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['siu']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['siv']      = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))

    for ll in np.arange(mtime):

        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print('2D interpolation open-sea and sea ice fields')
        print(datetime.now() - startTime)
      
        ## fileds interpolated from open-sea cells :
        triT, nanmaskT = compute_tri_and_mask(oce.WFOCORR.isel(time=ll), lonT*aa, latT, nandrop=True)

        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.WFOCORR.isel(time=ll), \
                                    weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['wfocorr'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.WFOAT.isel(time=ll), \
                                    weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['wfoat'][ll,:,:] = tzz

        #tri for sea ice variable SICONC, SIVOL
        triT, nanmaskT = compute_tri_and_mask(oce.SICONC.isel(time=ll), lonT*aa, latT, nandrop=True)

        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SICONC.isel(time=ll), \
                                    weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = np.nan # will be updated to missval later on
        dsmiso3d['siconc'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIVOL.isel(time=ll), \
                                    weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['sivol'][ll,:,:] = tzz

        tzz = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.FSITHERM.isel(time=ll), \
                                    weight=openseafracT, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['fsitherm'][ll,:,:] = tzz

      
        ## fileds interpolated from sea-ice cells :
      
        SIUX_notrot = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIUX.isel(time=ll), \
                                            weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=filmis, threshold=epsfr )
        SIVY_notrot = mp.horizontal_interp( latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.SIVY.isel(time=ll), \
                                            weight=oce.SICONC.isel(time=ll), skipna=True, filnocvx=filmis, threshold=epsfr )

        tzz = SIUX_notrot * np.cos(theT) - SIVY_notrot * np.sin(theT) # rotated to zonal
        tzz [ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['siu'][ll,:,:] = tzz

        tzz = SIVY_notrot * np.cos(theT) + SIUX_notrot * np.sin(theT) # rotated to meridional
        tzz [ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['siv'][ll,:,:] = tzz

    ## fileds interpolated from ice-shelf cells :
    dsmiso3d['libmassbffl'] = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['dydrflf']     = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['thdrflf']     = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    dsmiso3d['hadrflf']     = (["time", cdlat, cdlon], np.float32(np.zeros((mtime,mlat,mlon)) + missval))
    for ll in np.arange(mtime):

        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print('2D interpolation open-sea and sea ice fields')
        print(datetime.now() - startTime)
      
        #tri for ice shelf variable LIBMASSBFFL, DYDRFLF, THDRFLF, HADRFLF
        triT, nanmaskT = compute_tri_and_mask(oce.LIBMASSBFFL.isel(time=ll), lonT*aa, latT, nandrop=True)
        
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.LIBMASSBFFL.isel(time=ll), \
                                    weight=oce.SFTFLF, skipna=True, filnocvx=filmis, threshold=epsfr )
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['libmassbffl'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.DYDRFLF.isel(time=ll), \
                                    weight=oce.SFTFLF, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['dydrflf'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri(  triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.THDRFLF.isel(time=ll), \
                                    weight=oce.SFTFLF, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['thdrflf'][ll,:,:] = tzz
      
        tzz = mp.horizontal_interp_tri(  triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, oce.HADRFLF.isel(time=ll), \
                                    weight=oce.SFTFLF, skipna=True, filnocvx=filmis, threshold=epsfr ) 
        tzz[ np.isnan(tzz) | domcond ] = missval
        dsmiso3d['hadrflf'][ll,:,:] = tzz

    print('3D interpolation T fields')
    print(datetime.now() - startTime)
     # vertical interpolation of T, S, U, V to common vertical grid :
    
    tmpxS = LEVOFT*oce.SO
    tmp_SS = tmpxS.interp(z=dep_miso) / tmp_LEVT
    del tmpxS
    
    dsmiso3d['so'] = (["time", cdlev, cdlat, cdlon], np.float32(np.zeros((mtime,mdep,mlat,mlon)) + missval))
    for ll in np.arange(mtime):

        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print(datetime.now() - startTime)

        for kk in np.arange(mdep):
      
            # horizontal interpolations of time-varying 3d fields to common horizontal grid : 
            print(f'{kk+1}/{mdep}')
            
            triT, nanmaskT = compute_tri_and_mask(tmp_SS.isel(time=ll,z=kk), lonT*aa, latT, nandrop=True)

            condkk = ( (LEVOF_miso[kk,:,:] < epsfr) | domcond )
            
            tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_SS.isel(time=ll,z=kk), \
                                        weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
            tzz[ condkk | (np.isnan(tzz)) ] = missval
            dsmiso3d['so'][ll,kk,:,:] = tzz
    del tmp_SS
    
    #--------- 3D temperature
    tmpxT = LEVOFT*oce.THETAO
    tmp_TT = tmpxT.interp(z=dep_miso) / tmp_LEVT
    del tmpxT

    dsmiso3d['thetao']   = (["time", cdlev, cdlat, cdlon], np.float32(np.zeros((mtime,mdep,mlat,mlon)) + missval))

    for ll in np.arange(mtime):

        ## fields interpolated from sea cells ( C/T grid ):
        print(f'{ll+1}/{mtime}')
        print(datetime.now() - startTime)

        for kk in np.arange(mdep):

            # horizontal interpolations of time-varying 3d fields to common horizontal grid : 
            print(f'{kk+1}/{mdep}')

            triT, nanmaskT = compute_tri_and_mask(tmp_TT.isel(time=ll,z=kk), lonT*aa, latT, nandrop=True)

            condkk = ( (LEVOF_miso[kk,:,:] < epsfr) | domcond )

            tzz = mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_TT.isel(time=ll,z=kk), \
                                        weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )
            tzz[ condkk | (np.isnan(tzz)) ] = missval
            dsmiso3d['thetao'][ll,kk,:,:] = tzz
    del tmp_TT

    #--------- velocity
    tmpxU = LEVOFU*oce.UX
    tmp_UX = tmpxU.interp(z=dep_miso) / tmp_LEVU
    del tmpxU
    dsmiso3d['uo']   = (["time", cdlev, cdlat, cdlon], np.float32(np.zeros((mtime,mdep,mlat,mlon)) + missval))
    
    tmpxV = LEVOFV*oce.VY
    tmp_VY = tmpxV.interp(z=dep_miso) / tmp_LEVV
    del tmpxV
    dsmiso3d['vo']   = (["time", cdlev, cdlat, cdlon], np.float32(np.zeros((mtime,mdep,mlat,mlon)) + missval))

    print('3D interpolation UV fields')
    print(datetime.now() - startTime)
    for ll in np.arange(mtime):

        print(f'{ll+1}/{mtime}')
        print(datetime.now() - startTime)

        for kk in np.arange(mdep):

            # horizontal interpolations of time-varying 3d fields to common horizontal grid :
            print(f'{kk+1}/{mdep}')

            triU, nanmaskU = compute_tri_and_mask(tmp_UX.isel(time=ll,z=kk), lonU*aa, latU, nandrop=True)
            triV, nanmaskV = compute_tri_and_mask(tmp_VY.isel(time=ll,z=kk), lonV*aa, latU, nandrop=True)

            condkk = ( (LEVOF_miso[kk,:,:] < epsfr) | domcond )
        
            UX_notrot = mp.horizontal_interp_tri( triU, nanmaskU, latU, lonU*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_UX.isel(time=ll,z=kk), \
                                              weight=tmp_LEVU.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )

            VY_notrot = mp.horizontal_interp_tri( triV, nanmaskV, latV, lonV*aa, mlat, mlon, lat_miso1d, lon_miso1d*aa, tmp_VY.isel(time=ll,z=kk), \
                                              weight=tmp_LEVV.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr )

            tzz = UX_notrot * np.cos(theU) - VY_notrot * np.sin(theV) # rotated to zonal
            tzz[ condkk | (np.isnan(tzz)) ] = missval
            dsmiso3d['uo'][ll,kk,:,:] = tzz

            tzz = VY_notrot * np.cos(theV) + UX_notrot * np.sin(theU) # rotated to meridional
            tzz[ condkk | (np.isnan(tzz)) ] = missval
            dsmiso3d['vo'][ll,kk,:,:] = tzz
    

    del tmp_UX, tmp_VY
   
    #--------------------------------------------------------------------------
    # 3c- Create new xarray dataset and save to netcdf
    #--------------------------------------------------------------------------
    
    mp.add_standard_attributes_oce(dsmiso3d,miss=missval)
    put_global_attrs(dsmiso3d, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name, \
                     ocean_model=model, institute=institute, \
                     original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                     original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon'),\
                     metadata=metadata)
    
    file_miso3d = 'Oce3d_'+institute+'_'+model+'_'+abc+'_'+exp+'_'+period+'.nc'
    
    print('Creating ',file_miso3d)
    dsmiso3d.to_netcdf(file_miso3d,unlimited_dims="time")
    
    del dsmiso3d
    
    print('  Execution time: ',datetime.now() - startTime)

    #--------------------------------------------------------------------------
    # 4a- Interpolate to common section :
    for sec in task['sections']:
        print('')
        print(f'    ---------------------')
        print(f'    process section: {sec}')
        print(f'    ---------------------')
        print('')

        # Characteristics of MISOMIP section grid:
        [lon_sect1d,lat_sect1d,dep_sect] = mp.generate_section_grid_oce(region=reg,section=sec)
        mlonlatsec = np.size(lon_sect1d)
        mdepsect = np.size(dep_sect)

        # Reduce input data size
        #   NB: for some reason, a small number of points (low wdeg below) may lead to the following error here:
        #       scipy.spatial.qhull.QhullError: QH6214 qhull input error: not enough points(1) to construct initial simplex (need 4) 
        #       If it happens with your dataset, increase it further...
        wdeg = 3.0 * res_73S / ( 6.37e6 * np.pi / 180. )
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
        
        # vertical then horizontal interpolation of constant 3d fields to common grid :
        tmp_LEVT = LEVOFT.interp(z=dep_sect)
        
        LEVOF_sect = np.zeros((mdepsect,mlonlatsec))

        print(datetime.now())
        # as skipnan is false, same tri for all levels
        tri = Delaunay(np.column_stack((lonT*aa, latT)))
        print('LEVT interpolation')
        for kk in np.arange(mdepsect):
            if (np.mod(kk+1,100) == 0) :
                print(f'  {kk+1}/{mdepsect}')
            nanmask = (~np.isnan(tmp_LEVT[kk,:])) & (~np.isinf(tmp_LEVT[kk,:]))
            tzz = np.squeeze( mp.horizontal_interp_tri( tri, nanmask, latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_LEVT[kk,:] ) )
            tzz[ domcond ] = missval
            LEVOF_sect[kk,:] = tzz
        print(datetime.now())
        
        dssect = xr.Dataset(
            {
            "levof":     ([cdlev, "x"], np.float32(LEVOF_sect)),
            cdlon: (["x"], np.float32(lon_sect1d)),
            cdlat: (["x"], np.float32(lat_sect1d)),
            },
            coords={
            cdlev: np.float32(dep_sect),
            "time": oce.time.values
            },
        )

        # vertical interpolation of T, S, U, V to common vertical grid :
        tmpxS = LEVOFT*oce.SO.where(cond_secT,drop=True)
        tmp_SS = tmpxS.interp(z=dep_sect) / tmp_LEVT
        del tmpxS
        dssect['so']     = (["time", cdlev, "x"], np.zeros((mtime,mdepsect,mlonlatsec)) + missval)

        # horizontal interpolation
        print('Horizontal interpolation')
        for ll in np.arange(mtime):
            print(f'  {ll}/{mtime}')
            print(datetime.now())
            for kk in np.arange(mdepsect):
              if (np.mod(kk,100) == 0) :
                  print(f'    {kk+1}/{mdepsect}')

              triT, nanmaskT = compute_tri_and_mask(tmp_SS.isel(time=ll,z=kk), lonT*aa, latT, nandrop=True)

              condkk = ( (LEVOF_sect[kk,:] < epsfr) | domcond )
          
              tzz = np.squeeze(mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_SS.isel(time=ll,z=kk), \
                                                     weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr ) )
              tzz[ condkk | (np.isnan(tzz)) ] = missval
              dssect['so'][ll,kk,:] = tzz
        del tmp_SS

        # vertical interpolation of T, S, U, V to common vertical grid :
        tmpxT = LEVOFT*oce.THETAO.where(cond_secT,drop=True)
        tmp_TT = tmpxT.interp(z=dep_sect) / tmp_LEVT
        del tmpxT
        dssect['thetao'] = (["time", cdlev, "x"], np.zeros((mtime,mdepsect,mlonlatsec)) + missval)
        
        # horizontal interpolation
        for ll in np.arange(mtime):
            print(f'  {ll}/{mtime}')
            print(datetime.now())
            for kk in np.arange(mdepsect):
              if (np.mod(kk,100) == 0) :
                  print(f'    {kk+1}/{mdepsect}')

              triT, nanmaskT = compute_tri_and_mask(tmp_TT.isel(time=ll,z=kk), lonT*aa, latT, nandrop=True)

              condkk = ( (LEVOF_sect[kk,:] < epsfr) | domcond )

          
              tzz = np.squeeze(mp.horizontal_interp_tri( triT, nanmaskT, latT, lonT*aa, mlonlatsec, 1, lat_sect1d, lon_sect1d*aa, tmp_TT.isel(time=ll,z=kk), \
                                                     weight=tmp_LEVT.isel(z=kk), skipna=True, filnocvx=filmis, threshold=epsfr ) )
              tzz[ condkk | (np.isnan(tzz)) ] = missval
              dssect['thetao'][ll,kk,:] = tzz
        del tmp_TT

        print(datetime.now())
        
        #--------------------------------------------------------------------------
        # 4b- Create new xarray dataset and save to netcdf
        mp.add_standard_attributes_oce(dssect,miss=missval)
        put_global_attrs(dssect, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name,\
                         ocean_model=model, institute=institute, \
                         original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                         original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon'),\
                         metadata=metadata)
        
        file_sect = 'OceSec'+str(sec)+'_'+institute+'_'+model+'_'+abc+'_'+exp+'_'+period+'.nc'
        
        print('Creating ',file_sect)
        dssect.to_netcdf(file_sect,unlimited_dims="time")
        
        del dssect
        
        print('  Execution time: ',datetime.now() - startTime)
       
    #--------------------------------------------------------------------------
    # 5a- Interpolate to common mooring location :
    #
    #     For this mooring, we take the closest grid point that is not under 
    #     an ice shelf (observational mooring in front of a real ice shelf
    #     can be located in the cavity of some models).
    for moor in task['moorings']:
        
        print('')
        print(f'    ---------------------')
        print(f'    process mooring: {moor}')
        print(f'    ---------------------')
        print('')

        # Characteristics of MISOMIP mooring:
        [lon_moor0d,lat_moor0d,dep_moor] = mp.generate_mooring_grid_oce(region=reg,mooring=moor)
        mdepmoor = np.size(dep_moor)
        
        # Reduce input data size
        #   NB: for some reason, a small number of points (low wdeg below) may lead to the following error here:
        #       scipy.spatial.qhull.QhullError: QH6214 qhull input error: not enough points(1) to construct initial simplex (need 4) 
        #       If it happens with your dataset, increase it further...
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

        fracT = seafracT

        mtime = np.shape(oce.SO)[0]
        
        # Coordinates of closest point:
        condclo = ((lat_moor0d-latT.where(fracT >= 50.))**2+((lon_moor0d-lonT.where(fracT >= 50.))*aa)**2)
        is_all_nan = condclo.isnull().all().item()

        if is_all_nan:

            print('No points found for mooring',moor,' at (lon,lat) = ',lon_moor0d,lat_moor0d)
            print('   >>>>>>> file not saved')

        else:

            indxy = ((lat_moor0d-latT.where(fracT >= 50.))**2+((lon_moor0d-lonT.where(fracT >= 50.))*aa)**2).argmin()
            
            print('Closest point of open sea near mooring of (lon,lat) = ',lon_moor0d,lat_moor0d)
            print('   is at model point of (lon,lat) = ',lonT.where(fracT >= 50.).isel(sxy=indxy).values, latT.where(fracT >= 50.).isel(sxy=indxy).values )
            DOMMSK_moor = oce.DOMMSKT.where(cond_mooT,drop=True).where(fracT >= 50.).isel(sxy=indxy).values
            
            tmp_LEVT = LEVOFT.interp(z=dep_moor)
            LEVOF_moor = np.zeros((mdepmoor))
            for kk in np.arange(mdepmoor):
                tzz = tmp_LEVT[kk,:].where(fracT >= 50.).isel(sxy=indxy).values
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
            
                    tzz = tmp_SS[kk,:,ll].where(fracT >= 50.).isel(sxy=indxy).values
                    tzz[ np.isnan(DOMMSK_moor) | (LEVOF_moor[kk] < epsfr) ] = missval
                    SO_moor[ll,kk] = tzz
            
                    tzz = tmp_TT[kk,:,ll].where(fracT >= 50.).isel(sxy=indxy).values
                    tzz[ np.isnan(DOMMSK_moor) | (LEVOF_moor[kk] < epsfr) ] = missval
                    THETAO_moor[ll,kk] = tzz
            
            #--------------------------------------------------------------------------
            # 5b- Create new xarray dataset and save to netcdf
            
            dsmoor = xr.Dataset(
                {
                "so":        (["time", cdlev], np.float32(SO_moor)),
                "thetao":    (["time", cdlev], np.float32(THETAO_moor)),
                "levof":     ([cdlev], np.float32(LEVOF_moor)),
                },
                coords={
                cdlev: np.float32(dep_moor),
                "time": oce.time.values
                },
            )
            
            mp.add_standard_attributes_oce(dsmoor,miss=missval)
            put_global_attrs(dsmoor, experiment=exp, avg_hor_res_73S=res_73S, original_sim_name=original_sim_name,\
                            ocean_model=model, institute=institute, \
                            original_min_lat=oce.attrs.get('original_minlat'),original_max_lat=oce.attrs.get('original_maxlat'),\
                            original_min_lon=oce.attrs.get('original_minlon'),original_max_lon=oce.attrs.get('original_maxlon'),\
                            metadata=metadata)
            dsmoor.attrs['mooring_longitude'] = np.float32(lon_moor0d)
            dsmoor.attrs['mooring_latitude'] = np.float32(lat_moor0d)
            
            file_moor = 'OceMoor'+str(moor)+'_'+institute+'_'+model+'_'+abc+'_'+exp+'_'+period+'.nc'
            
            print('Creating ',file_moor)
            dsmoor.to_netcdf(file_moor,unlimited_dims="time")
            print('  Execution time: ',datetime.now() - startTime)
            
    print('  Execution time: ',datetime.now() - startTime)
