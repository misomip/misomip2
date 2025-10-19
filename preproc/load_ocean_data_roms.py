# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#====================================================================================

import numpy as np
import xarray as xr
import gsw
from .def_grids import grid_bounds_oce
from .interp_functions import calc_z
from datetime import datetime

#====================================================================================
def load_oce_mod_roms(files_T='ROMS_all.nc',\
                      files_S='dummy',\
                      files_U='dummy',\
                      files_V='dummy',\
                      files_I='dummy',\
                      files_SRF='dummy',\
                      files_M='dummy',\
                      rho0=1026.0, teos10=False, region='Amundsen', parallel=False ):
   """ Read ROMS outputs and define an xarray dataset containing 
       all variables required in MISOMIP2. It automatically detects
       whether coordinates are stereographic or lon-lat.

       Input:

       files_T: file or list of files containing the temperature and related variables [default='ROMS_all.nc']
       files_S: file or list of files containing the salinity variable [default=files_T]
       files_U: file or list of files containing the x-velocity and related variables [default=files_T]
       files_V: file or list of files containing the y-velocity and related variables [default=files_T]
       files_I: file or list of files containing the sea-ice variables [default=files_T]
       files_SRF: file or list of files containing the surface fluxes variables [default=files_T]
       files_M: file or list of files containing grid/mesh variables [default=files_T]

       rho0: volumic mass of seawater used in ocean model

       teos10=False -> assumes the nemo outputs are in potential temperature & practical salinity (EOS80)

             =True  -> assumes the nemo outputs are in CT and AS and convert to PT and PS
        
       parallel: If True, the open and preprocess steps of this function will be performed in parallel

       Output:
          xarray dataset of coordinates ("time", "z", "sxy") (sxy= one-dimensionalized horizontal space)

       Example1:
          ds = load_oce_mod_roms()

       Example2:
          dir= 'datadir/model/'
          ff = [ dir+'ROMS_y2009.nc', dir+'ROMS_y2010.nc', dir+'ROMS_y2011.nc' ]
          ds = load_oce_mod_roms(files_T=ff, files_M='grid.nc', rho0=1028.0, region='Weddell')

   """

   startTime = datetime.now()

   # if only files_T are specified, use it for all variables:
   if ( files_S == 'dummy' ):
      files_S = files_T
   if ( files_U == 'dummy' ):
      files_U = files_T
   if ( files_V == 'dummy' ):
      files_V = files_T
   if ( files_I == 'dummy' ):
      files_I = files_T
   if ( files_SRF == 'dummy' ):
      files_SRF = files_T
   if ( files_M == 'dummy' ):
      files_M = files_T

   ncT = xr.open_mfdataset(files_T, decode_coords=False, parallel=parallel) ; ncT=ncT.rename({'ocean_time':'time'})
   ncS = xr.open_mfdataset(files_S, decode_coords=False, parallel=parallel) ; ncS=ncS.rename({'ocean_time':'time'})
   ncU = xr.open_mfdataset(files_U, decode_coords=False, parallel=parallel) ; ncU=ncU.rename({'ocean_time':'time'})
   ncV = xr.open_mfdataset(files_V, decode_coords=False, parallel=parallel) ; ncV=ncV.rename({'ocean_time':'time'})
   ncI = xr.open_mfdataset(files_I, decode_coords=False, parallel=parallel) ; ncI=ncI.rename({'ocean_time':'time'})
   ncSRF = xr.open_mfdataset(files_SRF, decode_coords=False, parallel=parallel) ; ncSRF=ncSRF.rename({'ocean_time':'time'})
   ncM = xr.open_mfdataset(files_M, decode_coords=False, parallel=parallel)

   [mtime,ms,my,mx] = ncT.temp.shape

   # fill all nans with a non-nan value (useful to keep points in triangular interpolation, even if they will be masked)
   ncT=ncT.fillna(-99999.99)
   ncS=ncS.fillna(-99999.99)
   ncU=ncU.fillna(-99999.99)
   ncV=ncV.fillna(-99999.99)
   ncI=ncI.fillna(-99999.99)
   ncSRF=ncSRF.fillna(-99999.99)
   ncM=ncM.fillna(-99999.99)

   # longitude & latitude on U, V, T grids
   lonT = ncM.lon_rho
   latT = ncM.lat_rho
   lonU = ncM.lon_u
   latU = ncM.lat_u
   lonV = ncM.lon_v
   latV = ncM.lat_v

   # save original domain boundaries:
   domain_minlat = latT.min().values
   domain_maxlat = latT.max().values
   domain_minlon = lonT.min().values
   domain_maxlon = lonT.max().values

   # grid mesh widths along x and y at C/T, U and V points [m]:
   dxT = 1./ncM.pm
   dyT = 1./ncM.pn

   # local C/T, U, V grid rotation angle compared to the (zonal,meridional) direction [rad]
   thetaT = ncM.angle # angle between xi_rho and EAST
   thetaU = xr.DataArray( 0.500000000 * (ncM.angle.values+ncM.angle.shift(xi_rho=1).values), dims=['eta_u', 'xi_u'] )
   thetaV = xr.DataArray( 0.500000000 * (ncM.angle.values+ncM.angle.shift(eta_rho=1).values), dims=['eta_v', 'xi_v'] )
   print('    Minimum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.min().values*180./np.pi)
   print('    Maximum local grid angle in degrees w.r.t. (zonal,meridional):',thetaU.max().values*180./np.pi)

   # Domain mask (ones with a halo of nans), used not to interpolate beyond the original domain:
   halonan = np.ones((my,mx))
   halonan[0,:] = np.nan ; halonan[-1,:] = np.nan
   halonan[:,0] = np.nan ; halonan[:,-1] = np.nan
   DOMMSKT = xr.DataArray( halonan, dims=['eta_rho', 'xi_rho'] )
   DOMMSKU = xr.DataArray( halonan, dims=['eta_u', 'xi_u'] )
   DOMMSKV = xr.DataArray( halonan, dims=['eta_v', 'xi_v'] )

   # ice-shelf fraction [0-100]:
   SFTFLF = ncM.mask_rho.where( ncM.zice<-1., 0.0 )*100.

   # Bathymetry (including under ice shelves) [m, positive]
   if ( "h" in ncM.data_vars ):
     DEPTHO = ncM.h
   else:
     print('    WARNING :   No data found for DEPTHO  -->  filled with NaNs')
     DEPTHO = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['eta_rho', 'xi_rho'] ) 
 
   # Depth of ice shelf draft [m, positive]:
   if ( "zice" in ncM.data_vars ):
     DEPFLF = ncM.zice*(-1)
     DEPFLF = DEPFLF.where( ncM.mask_rho > 0.5, 0.e0)
   else:
     print('    WARNING :   No data found for DEPFLF  -->  filled with NaNs')
     DEPFLF = xr.DataArray( np.zeros((my,mx))*np.nan, dims=['eta_rho', 'xi_rho'] )
 
   # ocean temperature [degC]
   if ( "temp" in ncT.data_vars ):
     TT = ncT.temp
   elif ( "thetao" in ncT.data_vars ):
     TT = ncT.thetao
   else:
     print('    WARNING :   No data found for TT  -->  filled with NaNs')
     TT = xr.DataArray( np.zeros((mtime,ms,my,mx))*np.nan, dims=['time', 's_rho', 'eta_rho', 'xi_rho'] )

   # ocean salinity [1.e-3]
   if ( "salt" in ncS.data_vars ):
     SS = ncS.salt
   elif ( "so" in ncS.data_vars ):
     SS = ncS.so
   else:
     print('    WARNING :   No data found for SS  -->  filled with NaNs')
     SS = xr.DataArray( np.zeros((mtime,ms,my,mx))*np.nan, dims=['time', 's_rho', 'eta_rho', 'xi_rho'] )

   # sea bottom ocean temperature [degC]
   if ( "temp" in ncT.data_vars ):
     TTB = ncT.temp.isel(s_rho=ms-1) # deepest level
   elif ( "thetao" in ncT.data_vars ):
     TTB = ncT.thetao.isel(s_rho=ms-1)
   else:
     print('    WARNING :   No data found for TTB  -->  filled with NaNs')
     TTB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # sea bottom ocean salinity [1.e-3]
   if ( "salt" in ncT.data_vars ):
     SSB = ncT.salt.isel(s_rho=ms-1) # deepest level
   elif ( "so" in ncT.data_vars ):
     SSB = ncT.so.isel(s_rho=ms-1)
   else:
     print('    WARNING :   No data found for SSB  -->  filled with NaNs')
     SSB = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

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
   if ( "u" in ncU.data_vars ):
     UX = ncU.u
   elif ( "uo" in ncU.data_vars ):
     UX = ncU.uo
   else:
     print('    WARNING :   No data found for UX  -->  filled with NaNs')
     UX = xr.DataArray( np.zeros((mtime,ms,my,mx))*np.nan, dims=['time', 's_rho', 'eta_u', 'xi_u'] )

   # ocean y-ward velocity [m s-1]
   if ( "v" in ncV.data_vars ):
     VY = ncV.v
   elif ( "vo" in ncV.data_vars ):
     VY = ncV.vo
   else:
     print('    WARNING :   No data found for VY  -->  filled with NaNs')
     VY = xr.DataArray( np.zeros((mtime,ms,my,mx))*np.nan, dims=['time', 's_rho', 'eta_v', 'xi_v'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "utau" in ncU.data_vars ):
     TAUX = ncU.utau
   elif ( "TAUX" in ncU.data_vars ):
     TAUX = ncU.TAUX
   else:
     print('    WARNING :   No data found for TAUX  -->  filled with NaNs')
     TAUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_u', 'xi_u'] )

   # surface stress received by the ocean along x [W m-1]
   if ( "vtau" in ncV.data_vars ):
     TAUY = ncV.vtau
   elif ( "TAUY" in ncV.data_vars ):
     TAUY = ncV.TAUY
   else:
     print('    WARNING :   No data found for TAUY  -->  filled with NaNs')
     TAUY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_v', 'xi_v'] )

   # Sea surface height [m]
   if ( "zos" in ncT.data_vars ):
     ZOS = ncT.zos
   elif ( "zeta" in ncT.data_vars ):
     ZOS = ncT.zeta
   elif ( "ssh" in ncT.data_vars ):
     ZOS = ncT.ssh
   else:
     print('    WARNING :   No data found for ZOS  -->  filled with NaNs')
     ZOS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # mass barotropic streamfunction
   # see Griffies et al. (2016, section H26): d(psi)/dy=-U (U: x-ward mass transport), d(psi)/dx=V (V: yward mass transport)
   if ( "sobarstf" in ncU.data_vars ):
     MSFTBAROT = ncU.sobarstf * rho0
   else:
     print('    WARNING :   No data found for MSFTBAROT  -->  filled with NaNs')
     MSFTBAROT = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_u', 'xi_u'] )

   # ice shelf melt [kg m-2 s-1, positive for actual melting] :
   if ( "m" in ncT.data_vars ):
     FICESHELF = ncT.m*920.0
   else:
     print('    WARNING :   No data found for FICESHELF  -->  filled with NaNs')
     FICESHELF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # ice shelf dynamical driving (heat exchange velocity) [m s-1]:
   if ( "isfgammat" in ncT.data_vars ):
     DYDRFLF = ncT.isfgammat
   else:
     print('    WARNING :   No data found for DYDRFLF  -->  filled with NaNs')
     DYDRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # ice shelf thermal driving [degC]:
   if ( "isfthermdr" in ncT.data_vars ):
     THDRFLF = ncT.isfthermdr
   else:
     print('    WARNING :   No data found for THDRFLF  -->  filled with NaNs')
     THDRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # ice shelf haline driving [0.001]:
   if ( "isfhalindr" in ncT.data_vars ):
     HADRFLF = ncT.isfhalindr
   else:
     print('    WARNING :   No data found for HADRFLF  -->  filled with NaNs')
     HADRFLF = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # sea-ice concentration [0-100]
   if ( "siconc" in ncI.data_vars ):
     SICONC = ncI.siconc*100.0
     SICONC = SICONC.where( (~np.isnan(SICONC)) & (~np.isinf(SICONC)), 0.e0 )
   else:
     print('    WARNING :   No data found for SICONC  -->  filled with NaNs')
     SICONC = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )   

   # sea-ice volume per area [m]
   if ( "sivolu" in ncI.data_vars ):
     SIVOL = ncI.sivolu
   elif ( "sivol" in ncI.data_vars ):
     SIVOL = ncI.sivol
   else:
     print('    WARNING :   No data found for SIVOL  -->  filled with NaNs')
     SIVOL = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # sea-ice x-ward velocity [m/s]
   if ( "sivelu" in ncI.data_vars ):
     SIUX = ncI.sivelu
   elif ("siu" in ncI.data_vars ):
     SIUX = ncI.siu
   else:
     print('    WARNING :   No data found for SIUX  -->  filled with NaNs')
     SIUX = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # sea-ice y-ward velocity [m/s]
   if ( "sivelv" in ncI.data_vars ):
     SIVY = ncI.sivelv
   elif ("siv" in ncI.data_vars ):
     SIVY = ncI.siv
   else:
     print('    WARNING :   No data found for SIUY  -->  filled with NaNs')
     SIVY = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # Total heat flux received by the ocean surface (including ice-shelf/ocean interface) [W m-2] 
   # see Griffies et al. (2016, section K4-K5) NB: here, including correction if any unlike Griffies (to avoid 2 variables)
   if ( "qt_oce" in ncSRF.data_vars ):
     HFDS = ncSRF.qt_oce
   else:
     print('    WARNING :   No data found for HFDS  -->  filled with NaNs')
     HFDS = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # Water flux entering the ocean due to sea-ice (melting-freezing) and surface correction (SSS restoring)
   # (= fsitherm + wfocorr in Griffies 2016 section K2) [kg m-2 s-1]
   if ( "wfocorr" in ncSRF.data_vars ):
     WFOCORR = - ncSRF.wfocorr
   else:
     WFOCORR = xr.DataArray( np.zeros((mtime,my,mx)), dims=['time', 'eta_rho', 'xi_rho'] )
   if ( "fsitherm" in ncSRF.data_vars ):
     WFOSICOR = WFOCORR - ncSRF.fsitherm
   else:
     print('    WARNING :   No data found for WFOSICOR  -->  filled with NaNs')
     WFOSICOR = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )

   # Water flux entering the ocean due to rainfall, snowfall, condensation - evap, 
   # river runoff, iceberg and ice-shelf melt [kg m-2 s-1]  
   # (= pr+prs+evs+ficeberg+friver+ficeshelf in Griffies 2016, section K2)
   if ( "empmr" in ncSRF.data_vars ):
     WFOATRLI = - ncSRF.empmr + FICESHELF
   else:
     print('    WARNING :   No data found for WFOATRLI  -->  filled with NaNs')
     WFOATRLI = xr.DataArray( np.zeros((mtime,my,mx))*np.nan, dims=['time', 'eta_rho', 'xi_rho'] )
  
   #----------
   # Reduce the size of ocean dataset
 
   [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region=region)
   wdeg = 1.5 * np.max([dxT.where(latT <= latmax).max(),dyT.where(latT <= latmax).max()]) / 6.37e6 * 180. / np.pi / np.cos(0.5*(latmin+latmax)*np.pi/180.)
   lonmin=lonmin - wdeg # take a bit more for interpolation
   lonmax=lonmax + wdeg
   latmin=latmin - wdeg
   latmax=latmax + wdeg

   condT2d = ( (latT >= latmin) & (latT <= latmax) & (lonT >= lonmin) & (lonT <= lonmax) )

   for ii in np.arange(latT.shape[1]):
      if ( np.sum(condT2d.isel(xi_rho=ii).values) == 0 ):
        imin=ii
      else:
        imin=ii
        break
   for ii in np.arange(latT.shape[1]-1,0,-1):
      if ( np.sum(condT2d.isel(xi_rho=ii).values) == 0 ):
        imax=ii
      else:
        imax=ii
        break
   for jj in np.arange(latT.shape[0]):
      if ( np.sum(condT2d.isel(eta_rho=jj).values) == 0 ):
        jmin=jj
      else:
        jmin=jj
        break
   for jj in np.arange(latT.shape[0]-1,0,-1):
      if ( np.sum(condT2d.isel(eta_rho=jj).values) == 0 ):
        jmax=jj
      else:
        jmax=jj
        break

   print('    Reducing domain size to useful area, i.e.: ',[imin,imax,jmin,jmax])

   #----------
   # Interpolate 3d fileds from sigma to Z

   print('    Interpolating to Z levels...')

   newdepth = np.append( np.arange(0.,70.,10.), np.arange(80.,1520.,20.) )

   mz = np.size(newdepth)

   ztmp = calc_z( ncM.h.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, \
                  ncM.zice.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, \
                  ncT.theta_s[0].values, ncT.theta_b[0].values, ncT.hc[0].values, ms, \
                  zeta=None, Vstretching=ncT.Vstretching[0].values )
   ztmp=ztmp*(-1)

   # sea fraction in Z space:
   DepTUV = xr.DataArray( newdepth, dims=['z'], coords=dict( z=(['z'], newdepth) ) )
   DEPTHOT = DEPTHO.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1))
   DEPFLFT = DEPFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1))
   DEPTHOU = xr.DataArray( 0.5 * (DEPTHOT.values+DEPTHOT.shift(xi_rho=1).values), dims=['eta_u', 'xi_u'] )
   DEPFLFU = xr.DataArray( 0.5 * (DEPFLFT.values+DEPFLFT.shift(xi_rho=1).values), dims=['eta_u', 'xi_u'] )
   DEPTHOV = xr.DataArray( 0.5 * (DEPTHOT.values+DEPTHOT.shift(eta_rho=1).values), dims=['eta_v', 'xi_v'] )
   DEPFLFV = xr.DataArray( 0.5 * (DEPFLFT.values+DEPFLFT.shift(eta_rho=1).values), dims=['eta_v', 'xi_v'] )

   LEVOFT = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'eta_rho', 'xi_rho'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFT = LEVOFT.where( ( (DepTUV==0.) & (DEPFLFT==0.) ) | ( (DEPTHOT>DepTUV) & (DEPFLFT<DepTUV) ), 0.e0 ) \
            * ncM.mask_rho.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1))
   LEVOFU = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'eta_u', 'xi_u'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFU = LEVOFU.where( ( (DepTUV==0.) & (DEPFLFU==0.) ) | ( (DEPTHOU>DepTUV) & (DEPFLFU<DepTUV) ), 0.e0 ) \
            * ncM.mask_u.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1))
   LEVOFV = xr.DataArray( np.ones((mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))*100., dims=['z', 'eta_v', 'xi_v'], coords=dict( z=(['z'], newdepth) ) )
   LEVOFV = LEVOFV.where( ( (DepTUV==0.) & (DEPFLFV==0.) ) | ( (DEPTHOV>DepTUV) & (DEPFLFV<DepTUV) ), 0.e0 ) \
            * ncM.mask_v.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1))

   ZMODT = xr.DataArray( ztmp, dims=['s_rho', 'eta_rho', 'xi_rho'] )
   ZMODU = xr.DataArray( ztmp, dims=['s_rho', 'eta_u', 'xi_u'] )
   ZMODV = xr.DataArray( ztmp, dims=['s_rho', 'eta_v', 'xi_v'] )
   SO_red = SO.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1))
   THETAO_red = THETAO.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1))
   UX_red = UX.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1))
   VY_red = VY.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1))
   SO_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   THETAO_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   UX_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))
   VY_z = np.zeros((mtime,mz,np.shape(ztmp)[1],np.shape(ztmp)[2]))

   for kk in np.arange(np.size(newdepth)):

     tmpaT = ZMODT - newdepth[kk]
     tmpaT = tmpaT.where( tmpaT >= 0., 1.e20 )
     msksupT = tmpaT.where( tmpaT==tmpaT.min('s_rho'), 0.e0 )
     msksupT = msksupT.where( msksupT==0.e0, 1.e0 )
     dzsupT = tmpaT.where( tmpaT==tmpaT.min('s_rho'), 0.e0 ).sum('s_rho')
     tmpbT = ZMODT - newdepth[kk]
     tmpbT = tmpbT.where( tmpbT < 0., -1.e20 )
     mskinfT = tmpbT.where( tmpbT==tmpbT.max('s_rho'), 0.e0 )
     mskinfT = mskinfT.where( mskinfT==0.e0, 1.e0 )
     dzinfT = -tmpbT.where( tmpbT==tmpbT.max('s_rho'), 0.e0 ).sum('s_rho')

     tmpaU = ZMODU - newdepth[kk]
     tmpaU = tmpaU.where( tmpaU >= 0., 1.e20 )
     msksupU = tmpaU.where( tmpaU==tmpaU.min('s_rho'), 0.e0 )
     msksupU = msksupU.where( msksupU==0.e0, 1.e0 )
     dzsupU = tmpaU.where( tmpaU==tmpaU.min('s_rho'), 0.e0 ).sum('s_rho')
     tmpbU = ZMODU - newdepth[kk]
     tmpbU = tmpbU.where( tmpbU < 0., -1.e20 )
     mskinfU = tmpbU.where( tmpbU==tmpbU.max('s_rho'), 0.e0 )
     mskinfU = mskinfU.where( mskinfU==0.e0, 1.e0 )
     dzinfU = -tmpbU.where( tmpbU==tmpbU.max('s_rho'), 0.e0 ).sum('s_rho')
     
     tmpaV = ZMODV - newdepth[kk]
     tmpaV = tmpaV.where( tmpaV >= 0., 1.e20 )
     msksupV = tmpaV.where( tmpaV==tmpaV.min('s_rho'), 0.e0 )
     msksupV = msksupV.where( msksupV==0.e0, 1.e0 )
     dzsupV = tmpaV.where( tmpaV==tmpaV.min('s_rho'), 0.e0 ).sum('s_rho')
     tmpbV = ZMODV - newdepth[kk]
     tmpbV = tmpbV.where( tmpbV < 0., -1.e20 )
     mskinfV = tmpbV.where( tmpbV==tmpbV.max('s_rho'), 0.e0 )
     mskinfV = mskinfV.where( mskinfV==0.e0, 1.e0 )
     dzinfV = -tmpbV.where( tmpbV==tmpbV.max('s_rho'), 0.e0 ).sum('s_rho')
     
     # time varying 3d fields:
     for ll in np.arange(mtime):

        Saaa = msksupT*SO_red.isel(time=ll)
        Sbbb = mskinfT*SO_red.isel(time=ll)
        Sccc = ( Saaa.sum('s_rho')*dzinfT + Sbbb.sum('s_rho')*dzsupT ) / (dzinfT+dzsupT)
        SO_z[ll,kk,:,:] = Sccc.values

        Taaa = msksupT*THETAO_red.isel(time=ll)
        Tbbb = mskinfT*THETAO_red.isel(time=ll)
        Tccc = ( Taaa.sum('s_rho')*dzinfT + Tbbb.sum('s_rho')*dzsupT ) / (dzinfT+dzsupT)
        THETAO_z[ll,kk,:,:] = Tccc.values

        Uaaa = msksupU*UX_red.isel(time=ll)
        Ubbb = mskinfU*UX_red.isel(time=ll)
        Uccc = ( Uaaa.sum('s_rho')*dzinfU + Ubbb.sum('s_rho')*dzsupU ) / (dzinfU+dzsupU)
        UX_z[ll,kk,:,:] = Uccc.values

        Vaaa = msksupV*VY_red.isel(time=ll)
        Vbbb = mskinfV*VY_red.isel(time=ll)
        Vccc = ( Vaaa.sum('s_rho')*dzinfV + Vbbb.sum('s_rho')*dzsupV ) / (dzinfV+dzsupV)
        VY_z[ll,kk,:,:] = Vccc.values

   #----------
   # Create new xarray dataset including all useful variables:
   # reshaping (x,y) as 1-dimensional (sxy)

   nxy=(jmax-jmin+1)*(imax-imin+1)

   ds = xr.Dataset(
      {
       "SO":        (["time", "z", "sxy"], np.reshape( SO_z, (mtime,mz,nxy)) ),
       "THETAO":    (["time", "z", "sxy"], np.reshape( THETAO_z, (mtime,mz,nxy)) ),
       "UX":        (["time", "z", "sxy"], np.reshape( UX_z, (mtime,mz,nxy)) ),
       "VY":        (["time", "z", "sxy"], np.reshape( VY_z, (mtime,mz,nxy)) ),
       "TAUX":      (["time", "sxy"], np.reshape( TAUX.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TAUY":      (["time", "sxy"], np.reshape( TAUY.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "ZOS":       (["time", "sxy"], np.reshape( ZOS.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "TOB":       (["time", "sxy"], np.reshape( TOB.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SOB":       (["time", "sxy"], np.reshape( SOB.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "FICESHELF": (["time", "sxy"], np.reshape( FICESHELF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "DYDRFLF":   (["time", "sxy"], np.reshape( DYDRFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "THDRFLF":   (["time", "sxy"], np.reshape( THDRFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HADRFLF":   (["time", "sxy"], np.reshape( HADRFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "MSFTBAROT": (["time", "sxy"], np.reshape( MSFTBAROT.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "HFDS":      (["time", "sxy"], np.reshape( HFDS.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOATRLI":  (["time", "sxy"], np.reshape( WFOATRLI.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "WFOSICOR":  (["time", "sxy"], np.reshape( WFOSICOR.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SICONC":    (["time", "sxy"], np.reshape( SICONC.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVOL":     (["time", "sxy"], np.reshape( SIVOL.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIUX":      (["time", "sxy"], np.reshape( SIUX.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "SIVY":      (["time", "sxy"], np.reshape( SIVY.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, (mtime,nxy)) ),
       "LEVOFT":    (["z", "sxy"], np.reshape( LEVOFT.values, (mz,nxy)) ),
       "LEVOFU":    (["z", "sxy"], np.reshape( LEVOFU.values, (mz,nxy)) ),
       "LEVOFV":    (["z", "sxy"], np.reshape( LEVOFV.values, (mz,nxy)) ),
       "SFTFLF":    (["sxy"], np.reshape( SFTFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "DEPFLF":    (["sxy"], np.reshape( DEPFLF.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "DEPTHO":    (["sxy"], np.reshape( DEPTHO.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKT":   (["sxy"], np.reshape( DOMMSKT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKU":   (["sxy"], np.reshape( DOMMSKU.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, nxy) ),
       "DOMMSKV":   (["sxy"], np.reshape( DOMMSKV.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1)).values, nxy) ),
       "lonT":      (["sxy"], np.reshape( lonT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "lonU":      (["sxy"], np.reshape( lonU.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, nxy) ),
       "lonV":      (["sxy"], np.reshape( lonV.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1)).values, nxy) ),
       "latT":      (["sxy"], np.reshape( latT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "latU":      (["sxy"], np.reshape( latU.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, nxy) ),
       "latV":      (["sxy"], np.reshape( latV.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1)).values, nxy) ),
       "dxT":       (["sxy"], np.reshape( dxT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "dyT":       (["sxy"], np.reshape( dyT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "thetaT":    (["sxy"], np.reshape( thetaT.isel(xi_rho=slice(imin,imax+1),eta_rho=slice(jmin,jmax+1)).values, nxy) ),
       "thetaU":    (["sxy"], np.reshape( thetaU.isel(xi_u=slice(imin,imax+1),eta_u=slice(jmin,jmax+1)).values, nxy) ),
       "thetaV":    (["sxy"], np.reshape( thetaV.isel(xi_v=slice(imin,imax+1),eta_v=slice(jmin,jmax+1)).values, nxy) ),
       "depTUV":    (['z'], newdepth)
      },
      coords={
      "time": ncT.time.values,
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
