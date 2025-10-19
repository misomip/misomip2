# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]

#====================================================================================================
def add_standard_attributes_oce(ds,miss=9.969209968386869e36,verbose=False):
  """Define standard netcdf attributes for ocean variables already present in the ds xarray dataset.

     ds : xarray ocean dataset
     miss : missing value (default=9.969209968386869e36)
     verbose : True or False
  
  """

  ## coordinates :

  if ( "longitude" in ds.coords ):
    if verbose: print('define attributes for coordinate longitude')
    ds.longitude.encoding['_FillValue'] = None
    ds.longitude.attrs['units'] = 'degrees_east'
    ds.longitude.attrs['long_name'] = 'Longitude'
    ds.longitude.attrs['standard_name'] = 'longitude'

  if ( "latitude" in ds.coords ):
    if verbose: print('define attributes for coordinate latitude')
    ds.latitude.encoding['_FillValue'] = None
    ds.latitude.attrs['units'] = 'degrees_north'
    ds.latitude.attrs['long_name'] = 'Latitude'
    ds.latitude.attrs['standard_name'] = 'latitude'

  if ( "lev" in ds.coords ):
    if verbose: print('define attributes for coordinate lev')
    ds.depth.encoding['_FillValue'] = None
    ds.depth.attrs['units'] = 'm'
    ds.depth.attrs['positive'] = 'down'
    ds.depth.attrs['long_name'] = 'depth'
    ds.depth.attrs['standard_name'] = 'depth'
    ds.depth.attrs['comment'] = 'common MISOMIP2 grid; increases from the sea surface to the sea floor'

  if ( "time" in ds.coords ):
    if verbose: print('define attributes for coordinate time')
    ds.time.encoding['units'] = 'days since 1900-01-01'
    ds.time.encoding['_FillValue'] = None
    ds.time.attrs['standard_name'] = 'time'
    ds.time.attrs['long_name'] = 'time'


  ## variables :

  if ( ( "longitude" in ds.data_vars ) & ~( "longitude" in ds.coords ) ):
    if verbose: print('define attributes for variable longitude')
    ds.longitude.attrs['units'] = 'degrees_east'
    ds.longitude.attrs['long_name'] = 'Longitude'
    ds.longitude.attrs['standard_name'] = 'longitude'
    ds.longitude.attrs['comment'] = 'common MISOMIP2 grid'

  if ( ( "latitude" in ds.data_vars ) & ~( "latitude" in ds.coords ) ):
    if verbose: print('define attributes for variable latitude')
    ds.latitude.attrs['units'] = 'degrees_north'
    ds.latitude.attrs['long_name'] = 'Latitude'
    ds.latitude.attrs['standard_name'] = 'latitude'
    ds.latitude.attrs['comment'] = 'common MISOMIP2 grid'

  if ( "sftflf" in ds.data_vars ):
     if verbose: print('define attributes for variable sftflf')
     ds.sftflf.attrs['_FillValue'] = miss
     ds.sftflf.attrs['units'] = '%'
     ds.sftflf.attrs['long_name'] = 'Floating Ice Shelf Area Percentage'
     ds.sftflf.attrs['standard_name'] = 'floating\_ice\_shelf\_area\_fraction'
     ds.sftflf.attrs['cell_method'] = 'area: mean'

  if ( "sftof" in ds.data_vars ):
     if verbose: print('define attributes for variable sftof')
     ds.sftof.attrs['_FillValue'] = miss
     ds.sftof.attrs['units'] = '%'
     ds.sftof.attrs['long_name'] = 'Sea Area Percentage at the Surface'
     ds.sftof.attrs['standard_name'] = 'sea\_area\_fraction'
     ds.sftof.attrs['cell_method'] = 'area: mean'
     ds.sftof.attrs['comment'] = 'Equals zero if sftflf is 100% (ice shelf) or if the cell is fully occupied by land, but 100% for a cell fully or partially covered with sea ice'

  if ( "levof" in ds.data_vars ):
     if verbose: print('define attributes for variable levof')
     ds.levof.attrs['_FillValue'] = miss
     ds.levof.attrs['units'] = '%'
     ds.levof.attrs['long_name'] = 'Sea area fraction at each vertical level'
     ds.levof.attrs['standard_name'] = 'TBD'
     ds.levof.attrs['cell_method'] = 'area: mean'

  if ( "deptho" in ds.data_vars ):
     if verbose: print('define attributes for variable deptho')
     ds.deptho.attrs['_FillValue'] = miss
     ds.deptho.attrs['units'] = 'm'
     ds.deptho.attrs['long_name'] = 'Sea Floor Depth Below Geoid' 
     ds.deptho.attrs['standard_name'] = 'sea_floor_depth_below_geoid'
     ds.deptho.attrs['cell_method'] = 'area: mean where 3d ocean'

  if ( "depflf" in ds.data_vars ):
     if verbose: print('define attributes for variable depflf')
     ds.depflf.attrs['_FillValue'] = miss
     ds.depflf.attrs['units'] = 'm'
     ds.depflf.attrs['long_name'] = 'Depth of Floating Ice Base Below Geoid' 
     ds.depflf.attrs['standard_name'] = 'TBD'
     ds.depflf.attrs['cell_method'] = 'area: mean where ice shelf'

  if ( "thetao" in ds.data_vars ):
     if verbose: print('define attributes for variable thetao')
     ds.thetao.attrs['_FillValue'] = miss
     ds.thetao.attrs['units'] = 'degC'
     ds.thetao.attrs['long_name'] = 'Sea Water Potential Temperature'
     ds.thetao.attrs['standard_name'] = 'sea_water_potential_temperature'
     ds.thetao.attrs['cell_method'] = 'volume: mean where ocean; time: monthly mean'
     ds.thetao.attrs['comment'] = 'This is the quantity that approximates the practical salinity traditionally obtained through conductivity measurements (see appendix D of Griffies et al. (2016)'

  if ( "so" in ds.data_vars ):
     if verbose: print('define attributes for variable so')
     ds.so.attrs['_FillValue'] = miss
     ds.so.attrs['units'] = '0.001'
     ds.so.attrs['long_name'] = 'Sea Water Salinity (practical salinity)'
     ds.so.attrs['standard_name'] = 'sea_water_salinity' 
     ds.so.attrs['cell_method'] = 'volume: mean where ocean; time: monthly mean'

  if ( "tob" in ds.data_vars ):
     if verbose: print('define attributes for variable tob')
     ds.tob.attrs['_FillValue'] = miss
     ds.tob.attrs['units'] = 'degC'
     ds.tob.attrs['long_name'] = 'Sea Water Potential Temperature at Sea Floor'
     ds.tob.attrs['standard_name'] = 'sea_water_potential_temperature_at_sea_floor'
     ds.tob.attrs['cell_method'] = 'area: mean where bottom ocean; time: monthly mean'

  if ( "sob" in ds.data_vars ):
     if verbose: print('define attributes for variable sob')
     ds.sob.attrs['_FillValue'] = miss
     ds.sob.attrs['units'] = '0.001'
     ds.sob.attrs['long_name'] = 'Sea Water Salinity at Sea Floor (Practical Salinity)'
     ds.sob.attrs['standard_name'] = 'sea_water_salinity_at_sea_floor'
     ds.sob.attrs['cell_method'] = 'area: mean where bottom ocean; time: monthly mean'
     ds.sob.attrs['comment'] = 'This is practical salinity'

  if ( "uo" in ds.data_vars ):
     if verbose: print('define attributes for variable uo')
     ds.uo.attrs['_FillValue'] = miss
     ds.uo.attrs['units'] = 'm s-1'
     ds.uo.attrs['long_name'] = 'Sea Water X Velocity (Zonal)'
     ds.uo.attrs['standard_name'] = 'sea_water_x_velocity'
     ds.uo.attrs['cell_method'] = 'volume: mean where ocean; time: monthly mean'
     ds.uo.attrs['comment'] = 'This is zonal velocity on the common grid, positive eastward'

  if ( "vo" in ds.data_vars ):
     if verbose: print('define attributes for variable vo')
     ds.vo.attrs['_FillValue'] = miss
     ds.vo.attrs['units'] = 'm s-1'
     ds.vo.attrs['long_name'] = 'Sea Water Y Velocity (Meridional)'
     ds.vo.attrs['standard_name'] = 'sea_water_y_velocity'
     ds.vo.attrs['cell_method'] = 'volume: mean where ocean; time: monthly mean'
     ds.vo.attrs['comment'] = 'This is meridional velocity on the common grid, positive northward'

  if ( "tauuo" in ds.data_vars ):
     if verbose: print('define attributes for variable tauuo')
     ds.tauuo.attrs['_FillValue'] = miss
     ds.tauuo.attrs['units'] = 'N m-2'
     ds.tauuo.attrs['long_name'] = 'Sea Water Downward X Stress'
     ds.tauuo.attrs['standard_name'] = 'downward_x_stress_at_sea_water_surface'
     ds.tauuo.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.tauuo.attrs['comment'] = 'This is the zonal stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf (expressed as a 2D variable) and possibly icebergs and any momentum flux correction.'

  if ( "tauvo" in ds.data_vars ):
     if verbose: print('define attributes for variable tauvo')
     ds.tauvo.attrs['_FillValue'] = miss
     ds.tauvo.attrs['units'] = 'N m-2'
     ds.tauvo.attrs['long_name'] = 'Sea Water Downward Y Stress'
     ds.tauvo.attrs['standard_name'] = 'downward_y_stress_at_sea_water_surface'
     ds.tauvo.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.tauvo.attrs['comment'] = 'This is the meridional stress on the liquid ocean from overlying atmosphere, sea ice, ice shelf (expressed as a 2D variable) and possibly icebergs and any momentum flux correction.'

  if ( "msftbarot" in ds.data_vars ):
     if verbose: print('define attributes for variable vsftbarot')
     ds.msftbarot.attrs['_FillValue'] = miss
     ds.msftbarot.attrs['units'] = 'kg s-1'
     ds.msftbarot.attrs['long_name'] = 'Ocean Barotropic Mass Streamfunction'
     ds.msftbarot.attrs['standard_name'] = 'ocean_barotropic_mass_streamfunction'
     ds.msftbarot.attrs['cell_method'] = 'area: mean; time: monthly mean'
     ds.msftbarot.attrs['comment'] = 'Quasi-barotropic streamfunction as discussed in appendix H26 of Griffies et al. (2016); for Boussinesq models, this is simply the volume barotropic streamfunction times the reference seawater volumic mass; the streamfunction $\Psi$ is computed so that $\partial_y \Psi = U^\rho$ and $\partial_x \Psi = -V^\rho$, where $U^\rho$ and $V^\rho$ are the zonal and meridional barotropic mass transports.'

  if ( "zos" in ds.data_vars ):
     if verbose: print('define attributes for variable zos')
     ds.zos.attrs['_FillValue'] = miss
     ds.zos.attrs['units'] = 'm'
     ds.zos.attrs['long_name'] = 'Sea Surface Height Above Geoid'
     ds.zos.attrs['standard_name'] = 'sea_surface_height_above_geoid'
     ds.zos.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.zos.attrs['comment'] = 'This is the dynamic sea surface height above geoid, i.e. not including steric sea-level changes (see appendix H7 of Griffies et al. (2016)'

  if ( "wfoat" in ds.data_vars ):
     if verbose: print('define attributes for variable wfoat')
     ds.wfoat.attrs['_FillValue'] = miss
     ds.wfoat.attrs['units'] = 'kg m-2 s-1'
     ds.wfoat.attrs['long_name'] = 'Water Flux into Sea Water from Atmosphere'
     ds.wfoat.attrs['standard_name'] = 'TBD'
     ds.wfoat.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.wfoat.attrs['positive'] = 'downward'
     ds.wfoat.attrs['comment'] = 'This is calculated as condensation minus evaporation plus solid and liquid precipitation, only considering the part of these fluxes that enters the sea-ice free portion of the cell, but expressed per area of sea and sea-ice; considering appendix K2-K3 of Griffies et al. (2016), wfoat=pr+prsn+evs; models using virtual salt fluxes are invited to calculate an equivalent freshwater mass flux'

  if ( "flandice" in ds.data_vars ):
     if verbose: print('define attributes for variable flandice')
     ds.flandice.attrs['_FillValue'] = miss
     ds.flandice.attrs['units'] = 'kg m-2 s-1'
     ds.flandice.attrs['long_name'] = 'Water Mass Flux Into Sea Water From Land Ice'
     ds.flandice.attrs['standard_name'] = 'water_flux_into_sea_water_from_land_ice'
     ds.flandice.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.flandice.attrs['positive'] = 'downward'
     ds.flandice.attrs['comment'] = 'This is calculated as runoff from rivers or surface ice-sheet melting, plus iceberg melt, plus ice-shelf melt minus refreezing; considering appendix K2-K3 of Griffies et al. (2016), flandice=friver+ficeberg+ficeshelf; models using virtual salt fluxes are invited to calculate an equivalent freshwater mass flux'

  if ( "fsitherm" in ds.data_vars ):
     if verbose: print('define attributes for variable fsitherm')
     ds.fsitherm.attrs['_FillValue'] = miss
     ds.fsitherm.attrs['units'] = 'kg m-2 s-1'
     ds.fsitherm.attrs['long_name'] = 'Water Flux into Sea Water Due to Sea Ice Thermodynamics'
     ds.fsitherm.attrs['standard_name'] = 'water_flux_into_sea_water_due_to_sea_ice_thermodynamics'
     ds.fsitherm.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.fsitherm.attrs['positive'] = 'downward'
     ds.fsitherm.attrs['comment'] = 'This is the net flux, calculated as sea-ice melt minus sea-ice formation/freezing; this is the flux into the total sea cell (open + sea-ice covered); models using virtual salt fluxes are invited to calculate an equivalent freshwater mass flux'

  if ( "wfocorr" in ds.data_vars ):
     if verbose: print('define attributes for variable wfocorr')
     ds.wfocorr.attrs['_FillValue'] = miss
     ds.wfocorr.attrs['units'] = 'kg m-2 s-1'
     ds.wfocorr.attrs['long_name'] = 'Water Mass Flux Into Sea Water From Salinity Correction'
     ds.wfocorr.attrs['standard_name'] = 'TBD'
     ds.wfocorr.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.wfocorr.attrs['positive'] = 'downward'
     ds.wfocorr.attrs['comment'] = 'This is the flux corresponding to the sea surface salinity restoring/adjustment that is common in global ocean models; it should be set to zero for models with no correction; models using virtual salt fluxes are invited to calculate an equivalent freshwater mass flux; this variable is not officially part of CMIP6 but was used in OMIP'

  if ( "hfds" in ds.data_vars ):
     if verbose: print('define attributes for variable hfds')
     ds.hfds.attrs['_FillValue'] = miss
     ds.hfds.attrs['units'] = 'W m-2'
     ds.hfds.attrs['long_name'] = 'Downward Heat Flux at Sea Water Surface'
     ds.hfds.attrs['standard_name'] = 'TBD'
     ds.hfds.attrs['cell_method'] = 'area: mean where 3d ocean; time: monthly mean'
     ds.hfds.attrs['comment'] = 'This is calculated from the net shortwave and longwave radiative fluxes penetrating into the liquid water, the sensible and latent heat fluxes at the atmosphere--ocean, sea-ice--ocean, ice-shelf--ocean (expressed as a 2D variable) and iceberg--ocean interfaces, including those related to the heat content of runoff or precipitation, and any heat flux correction at the ocean surface; see list of individual fluxes in appendix K4 of Griffies et al. (2016); this variable is similar to the hfds variable in CMIP/OMIP, except that it includes potential heat flux correction'

  if ( "libmassbffl" in ds.data_vars ):
     if verbose: print('define attributes for variable libmassbffl')
     ds.libmassbffl.attrs['_FillValue'] = miss
     ds.libmassbffl.attrs['units'] = 'kg m-2 s-1'
     ds.libmassbffl.attrs['long_name'] = 'Basal Specific Mass Balance of Floating Ice Shelf'
     ds.libmassbffl.attrs['standard_name'] = 'land\_ice\_basal\_specific\_mass\_balance\_flux'
     ds.libmassbffl.attrs['cell_method'] = 'area: mean where ice shelf; time: monthly mean'
     ds.libmassbffl.attrs['comment'] = 'This differs from the ficeshelf term in Griffies et al. (2016), which was the net water mass flux into sea water from ice shelf, i.e. per unit of ocean area, while libmassbffl is per unit of ice-shelf area; positive for melting, negative for refreezing'

  if ( "dydrflf" in ds.data_vars ):
     if verbose: print('define attributes for variable dydrflf')
     ds.dydrflf.attrs['_FillValue'] = miss
     ds.dydrflf.attrs['units'] = 'm s-1'
     ds.dydrflf.attrs['long_name'] = 'Dynamical Driving at the Base of Floating Ice Shelf'
     ds.dydrflf.attrs['standard_name'] = 'TBD'
     ds.dydrflf.attrs['cell_method'] = 'area: mean where ice shelf; time: monthly mean'
     ds.dydrflf.attrs['comment'] = 'This is also referred to as the heat exchange velocity, i.e. friction velocity times heat exchange coefficient'

  if ( "thdrflf" in ds.data_vars ):
     if verbose: print('define attributes for variable thdrflf')
     ds.thdrflf.attrs['_FillValue'] = miss
     ds.thdrflf.attrs['units'] = 'degC'
     ds.thdrflf.attrs['long_name'] = 'Thermal Driving at the Base of Floating Ice Shelf'
     ds.thdrflf.attrs['standard_name'] = 'TBD'
     ds.thdrflf.attrs['cell_method'] = 'area: mean where ice shelf; time: monthly mean'
     ds.thdrflf.attrs['comment'] = 'This is calculated as the potential temperature in the top ocean boundary layer beneath the ice shelf, minus the freezing potential temperature at the ice–ocean interface'

  if ( "hadrflf" in ds.data_vars ):
     if verbose: print('define attributes for variable hadrflf')
     ds.hadrflf.attrs['_FillValue'] = miss
     ds.hadrflf.attrs['units'] = '0.001'
     ds.hadrflf.attrs['long_name'] = 'Haline Driving at the Base of Floating Ice Shelf'
     ds.hadrflf.attrs['standard_name'] = 'TBD'
     ds.hadrflf.attrs['cell_method'] = 'area: mean where ice shelf; time: monthly mean'
     ds.hadrflf.attrs['comment'] = 'This is calculated as the practical salinity in the top ocean boundary layer beneath the ice shelf minus the salinity at the ice–ocean interface'

  if ( "siconc" in ds.data_vars ):
     if verbose: print('define attributes for variable siconc')
     ds.siconc.attrs['_FillValue'] = miss
     ds.siconc.attrs['units'] = '%'
     ds.siconc.attrs['long_name'] = 'Sea-Ice Area Percentage'
     ds.siconc.attrs['standard_name'] = 'sea_ice_area_fraction'
     ds.siconc.attrs['cell_method'] = 'area: mean where sea; time: monthly mean'

  if ( "sivol" in ds.data_vars ):
     if verbose: print('define attributes for variable sivol')
     ds.sivol.attrs['_FillValue'] = miss
     ds.sivol.attrs['units'] = 'm'
     ds.sivol.attrs['long_name'] = 'Sea-Ice Volume per Area'
     ds.sivol.attrs['standard_name'] = 'sea_ice_thickness'
     ds.sivol.attrs['cell_method'] = 'area: mean where sea; time: monthly mean'

  if ( "siu" in ds.data_vars ):
     if verbose: print('define attributes for variable siu')
     ds.siu.attrs['_FillValue'] = miss
     ds.siu.attrs['units'] = 'm s-1'
     ds.siu.attrs['long_name'] = 'X-Component of Sea-Ice Velocity (Zonal)'
     ds.siu.attrs['standard_name'] = 'sea_ice_x_velocity'
     ds.siu.attrs['cell_method'] = 'area: mean where sea-ice; time: monthly mean'
     ds.siu.attrs['comment'] = 'Zonal velocity on the MISOMIP2 grid'

  if ( "siv" in ds.data_vars ):
     if verbose: print('define attributes for variable siv')
     ds.siv.attrs['_FillValue'] = miss
     ds.siv.attrs['units'] = 'm s-1'
     ds.siv.attrs['long_name'] = 'Y-Component of Sea-Ice Velocity (Meridional)'
     ds.siv.attrs['standard_name'] = 'sea_ice_y_velocity'
     ds.siv.attrs['cell_method'] = 'area: mean where sea-ice; time: monthly mean'
     ds.siv.attrs['comment'] = 'Meridional velocity on the MISOMIP2 grid'


#====================================================================================================
def rename_dimensions(ds):
   """ rename dimensions of xarray dataset ds to standard (x,y,z,time) if needed

       usage: 
              dds = rename_dimensions(dds)
   """

   if ( "time_counter" in ds.dims ):
     ds=ds.rename({'time_counter':'time'})
   if ( "t" in ds.dims ):
     ds=ds.rename({'t':'time'})

   if ( "depth" in ds.dims ):
     ds=ds.rename({'depth':'z'})
   if ( "level" in ds.dims ):
     ds=ds.rename({'level':'z'})
   if ( "nav_lev" in ds.dims ):
     ds=ds.rename({'nav_lev':'z'})
   if ( "deptht" in ds.dims ):
     ds=ds.rename({'deptht':'z'})
   if ( "depthu" in ds.dims ):
     ds=ds.rename({'depthu':'z'})
   if ( "depthv" in ds.dims ):
     ds=ds.rename({'depthv':'z'})
   if ( "Z" in ds.dims ):
     ds=ds.rename({'Z':'z'})

   if ( "x_grid_T" in ds.dims ):
     ds=ds.rename({'x_grid_T':'x'})
   if ( "x_grid_U" in ds.dims ):
     ds=ds.rename({'x_grid_U':'x'})
   if ( "x_grid_V" in ds.dims ):
     ds=ds.rename({'x_grid_V':'x'})
   if ( "XC" in ds.dims ):
     ds=ds.rename({'XC':'x'})
   if ( "XG" in ds.dims ):
     ds=ds.rename({'XG':'x'})

   if ( "y_grid_T" in ds.dims ):
     ds=ds.rename({'y_grid_T':'y'})
   if ( "y_grid_U" in ds.dims ):
     ds=ds.rename({'y_grid_U':'y'})
   if ( "y_grid_V" in ds.dims ):
     ds=ds.rename({'y_grid_V':'y'})
   if ( "YC" in ds.dims ):
     ds=ds.rename({'YC':'y'}) 
   if ( "YG" in ds.dims ):
     ds=ds.rename({'YG':'y'})

   return ds
