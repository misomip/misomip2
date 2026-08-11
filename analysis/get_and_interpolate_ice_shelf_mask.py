import numpy as np
import zenodo_get as zg
import xarray as xr
from scipy import interpolate
from grid_utils import generate_3d_grid_oce

#====================================================================================================

# Download latest version of ice shelf mask data if not already there:
doi = '10.5281/zenodo.15863351' # latest version
zg.download(record_or_doi=doi,output_dir='./DATA',file_glob="Mask_Iceshelf_IMBIE2.nc")


#====================================================================================================
# interpolate to MISMOMIP2 regional grids:

ds=xr.open_dataset('./DATA/Mask_Iceshelf_IMBIE2.nc')
mx=ds.x.size 
my=ds.y.size
msk1d = np.reshape( ds.Iceshelf_extrap.values, mx*my )
lon1d = np.reshape( ds.lon.values, mx*my )
lat1d = np.reshape( ds.lat.values, mx*my )

for reg in ['Amundsen', 'Weddell']:

   # define MISOMIP2 grid
   [lon_miso,lat_miso,dep_miso] = generate_3d_grid_oce(region=reg)
   mlon = np.size(lon_miso)
   mlat = np.size(lat_miso)
   lon2d_miso, lat2d_miso = np.meshgrid( lon_miso, lat_miso )
   lon_miso1d = np.reshape( lon2d_miso, mlon*mlat )
   lat_miso1d = np.reshape( lat2d_miso, mlon*mlat )

   # interpolate to nearest neighbour
   tmpx = interpolate.griddata( (lon1d,lat1d), msk1d, (lon_miso1d,lat_miso1d), method='nearest' )
   msk_miso = np.reshape( tmpx, (mlat, mlon) )

   # save to netcdf
   dsmiso = xr.Dataset(
       {
        "mask_ice_shelf":          (["lat", "lon"], np.float32(msk_miso)),
        "NAME":                    (["ID"], ds.NAME.values),
        "MeltRignot":              (["ID"], np.float32(ds.MeltRignot.values)),
        "UncertaintiesRignot":     (["ID"], np.float32(ds.UncertaintiesRignot.values)),
        "AreaRignot":              (["ID"], np.float32(ds.AreaRignot.values)),
        "MeltAdusumilli":          (["ID"], np.float32(ds.MeltAdusumilli.values)),
        "UncertaintiesAdusumilli": (["ID"], np.float32(ds.UncertaintiesAdusumilli.values)),
        "AreaAdusumilli":          (["ID"], np.float32(ds.AreaAdusumilli.values)),
        "MeltPaolo":               (["ID"], np.float32(ds.MeltPaolo.values)),
        "UncertaintiesPaolo":      (["ID"], np.float32(ds.UncertaintiesPaolo.values)),
        "AreaPaolo":               (["ID"], np.float32(ds.AreaPaolo.values)),
        "MeltDavison":             (["ID"], np.float32(ds.MeltDavison.values)),
        "UncertaintiesDavison":    (["ID"], np.float32(ds.UncertaintiesDavison.values)),
       },  
       coords={
          "lon": np.float32(lon_miso),
          "lat": np.float32(lat_miso),
          "ID": np.int64(ds.ID.values),
           },
   )

   dsmiso.lon.encoding['_FillValue'] = None
   dsmiso.lon.attrs['units'] = 'degrees_east'
   dsmiso.lon.attrs['long_name'] = 'Longitude'
   dsmiso.lon.attrs['standard_name'] = 'longitude'

   dsmiso.lat.encoding['_FillValue'] = None
   dsmiso.lat.attrs['units'] = 'degrees_north'
   dsmiso.lat.attrs['long_name'] = 'Latitude'
   dsmiso.lat.attrs['standard_name'] = 'latitude'

   dsmiso.ID.encoding['_FillValue'] = None
   dsmiso.ID.attrs['long_name'] = 'Ice Shelf ID'

   dsmiso["mask_ice_shelf"].attrs['long_name'] = 'Mask of individual ice shelves and their drainage basin'
   dsmiso["NAME"].attrs['long_name'] = 'Ice shelf names'
   dsmiso["MeltRignot"].attrs['long_name'] = 'Ice Shelf Basal Mass Balance from Rignot et al. (2013)'
   dsmiso["UncertaintiesRignot"].attrs['long_name'] = 'Uncertainty of Ice Shelf Basal Mass Balance from Rignot et al. (2013)'
   dsmiso["AreaRignot"].attrs['long_name'] = 'Ice Shelf Area used in Rignot et al. (2013)'
   dsmiso["MeltAdusumilli"].attrs['long_name'] = 'Ice Shelf Basal Mass Balance from Adusumilli et al. (2020)'
   dsmiso["UncertaintiesAdusumilli"].attrs['long_name'] = 'Uncertainty of Ice Shelf Basal Mass Balance from Adusumilli et al. (2020)'
   dsmiso["AreaAdusumilli"].attrs['long_name'] = 'Ice Shelf Area used in Adusumilli et al. (2020)'
   dsmiso["MeltPaolo"].attrs['long_name'] = 'Ice Shelf Basal Mass Balance from Paolo et al. (2023)'
   dsmiso["UncertaintiesPaolo"].attrs['long_name'] = 'Uncertainty of Ice Shelf Basal Mass Balance from Paolo et al. (2023)'
   dsmiso["AreaPaolo"].attrs['long_name'] = 'Ice Shelf Area used in Paolo et al. (2023)'
   dsmiso["MeltDavison"].attrs['long_name'] = 'Ice Shelf Basal Mass Balance from Davison et al. (2023)'
   dsmiso["UncertaintiesDavison"].attrs['long_name'] = 'Uncertainty of Ice Shelf Basal Mass Balance from Davison et al. (2023)'

   for var in ["MeltRignot","UncertaintiesRignot","MeltAdusumilli","UncertaintiesAdusumilli","MeltPaolo","UncertaintiesPaolo","MeltDavison","UncertaintiesDavison"]:
      dsmiso[var].attrs['units'] = 'Gt/yr'
   for var in ["AreaRignot","AreaAdusumilli","AreaPaolo"]:
      dsmiso[var].attrs['units'] = 'km2'

   dsmiso.attrs['project'] = 'MISOMIP2'
   dsmiso.attrs['history'] = 'Generated using get_and_interpolate_ice_shelf_mask.py'
   dsmiso.attrs['references'] = 'Mouginot et al. (2017) : https://doi.org/10.5067/AXE4121732AD ; IMBIE2: https://doi.org/10.5281/zenodo.15863352'

   dsmiso.to_netcdf('DATA/Mask_Ice_Shelves_'+reg+'.nc')
