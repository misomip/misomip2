import numpy as np
import zenodo_get as zg
import xarray as xr
from scipy import interpolate
from pyproj import Geod

#====================================================================================================
def grid_bounds_oce(region='Amundsen'):
   """ Gives minimum and maximum longitude and latitude for the common MISOMIP2 ocean grid

       region: 'Amundsen' (default) or 'Weddell'

       exemple: [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region='Amundsen')
   """
   if ( region == 'Amundsen' ):
     longitude_min = -140.0
     longitude_max =  -90.0
     latitude_min  =  -76.0
     latitude_max  =  -69.0
   elif ( region == 'Weddell' ):
     longitude_min = -90.0
     longitude_max =   0.0
     latitude_min  = -85.0
     latitude_max  = -60.0
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude_min,longitude_max,latitude_min,latitude_max]

#====================================================================================================
def generate_3d_grid_oce(region='Amundsen'):
   """Generates (longitude, latitude, depth) of the common MISOMIP2 3d ocean grid

      region: 'Amundsen' (default) or 'Weddell'

      exemple: [lon,lat,depth]=generate_3d_grid_oce(region='Amundsen')
   """

   [lonmin,lonmax,latmin,latmax] = grid_bounds_oce(region=region)

   if ( region == 'Amundsen' ):
     longitude=np.arange(lonmin,lonmax+0.1,0.1)
     latitude=np.arange(latmin,latmax+1./30.,1./30.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   elif ( region == 'Weddell' ):
     longitude=np.arange(lonmin,lonmax+1./3.,1./3.)
     latitude=np.arange(latmin,latmax+1./10.,1./10.)
     depth=np.array([0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1500.])
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")

   return [longitude,latitude,depth]

#====================================================================================================
def cell_area(region='Amundsen'):
   """Generates cell area of MISOMIP2 ocean grid [m2]

      region: 'Amundsen' (default) or 'Weddell'

      exemple: areacello=cell_area(region='Amundsen')
   """

   [lon,lat,depth]=generate_3d_grid_oce(region=region) 
   
   if ( region == 'Amundsen' ):
     lon_u = lon+0.5*0.1
     lon_l = lon-0.5*0.1
     lat_u = lat+0.5*1./30.
     lat_l = lat-0.5*1./30.
   elif ( region == 'Weddell' ):
     lon_u = lon+0.5*1./3.
     lon_l = lon-0.5*1./3.
     lat_u = lat+0.5*1./10.
     lat_l = lat-0.5*1./10.
   else:
     sys.exit("~!@#$%^* error : region is not defined, choose either 'Amundsen' or 'Weddell'")
   
   mlon = np.size(lon)
   mlat = np.size(lat)

   areacello = np.zeros((mlat,mlon))

   # Initialize WGS84 ellipsoid model
   geod = Geod(ellps="WGS84")

   for klon in range(mlon):

      longitudes = [ lon_l[klon], lon_u[klon], lon_u[klon], lon_l[klon] ]

      for klat in range(mlat):
      
         latitudes  = [ lat_l[klat], lat_l[klat], lat_u[klat], lat_u[klat] ] 

         area_m2, perimeter_m = geod.polygon_area_perimeter(longitudes, latitudes)
         areacello[klat,klon] = abs(area_m2)

   return areacello   

