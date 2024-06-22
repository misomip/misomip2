# 2021-03 : Initial code [N. Jourdain, IGE-CNRS]
#============================================================================================
import numpy as np
from scipy import interpolate

#============================================================================================
def vertical_interp(original_depth,interpolated_depth):
   """ Find upper and lower bound indices for simple vertical 1d interpolation

   """
   if ( original_depth[1] < original_depth[0] ):
     ll_kupward = True
   else:
     ll_kupward = False
   nn = np.size(interpolated_depth)
   kinf=np.zeros(nn,dtype='int')
   ksup=np.zeros(nn,dtype='int')
   for k in np.arange(nn):
      knear = np.argmin( np.abs( original_depth - interpolated_depth[k] ) )
      if (original_depth[knear] > interpolated_depth[k]):
        ksup[k] = knear
        if ll_kupward:
          kinf[k] = np.min([ np.size(original_depth)-1, knear+1 ])
        else:
          kinf[k] = np.max([ 0, knear-1 ])
      else:
        kinf[k] = knear
        if ll_kupward:
          ksup[k] = np.max([ 0, knear-1 ])
        else:
          ksup[k] = np.min([ np.size(original_depth)-1, knear+1 ])
   return (kinf,ksup)

#============================================================================================
def calc_z (h, zice, theta_s, theta_b, hc, N, zeta=None, Vstretching=4):
    """
    Given ROMS grid variables, calculate s-coordinates, stretching curves, and
    z-coordinates. Assumes Vtransform = 2.

    Input (all straight out of grid file and *.in file):
    h, zice = 2D arrays containing values for bathymetry and ice shelf draft.
              Both have dimension latitude x longitude.
    theta_s, theta_b, hc, N = scalar parameters
    zeta = optional 2D array containing values for sea surface height
    Vstretching = optional integer showing stretching transfomration, 2 or 4

    Output:
    z = 3D array containing negative z-coordinate values for depth on the rho 
        grid; dimension depth x latitude x longitude
    s = 1D array of s-coordinate values
    C = 1D array of stretching curve values

    source: https://github.com/kuechenrole/antarctic_melting/blob/master/src/tools/calc_z.py
            by Ole Richter

    Follows the method of scoord_zice.m and stretching.m on katabatic
    (in /ds/projects/iomp/matlab_scripts/ROMS_NetCDF/iomp_IAF/)
    which is also explained on the ROMS wiki:
    https://www.myroms.org/wiki/Vertical_S-coordinate.
    """

    alpha = 1.0
    beta = 1.0
    ds = 1.0/N
    lev = np.arange(1,N+1)-0.5
    s = (lev-N)*ds

    if Vstretching == 2:
        Csur = (-np.cosh(theta_s*s) + 1)/(np.cosh(theta_s) - 1)
        Cbot = np.sinh(theta_b*(s+1))/np.sinh(theta_b) - 1
        weight = (s+1)**alpha*(1 + (alpha/beta)*(1 - (s+1)**beta))
        C = weight*Csur + (1-weight)*Cbot
    elif Vstretching == 4:
        C = (1.-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.)
        C = (np.exp(theta_b*C)-1.)/(1.-np.exp(-theta_b))
        
    h = h - abs(zice)

    num_lon = np.size(h, 1)
    num_lat = np.size(h, 0)
    z = np.zeros((N, num_lat, num_lon))
    for k in range(N):
        z0 = (h*C[k] + hc*s[k])/(h + hc)
        if zeta is None:
            z[k,:,:] = h*z0 - abs(zice)
        else:
            z[k,:,:] = (zeta+h)*z0 + zeta - abs(zice)

    return z

#============================================================================================
def horizontal_interp( lat_in_1d, lon_in_1d, mlat_misomip, mlon_misomip, lat_out_1d, lon_out_1d, \
                       var_in_1d, weight=[], threshold=1.e20, skipna=False, filnocvx=False ):
   """ Interpolates one-dimension data horizontally to a 2d numpy array reshaped to the misomip standard (lon,lat) format.

       Method: triangular linear barycentryc interpolation, using nans (i.e. gives nan if any nan in the triangle)

       Input:

       * lon\_in\_1d, lat\_in\_1d: 1d longitude and latitude of data to interpolate [xarray 1d data array]
 
       * mlat\_misomip, mlon\_misomip: misomip grid size (nb points) alond latitude and longitude dimensions
 
       * lon\_out\_1d, lat\_out\_1d: 1d longitude and latitude of the target misomip grid [numpy 1d data array]

       * var\_in\_1d: 1d input data (same dimension as lon\_in\_1d and lat\_in\_1d) [xarray 1d data array]

       * skipna = False to keep nans in interpolation, i.e. gives nan if any triangle node is nan [default]

                = True to find interpolation triangle nodes with non-nan values

       * filnocvx = True to use nearest-neighbor to fill non-convex areas, i.e. for which no triangulation is possible [default]

                  = False to fill non-convex areas with nans 

       * weight = weights used for interpolation [optional, xarray data array]

       * threshold = threshold below which weight value indicates a masked point [default=1.e20]

       Output:

       * numpy data array of dimension (mlat_misomip, mlon_misomip)

   """
   miss=-999999.99 # local variable, output missing values will be nan
   var1d_nonan = var_in_1d.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
   if ( var1d_nonan.size ==0 ):
     out = np.zeros((mlat_misomip, mlon_misomip))*np.nan
   else:
     if skipna:
       lon_in_1d_nonan = lon_in_1d.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
       lat_in_1d_nonan = lat_in_1d.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
       if ( np.size(weight) == 0 ):
         txxxx = interpolate.griddata( (lon_in_1d_nonan.values,lat_in_1d_nonan.values), var1d_nonan.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
       else:
         wgt1d_nonan = weight.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
         wgt1d_nonan = wgt1d_nonan.where( ~np.isnan(wgt1d_nonan) & ~np.isinf(wgt1d_nonan), 0.e0 ) # if nan in mask but not in input data
         prod=var1d_nonan*wgt1d_nonan
         txxxx = interpolate.griddata( (lon_in_1d_nonan.values,lat_in_1d_nonan.values), prod.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         wgt   = interpolate.griddata( (lon_in_1d_nonan.values,lat_in_1d_nonan.values), wgt1d_nonan.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         txxxx = txxxx / wgt
         txxxx[ wgt < threshold ] = miss
       out = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
       if filnocvx:
         if ( np.size(weight) == 0 ):
           # fill non-convex areas with nearest point (whatever its value):
           tssss = interpolate.griddata( (lon_in_1d_nonan.values,lat_in_1d_nonan.values), var1d_nonan, (lon_out_1d,lat_out_1d), method='nearest' )
         else:
           # fill non-convex areas with nearest point having weight >= threshold :
           tmplon=lon_in_1d_nonan.where( wgt1d_nonan >= threshold, drop=True )
           tmplat=lat_in_1d_nonan.where( wgt1d_nonan >= threshold, drop=True )
           tmpvar=var1d_nonan.where( wgt1d_nonan >= threshold, drop=True )
           if ( tmpvar.size == 0 ):
             tssss = np.zeros((mlat_misomip*mlon_misomip)) * np.nan
           else:
             tssss = interpolate.griddata( (tmplon.values,tmplat.values), tmpvar.values, (lon_out_1d,lat_out_1d), method='nearest' )
         tmp   = np.reshape( tssss, (mlat_misomip, mlon_misomip) )
         out[ (np.isnan(out)) | (np.isinf(out)) ] = tmp[ (np.isnan(out)) | (np.isinf(out)) ] # points out of the convex area
       out[ out == miss ] = np.nan # points with weight below threshold
     else:
       if filnocvx:
         if ( np.size(weight) == 0 ):
           txxxx = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), var_in_1d.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         else:
           prod=var_in_1d*weight
           txxxx = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), prod.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           wgt   = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), weight.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=1.e0 )
           txxxx = txxxx / wgt
           txxxx[ wgt < threshold ] = miss
         out = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )
         lon_in_1d_nonan = lon_in_1d.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
         lat_in_1d_nonan = lat_in_1d.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True ) 
         if ( np.size(weight) == 0 ):
           # fill non-convex areas with nearest point:
           tssss = interpolate.griddata( (lon_in_1d_nonan.values,lat_in_1d_nonan.values), var1d_nonan.values, (lon_out_1d,lat_out_1d), method='nearest' )
         else:
           # fill non-convex areas with nearest point having weight >= threshold :
           wgt1d_nonan = weight.where( (~np.isnan(var_in_1d)) & (~np.isinf(var_in_1d)), drop=True )
           wgt1d_nonan = wgt1d_nonan.where( ~np.isnan(wgt1d_nonan) & ~np.isinf(wgt1d_nonan), 0.e0 ) # if nan in mask but not in input data
           tmplon=lon_in_1d_nonan.where( wgt1d_nonan >= threshold , drop=True )
           tmplat=lat_in_1d_nonan.where( wgt1d_nonan >= threshold , drop=True )
           tmpvar=var1d_nonan.where( wgt1d_nonan >= threshold , drop=True )
           tmpwgt=wgt1d_nonan.where( wgt1d_nonan >= threshold , drop=True )
           tssss = interpolate.griddata( (tmplon.values,tmplat.values), tmpvar.values, (lon_out_1d,lat_out_1d), method='nearest' )
           wgtss = interpolate.griddata( (tmplon.values,tmplat.values), tmpwgt.values, (lon_out_1d,lat_out_1d), method='nearest' )
           tssss[ wgtss < threshold ] = miss
         tmp = np.reshape( tssss, (mlat_misomip, mlon_misomip) )
         out[ (np.isnan(out)) | (np.isinf(out)) ] = tmp[ (np.isnan(out)) | (np.isinf(out)) ]
         out[ out == miss ] = np.nan
       else:
         # Simplest form of horizontal tirangular linear interpolation:
         if ( np.size(weight) == 0 ):
           txxxx = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), var_in_1d.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
         else:
           prod=var_in_1d*weight
           txxxx = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), prod.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           wgt   = interpolate.griddata( (lon_in_1d.values,lat_in_1d.values), weight.values, (lon_out_1d,lat_out_1d), method='linear', fill_value=np.nan )
           txxxx = txxxx / wgt
           txxxx[ wgt < threshold ] = np.nan
         out   = np.reshape( txxxx, (mlat_misomip, mlon_misomip) )  
   return out
