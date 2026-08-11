import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
from pathlib import Path
from def_mod_list import define_models
from grid_utils import cell_area 
            
# import model list, institute name and colors:
dir_MIPkit,dir_models,mod,inst,mcolor = define_models()
Nmod = len(mod)

alphabet='abcdefghijklmnopqrstuv'

for region in ['A', 'W']:

   if region == 'A':
      filemsk = glob.glob(dir_models+'/Mask_Ice_Shelves_Amundsen.nc')
      areacello = cell_area(region='Amundsen')
      # list of ice shelves to plot:
      ice_shelves_plot = [ "Getz", "Dotson/Philbin_Inlet",  "Crosson", "Thwaites", "Pine_Island", "Cosgrove" ]
   elif region == 'W':
      filemsk = glob.glob(dir_models+'/Mask_Ice_Shelves_Weddell.nc')
      areacello = cell_area(region='Weddell')
      # list of ice shelves to plot:
      ice_shelves_plot = [ "LarsenC", "LarsenD", "Ronne", "Filchner", "Brunt_Stancomb", "Riiser-Larsen" ]
   if not filemsk:
      print('ERROR: no ice shelf mask file >>>> you need to run get_and_interpolate_ice_shelf_mask.py before this script')
      exit()
   else:
      filemsk=filemsk[0]

   nc = 2
   nr = int(len(ice_shelves_plot)/nc)
   fig, axs = plt.subplots(nrows=nr,ncols=nc,figsize=(21.0,6.5*nr))
   axs = axs.ravel()

   dsmsk = xr.open_dataset(filemsk)

   xcell_area = xr.DataArray(data=areacello,dims=["lat", "lon"])

   for kmod in range(Nmod):

      file_mod=glob.glob(dir_models+'/Oce3d_'+inst[kmod]+'_'+mod[kmod]+'_a_Ocean'+region+'-hind_*.nc')
      if not file_mod:
         print('WARNING, NO FILE FOR Oce3d_'+inst[kmod]+'_'+mod[kmod]+'_a_Ocean'+region+'-hind  >>>> SKIPPED')
      else:
         file_mod=file_mod[0]
         print(file_mod)
         ds_mod = xr.open_dataset(file_mod)

         if ( "libmassbffl" in ds_mod.data_vars ):

            for kisf in range(len(ice_shelves_plot)):
   
               # get ice shelf ID and corresponding mask:
               IDisf = dsmsk.ID.where((dsmsk.NAME == ice_shelves_plot[kisf]),drop=True)[0].values
               msk = dsmsk.mask_ice_shelf.where((dsmsk.mask_ice_shelf == IDisf),other=0)/IDisf
               melt = ds_mod.libmassbffl * ds_mod.sftflf * 0.01 * xcell_area * msk
               ice_shelf_area_km2 = int((ds_mod.sftflf * 0.01 * xcell_area * msk).sum(["lon","lat"]).values * 1.e-6)
               bmb = melt.sum(["lon","lat"])*1.e-12*365.25*86400 # Gt/yr
   
               lab = inst[kmod]+'_'+mod[kmod]+'_a ('+str(ice_shelf_area_km2)+r' km^2)' 
               print(lab)
               axs[kisf].plot(ds_mod.time,bmb,color=mcolor[kmod],label=lab)

         else:

            print('WARNING: no "libmassbffl" variable in model '+inst[kmod]+'_'+mod[kmod]+'_a  >>>>>>> SKIPPING IT !!!!!')

   Path('./figures').mkdir(parents=True, exist_ok=True)
   fig.savefig('figures/Ice_shelf_melt_Ocean'+region+'-hind_MISOMIP2.jpg')
   fig.savefig('figures/Ice_shelf_melt_Ocean'+region+'-hind_MISOMIP2.pdf')
