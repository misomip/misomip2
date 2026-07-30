import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import glob
from pathlib import Path

#================================================================

dir_MIPkit = './MIPkit' # contains MIPkit-A and MIPkit-W
dir_models = './DATA' # contains all the model netcdf files

mod    = ['NEMO4.0'      , 'MITgcm'       ]
inst   = ['IGE-CNRS-UGA' , 'UCLA-UMD'     ]
mcolor = ['tab:blue'     , 'tab:orange'   ]
Nmod = len(mod)

#================================================================

for region in ['A', 'W']:

   if region == 'A':
      Nmoor = 8 # number of available moorings
      figA, axA = plt.subplots(nrows=3,ncols=3,figsize=(21.0,21.0))
      axs = axA.ravel()
   else:
      Nmoor = 3 # temporary for Weddell Sea (other ones missing from the MIPkit_W)
      figW, axW = plt.subplots(nrows=3,ncols=3,figsize=(21.0,21.0))
      axs = axW.ravel()
  
   plt.rcParams['xtick.top'] = True
   plt.rcParams['xtick.labeltop'] = True
   plt.rcParams['xtick.bottom'] = False
   plt.rcParams['xtick.labelbottom'] = False
   
   alphabet='abcdefghijklmnopqrstuv'
   
   for kmoor in range(Nmoor):
   
      IDmoor = str(kmoor+1)
      file_obs=glob.glob(dir_MIPkit+'/MIPkit-'+region+'/OceMoor'+IDmoor+'_MIPkit'+region+'_*.nc')[0]
      print(file_obs)
      ds_obs = xr.open_dataset(file_obs)
      Tobs_mean = ds_obs.thetao.mean('time')
 
      valid_times = ds_obs['time'].where(ds_obs['thetao'].min('lev').notnull(), drop=True)
      date_obs_beg = valid_times[0].dt.strftime("%Y-%m-%d").values
      year_obs_beg = valid_times[0].dt.strftime("%Y").values
      date_obs_end = valid_times[-1].dt.strftime("%Y-%m-%d").values
      year_obs_end = valid_times[-1].dt.strftime("%Y").values
      print('Obs from ',date_obs_beg,' to ',date_obs_end)

      axs[kmoor].scatter(Tobs_mean,ds_obs.lev,s=100,c='k',zorder=1.0)
      if kmoor in [0,3,6]:   
         axs[kmoor].set_ylabel('Depth (m)',fontsize=24)
      if kmoor < 3:
         axs[kmoor].set_xlabel('Temperature (°C)',fontsize=24)
      axs[kmoor].tick_params(axis='both',labelsize=20)
      axs[kmoor].xaxis.set_label_position('top')
      axs[kmoor].xaxis.tick_top()
      tit='('+alphabet[kmoor]+') Mooring '+IDmoor+'    '+year_obs_beg+'-'+year_obs_end
      if region == 'A':
         axs[kmoor].set_xlim(-2.0,2.0)
         axs[kmoor].set_ylim(990,0) # swap order to revert y axis   
      else:
         axs[kmoor].set_xlim(-2.2,-1.0)
         axs[kmoor].set_ylim(1150,0) # swap order to revert y axis   
      y_min, y_max = plt.ylim()
      axs[kmoor].set_title(tit,fontsize=22,fontweight='bold',y=y_min*1.1)

      for kmod in range(Nmod):
         
         file_mod=glob.glob(dir_models+'/OceMoor'+IDmoor+'_'+inst[kmod]+'_'+mod[kmod]+'_a_Ocean'+region+'-hind_*.nc')[0]
         print(file_mod)
         ds_mod = xr.open_dataset(file_mod)
         Tmod_mean = ds_mod.thetao.sel(time=slice(date_obs_beg,date_obs_end)).mean('time')
 
         axs[kmoor].plot(Tmod_mean,ds_mod.lev,linewidth=1.6,color=mcolor[kmod],zorder=kmoor*1.0/Nmoor)
   
   # customized legend:
   axs[Nmoor].fill([0,1,1,0,0],[0,0,1,1,0],color='white',edgecolor=None)
   axs[Nmoor].scatter(0.05,0.9,s=100,c='k')
   if region == 'A':
      axs[Nmoor].text(0.5,1.0,'Amundsen Sea Moorings',fontsize=24,va='center',ha='center',fontweight='bold')
   elif region == 'W':
       axs[Nmoor].text(0.5,1.0,'Weddell Sea Moorings',fontsize=24,va='center',ha='center',fontweight='bold')
   axs[Nmoor].text(0.15,0.9,'Observations',fontsize=20,va='center',ha='left')
   for kmod in range(Nmod):
      axs[Nmoor].plot([0,0.1],[0.8-0.1*kmod,0.8-0.1*kmod],linewidth=1.6,color=mcolor[kmod])
      axs[Nmoor].text(0.15,0.8-0.1*kmod,inst[kmod]+'_'+mod[kmod]+'_a',fontsize=20,va='center',ha='left')
   for spine in axs[Nmoor].spines.values():
       spine.set_visible(False)
   axs[Nmoor].tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                  labelbottom=False, labeltop=False, labelleft=False, labelright=False)
   
Path('./figures').mkdir(parents=True, exist_ok=True)
figA.savefig('figures/Moorings_A_mean_profiles_MISOMIP2.jpg')   
figA.savefig('figures/Moorings_A_mean_profiles_MISOMIP2.pdf') 
figW.savefig('figures/Moorings_W_mean_profiles_MISOMIP2.jpg')   
figW.savefig('figures/Moorings_W_mean_profiles_MISOMIP2.pdf') 
