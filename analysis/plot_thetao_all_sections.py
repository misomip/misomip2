import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.path import Path as mplpath
import glob
from pathlib import Path
import cmocean
import cmocean.cm as cmo
from def_mod_list import define_models

# import model list, institute name and colors:
dir_MIPkit,dir_models,mod,inst,mcolor = define_models()
Nmod = len(mod)

# min max of the colormap:
Tmin = -2.0
Tmax =  2.0

alphabet='abcdefghijklmnopqrstuv'

#for region in ['A', 'W']:
for region in ['A', 'W']:

   if region == 'A':
      Nsec = 2 # number of available sections
   else:
      Nsec = 4 # 
  
   for ksec in range(Nsec):

      IDsec = str(ksec+1) 
   
      file_tmp=glob.glob(dir_MIPkit+'/MIPkit-'+region+'/OceSec'+IDsec+'_MIPkit'+region+'_*.nc')
      Nfiles=len(file_tmp)

      for kf in range(Nfiles): # one file per observed year => one plot per observed year

         nc = 2 # or nc = int(np.sqrt(Nmod))
         nr = int(Nmod/nc)
         fig, axs = plt.subplots(nrows=nr,ncols=nc,figsize=(21.0,10.5*nr))
         axs = axs.ravel()

         file_obs = file_tmp[kf]
         print(file_obs)
         ds_obs = xr.open_dataset(file_obs)
         Tobs = ds_obs.thetao
 
         valid_times = ds_obs['time'].mean()
         date_obs = valid_times.dt.strftime("%Y-%m-%d").values
         year_obs = valid_times.dt.strftime("%Y").values
         print('Obs at ',date_obs,',  year ', year_obs)

         for kmod in range(Nmod):

            file_mod=glob.glob(dir_models+'/OceSec'+IDsec+'_'+inst[kmod]+'_'+mod[kmod]+'_a_Ocean'+region+'-hind_*.nc')
            if not file_mod:
               print('WARNING, NO FILE FOR OceSec'+IDsec+'_'+inst[kmod]+'_'+mod[kmod]+'_a_Ocean'+region+'-hind  >>>> SKIPPED')
            else:
               file_mod=file_mod[0]
               print(file_mod)
               ds_mod = xr.open_dataset(file_mod)
               Tmod = ds_mod.thetao.sel(time=date_obs,method='nearest')
               
               pax = axs[kmod].pcolormesh(ds_mod.x,ds_mod.lev,Tmod,vmin=Tmin,vmax=Tmax,cmap=cmo.thermal,zorder=0.1)
   
               # plot observations : customized rectangular marker scatter plot:
               rect_vertices = [(-0.5, -0.2), (0.5, -0.2), (0.5, 0.2), (-0.5, 0.2)] # to define rectangular markers
               rect_vertices2 = [(-0.5, -0.2), (-0.5, 0.2), (0.5, -0.2), (0.5, 0.2)] # to define rectangular markers
               lateral_codes = [ mplpath.MOVETO, mplpath.LINETO,   # Draw left line
                                 mplpath.MOVETO, mplpath.LINETO ]  # Draw right line
               lateral_edges_marker = mplpath(rect_vertices2, lateral_codes)
               for klon in range(ds_obs.lon.size):
                  # find matching x model coordinate:
                  dis = ( ds_obs.lon[klon] - ds_mod.lon )**2 + ( ds_obs.lat[klon] - ds_mod.lat )**2
                  indx = dis.argmin().values
                  xmatch = ds_mod.x.isel(x=indx).values
                  # Step A: Plot the filled square body with NO edges (linewidths=0)
                  Tplot = Tobs.isel(lon=klon)
                  xplot = np.zeros((ds_obs.lev.where( ((Tplot>-5.0)&(Tplot<30.0)), drop=True).size))+xmatch
                  yplot = ds_obs.lev.where( ((Tplot>-5.0)&(Tplot<30.0)), drop=True)
                  Tplot = Tplot.where( ((Tplot>-5.0)&(Tplot<30.0)), drop=True)
                  axs[kmod].scatter(xplot,yplot,c=Tplot,\
                                    marker=rect_vertices,s=100,cmap=cmo.thermal,vmin=Tmin,vmax=Tmax,\
                                    linewidths=0,zorder=1.0+klon*0.01)
                  # Step B: Overlay the lateral edges using the custom marker   
                  axs[kmod].scatter(xplot,yplot,marker=lateral_edges_marker,s=100,\
                                    edgecolors='gray',linewidths=1,facecolors='none',zorder=1.0+klon*0.01+0.005)    
   
               # general plot setting:
               xtick = np.linspace(0,ds_mod.x.max().values,5).astype('int')   
               xl1 = np.abs(ds_mod.lon.isel(x=xtick).values)
               xl2 = np.abs(ds_mod.lat.isel(x=xtick).values)
               xlab = np.array([])
               for kk in range(len(xl1)):
                  string = f"{xl1[kk]:.1f}"+'°W\n'+f"{xl2[kk]:.1f}"+'°S'
                  xlab = np.append(xlab,string)
               axs[kmod].set_xticks(xtick,xlab)
               axs[kmod].set_ylabel('Depth (m)',fontsize=24)
               axs[kmod].tick_params(axis='both',labelsize=20)
               tit = '('+alphabet[kmod]+') '+inst[kmod]+'_'+mod[kmod]+'_a'
               axs[kmod].set_title(tit,fontsize=22,fontweight='bold')
               ymin,ymax = axs[kmod].get_ylim()
               axs[kmod].set_ylim(ymax,0) # swap order to revert y axis   
   
               # colorbar:
               if ( ( region == 'W' ) & ( ksec == 0) ):
                  cax = axs[kmod].inset_axes([0.25, 0.20, 0.5, 0.05])
               elif ( ( region == 'W' ) & ( ksec == 1) ):
                  cax = axs[kmod].inset_axes([0.40, 0.09, 0.5, 0.05])
               elif ( ( region == 'W' ) & ( ksec == 2) ):
                  cax = axs[kmod].inset_axes([0.05, 0.09, 0.5, 0.05])
               elif ( ( region == 'W' ) & ( ksec == 3) ):
                  cax = axs[kmod].inset_axes([0.35, 0.09, 0.5, 0.05])
               else:
                  cax = axs[kmod].inset_axes([0.25, 0.09, 0.5, 0.05])
               cbar=fig.colorbar(pax, cax=cax, orientation="horizontal")
               if ( ( region == 'W' ) & ( ksec == 0) ):
                  cbar.ax.tick_params(labelsize=16,color='white',labelcolor='white') 
                  cbar.set_label('Potential Temperature (°C)',fontsize=16,color='white')
               else:
                  cbar.ax.tick_params(labelsize=16,color='black',labelcolor='black') 
                  cbar.set_label('Potential Temperature (°C)',fontsize=16,color='black')

         Path('./figures').mkdir(parents=True, exist_ok=True)
         fig.savefig('figures/Section'+IDsec+'_Ocean'+region+'-hind_'+year_obs+'_MISOMIP2.jpg')   
         fig.savefig('figures/Section'+IDsec+'_Ocean'+region+'-hind_'+year_obs+'_MISOMIP2.pdf')
