# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 12:41:43 2021

@author: noepi
"""

#Import=======================================================================
import numpy as np
import netCDF4 as nc


import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.colors as colors


#Load netCDF (data)===========================================================

access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:512,0,:,:])*365*86400
dset.close()

#positive for moisture divergence and neg for converging
access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:512,0,:,:])*365*86400
dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:512,0,:,:])*365*86400          #2592000 nécessaire car on a [kg/((m^2)*s)] et qu'on veut représenter sur notre map en [kg/((m^2)*month)]
dset.close()

access_pr_file = 'Total_column_water.nc'            #[kg/m^2]
dset = nc.Dataset(access_pr_file)
tcw = np.array(dset['tcw'][0:512,0,:,:])
dset.close()

access_pr_file = 'Mean_sea_level_pressure.nc'       #[P]
dset = nc.Dataset(access_pr_file)
mslp = np.array(dset['msl'][0:512,0,:,:])/100       #nothing
dset.close()

access_pr_file = 'SST.nc'                           #[K]
dset = nc.Dataset(access_pr_file)
sst = np.array(dset['sst'][0:512,0,:,:])
dset.close()

access_pr_file = 'Sea_ice_cover.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
sic = np.array(dset['siconc'][0:512,0,:,:])        
dset.close()



#Mean 1979-2020 vector========================================================
def sum_year(array):
    j = -1
    year_array = np.zeros((42,121,1440))
    for i in range(0,504):
        if i%12 == 0:
            j += 1
        year_array[j,:,:] += array[i,:,:]
    year_array /= 12                 #peut être que sur un cela donne le taux moyen de "x" pour l'année en question c'est donc une moyenne annuel de (i.e: mer) 
    return year_array


           
def annual_mean(array):
    annual_mean_array = np.zeros((121,1440))
    for i in range(0,42):
        annual_mean_array[:,:] += array[i,:,:]
    annual_mean_array /= 42
    return annual_mean_array
        
#Application over arrays======================================================
annual_mean_mtpr = annual_mean(sum_year(mtpr))
annual_mean_mer = annual_mean(sum_year(mer))
annual_mean_mvimd = annual_mean(sum_year(mvimd))
annual_mean_tcw = annual_mean(sum_year(tcw))
annual_mean_PE = annual_mean_mtpr - (-annual_mean_mer)
annual_mean_pressure = annual_mean(sum_year(mslp))
annual_mean_sst = annual_mean(sum_year(sst))
annual_mean_sic = annual_mean(sum_year(sic))

#Plots =======================================================================
def plot_map(lon,lat, lb, ub, step, units ,args, name, which, arrow, my_color, sic_ = None):
    # Determine hemisphere(s) to plot
    if np.max(lat) < 0:
      hemi = ("s")
    elif np.min(lat) > 0:
      hemi = ("n")
    else:
      hemi = ("n", "s")
      
    # Define a figure name
    colmap = my_color
    extmethod = arrow #both
    #varin='twc'                    
    hemisphere='s' #!!! 
    #units='$\mathregular{kg/m^2}$'
    #lb=-40
    #ub=100
    #filein = '1979-2021'               #moi c'est 1979-2021
    #filein1 = 'Mean evaporation rate'#!!!1979-2019  
    #fname = filein.split("/")[-1] 
    
    # -----------------------
    #  Create the colorbar |
    # -----------------------
    # If lower and upper bounds aren't defined, determine them from the field
    if lb is None:
      lb = np.max([1e9, np.min(field)])
    if ub is None:
      ub = np.max(field)
    
    nsteps=step#!!!
    clevs = np.arange(lb, ub, nsteps)
    
    # Load the colormap
    cmap = eval("plt.cm." + colmap)
    
    # Colors for values outside the colorbar
    # Get first and last colors of the colormap
    first_col = cmap(0)[0:3]
    last_col  = cmap(255)[0:3]
    first_col1 = cmap(10)[0:3]#added by Xia
    # Create "over" and "under" colors.
    # Take respectively the latest and first color and
    # darken them by 50%
    col_under = [i / 2.0 for i in first_col]
    col_over  = [i / 2.0 for i in last_col ]
     
    # --------------------
    # Create a projection |
    # --------------------
    # Determine if we are in the northern hemisphere or in the Southern hemisphere.
    if hemisphere == "n":
      boundlat = 45.
      l0 = 0.
    elif hemisphere == "s":
      boundlat = -60.
      l0 = 180.
    else:
      sys.exit("(map_polar_field) Hemisphere unkown")
    
    # Projection name
    projname = hemisphere + "plaea" #'npstere' #'ortho' #
    map = Basemap(projection = projname, boundinglat = boundlat, lon_0 = l0, resolution = 'l',round=True)
    lon, lat = np.meshgrid(lon, lat)
    x, y = map(lon, lat)
    yy = np.arange(0, y.shape[0], 8)
    xx = np.arange(0, x.shape[1], 20)
    points = np.meshgrid(yy, xx)
    # ----------------
    # Plot the figure |
    # ----------------
    #(args, name, number) :
    #for jt in range(1):
    #jt = number                                  #la fin y a pas car on a pas les données jusqu'a fin 2021 :) 
    hemisphere = hemi
    #fname = filein.split("/")[-1] + "_" + varin + "_" + hemisphere + "_t"  
    #title = filein1.split("/")[-1] 
    title = "Annual mean 1979-2020 {0}".format(name)
    
    field = np.squeeze(args[:, :])
    sic = np.squeeze(sic_[:,:])
    
    # Open figure
    fig=plt.figure("fig", figsize = (6, 6), dpi = 500)
  
    # Draw coastlines, country boundaries, fill continents, meridians & parallels
    map.drawcoastlines(linewidth = 0.25)
    map.drawmeridians(np.arange(12, 348, 72), linewidth=1,latmax=90, dashes=[1, 0], labels=[1, 1, 1, 1])
    map.drawparallels(np.arange(-90, 90, 10), linewidth=0.5)
    circle = map.drawmapboundary(linewidth=1, color='k')
    circle.set_clip_on(False)
    
    if which == "double":
        # Create a contourf object called "cs"
        norm = colors.DivergingNorm(vmin=lb, vcenter= 0, vmax=ub)#!!!
        
        si = map.contour(x,y, sic, [0.2,1], colors = "black" , linestyles = "dashed")
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    elif which == "simple":
        # Create a contourf object called "cs"
        norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
        
        si = map.contour(x,y, sic, [0.2,1], colors = "black" , linestyles = "dashed")
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    

    plt.title(title,loc='center', pad=20)
    #plt.text(cs, 2, 0.5, '(b)', fontsize=13, fontweight = 'bold')
    #plt.annotate('(a)', xy=(0.03, 0.92), xycoords='axes fraction',fontsize=13, fontweight = 'bold')
    cax = fig.add_axes([0.98, 0.14, 0.04, 0.7])
    cbar = fig.colorbar(cs, cax=cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(units,fontsize=13, fontweight = 'bold')
    
     

    

#plot_map(lon,lat, 0, 10 , 1 ,'$\mathregular{mm}$', annual_mean_tcw, "tcw", "simple", "max", "Blues")   #[kg/m^2] = [mm]
#plot_map(lon,lat, -100, 740, 70 ,'$\mathregular{mm/year}$', -annual_mean_mvimd, "mvimc", "double", "both", "BrBG")
#plot_map(lon,lat, -100, 740, 70 ,'$\mathregular{mm/year}$', annual_mean_PE, "P-E", "double", "both" , "BrBG")    
#plot_map(lon,lat, 0, 550, 50 ,'$\mathregular{mm/year}$', -annual_mean_mer, "Evaporation", "simple", "max" , "Reds", sic_ = annual_mean_sic) 
#plot_map(lon,lat, 0, 1100, 100 ,'$\mathregular{mm/year}$', annual_mean_mtpr, "Precipitation", "simple", "max" , "Blues")   
plot_map(lon,lat, 271 , 275, 0.5 ,'$\mathregular{K}$', annual_mean_sst, "Sea surface temperature", "simple", "both" , "coolwarm_r", sic_ = annual_mean_sic)   
#plot_map(lon,lat, 982 , 1016, 2 ,'$\mathregular{hPa}$', annual_mean_pressure, "Mean sea level pressure", "simple", "both" , "coolwarm")   

    
    
    
    
    
    
    
    
    
    
    