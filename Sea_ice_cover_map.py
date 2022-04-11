# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:56:55 2021

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


access_pr_file = 'Sea_ice_cover.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
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

annual_mean_sic = annual_mean(sum_year(sic))

#Selection if ice or not======================================================

def selection(array):
    """
    if sic > 0.2 ==> sic = 1 , means there is ice on this grid cell
    else sic = 0, means there is no ice on this grid cell
    """
    for j in range(len(lat)):
        for k in range(len(lon)):
            if array[j,k] > 0.2:
                array[j,k] = 1
            else:
                array[j,k] = 0
    return array

#Plots =======================================================================
def plot_map(lon,lat, lb, ub, step, units ,args, name, which, arrow, my_color):
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
    
    
    # Open figure
    fig=plt.figure("fig", figsize = (6, 6), dpi = 500)
  
    # Draw coastlines, country boundaries, fill continents, meridians & parallels
    map.drawcoastlines(linewidth = 0.25)
    map.drawmeridians(np.arange(0, 360, 45), linewidth=0.5,latmax=90)#,labels=[True,True,True,True],fontsize=8)
    map.drawparallels(np.arange(-90, 90, 10), linewidth=0.5)
    circle = map.drawmapboundary(linewidth=1, color='k')
    circle.set_clip_on(False)
    
    if which == "double":
        # Create a contourf object called "cs"
        norm = colors.DivergingNorm(vmin=lb, vcenter= 0 , vmax=ub)#!!!
        
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    elif which == "simple":
        # Create a contourf object called "cs"
        norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
        
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    

    plt.title(title,loc='center')
    #plt.text(cs, 2, 0.5, '(b)', fontsize=13, fontweight = 'bold')
    #plt.annotate('(a)', xy=(0.03, 0.92), xycoords='axes fraction',fontsize=13, fontweight = 'bold')
    cax = fig.add_axes([0.91, 0.14, 0.04, 0.7])
    cbar = fig.colorbar(cs, cax=cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(units,fontsize=13, fontweight = 'bold')


    
#plot_map(lon,lat, 0, 1, 0.1 ,'fraction of ice in the cell', annual_mean_sic, "sea ice cover", "simple", "max" , "Blues")
plot_map(lon,lat, 0, 1, 0.1 ,'fraction of ice in the cell', selection(annual_mean_sic), "sea ice cover", "simple", "max" , "Blues")

