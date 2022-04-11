# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 13:51:59 2021

@author: noepi
"""

import numpy as np
import netCDF4 as nc


import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.colors as colors
import xarray as xr
from scipy.stats import linregress

sec_year = 31536000
sec_month = 2592000


access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:512,0,:,:])*sec_year                       #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
dset.close()                                                                   #2592000 si on est en [kg/((m^2)*month)]

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:512,0,:,:])*sec_year
dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:512,0,:,:])*sec_year
dset.close()

access_pr_file = 'Mean_sea_level_pressure.nc'      #[Pa]
dset = nc.Dataset(access_pr_file)
mslp = np.array(dset['msl'][0:512,0,:,:])
dset.close()

access_pr_file = 'SST.nc'      #[K]
dset = nc.Dataset(access_pr_file)
sst = np.array(dset['sst'][0:512,0,:,:])
dset.close()

access_pr_file = 'Sea_ice_cover.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
sic = np.array(dset['siconc'][0:512,0,:,:])        
dset.close()

access_pr_file = '2m_temperature.nc'      #[K]
dset = nc.Dataset(access_pr_file)
sat = np.array(dset['t2m'][0:512,:,:])
dset.close()

#access_pr_file = 'Total_column_water.nc'
#dset = nc.Dataset(access_pr_file)
#tcw = np.array(dset['tcw'][0:512,0,:,:])
#dset.close()

#negative values indicate evaporation and positive values indicate deposition.
#access_pr_file = 'Mean_snow_evaporation_rate.nc'    
#dset = nc.Dataset(access_pr_file)
#mser = np.array(dset['mser'][0:512,0,:,:])*sec_year      
#dset.close()



def sum_year(array):
    somme_year = np.zeros((42,121,1440))
    for i in range(0,504):            #somme 1979-fin 2020
        somme_year[i//12,:,:] += (array[i,:,:])/12
    return somme_year

mer_by_year = sum_year(mer)     #[kg/(m^2*year)]
mtpr_by_year = sum_year(mtpr)
mvimd_by_year = sum_year(mvimd)
#mser_by_year = sum_year(mser)
PE_by_year = mtpr_by_year - (mer_by_year)
mslp_by_year = sum_year(mslp)
sst_by_year = sum_year(sst)
sic_by_year = sum_year(sic)
sat_by_year = sum_year(sat)

#Linear regression-----------------------------------------------------------------------------

def years(nbr):
    years = np.zeros(nbr)
    for i in range(nbr):
        years[i] = 1979 + i
    return years
    
    
def lin_reg (x, y):
    linreg = linregress(x, y)
    m = linreg.slope
    b = linreg.intercept
    p = linreg.pvalue
    return(m,b,p)


def linreg_grid(years, array):
    linreg_g =  np.full((121, 1440), None)
    p_val =  np.full((121, 1440), None)
    for i in range(121):
        for j in range(1440):
            lin_val = lin_reg(years, array[:,i,j])
            linreg_g[i,j] = lin_val[0]
            p_val[i,j] = lin_val[2]
    return linreg_g,p_val

#def linreg_grid_S(years, array):
#    linreg_g =  np.full((121, 1440), None)
#    for i in range(121):
#        for j in range(1440):
#            lin_val = lin_reg(years, array[:,i,j])
#            if lin_val[2] <= 0.05:
#                linreg_g[i,j] = lin_val[0]
#    return linreg_g
        

mer_linreg, mer_pval = linreg_grid(years(42), -mer_by_year) 
mtpr_linreg, mtpr_pval = linreg_grid(years(42), mtpr_by_year)
mvimd_linreg, mvimd_pval =  linreg_grid(years(42), -mvimd_by_year)
#mser_linreg = linreg_grid(years(42), -mser_by_year)[0]
PE_linreg, PE_pval = linreg_grid(years(42), PE_by_year)
mslp_linreg, mslp_pval = linreg_grid(years(42), mslp_by_year)
sst_linreg, sst_pval = linreg_grid(years(42), sst_by_year)
sic_linreg, sic_pval = linreg_grid(years(42), sic_by_year)
sat_linreg, sat_pval = linreg_grid(years(42), sat_by_year)


#Plots =======================================================================
def plot_map(lon,lat, lb, ub, step, units ,args, name, which, arrow, my_color, array_pval):
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
    first_col1 = cmap(10)[0:3]
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
    
    #zm = np.ma.masked_less(p_val, 0.05)
    
    # ----------------
    # Plot the figure |
    # ----------------
    #(args, name, number) :
    #for jt in range(1):
    #jt = number                                  #la fin y a pas car on a pas les données jusqu'a fin 2021 :) 
    hemisphere = hemi
    #fname = filein.split("/")[-1] + "_" + varin + "_" + hemisphere + "_t"  
    #title = filein1.split("/")[-1] 
    title = "{0}".format(name)
    
    field = np.squeeze(args[:, :])
    
    
    # Open figure
    fig=plt.figure("fig", figsize = (6, 6), dpi = 150)
  
    # Draw coastlines, country boundaries, fill continents, meridians & parallels
    map.drawcoastlines(linewidth = 0.75)
    map.drawmeridians(np.arange(12, 348, 72), linewidth=1,latmax=90, dashes=[1, 0], labels=[1, 1, 1, 1])
    map.drawparallels(np.arange(-90, 90, 10), linewidth=0.5)
    circle = map.drawmapboundary(linewidth=1, color='k')
    circle.set_clip_on(False)

    
    if which == "double":
        # Create a contourf object called "cs"
        center = (lb+ub)/2
        norm = colors.DivergingNorm(vmin=lb, vcenter= 0 , vmax=ub)#!!!
        
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        ds = map.contourf(x, y, array_pval, [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    elif which == "simple":
        # Create a contourf object called "cs"
        norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
        
        cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
        ds = map.contourf(x, y, array_pval, [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
        cs.cmap.set_under(col_under)
        cs.cmap.set_over(col_over)
    
    #à mettre pour titre figure
    #plt.title(title,loc='center', pad=20)
    
    
    
    #plt.text(cs, 2, 0.5, '(b)', fontsize=13, fontweight = 'bold')
    #plt.annotate('(a)', xy=(0.03, 0.92), xycoords='axes fraction',fontsize=13, fontweight = 'bold')
    cax = fig.add_axes([0.98, 0.14, 0.04, 0.7])
    cbar = fig.colorbar(cs, cax=cax)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(units,fontsize=13, fontweight = 'bold')

    

plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year/year}$', mer_linreg , "Trends of evaporation (1979-2020)", "double", "both" , "coolwarm", mer_pval)   


plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year/year}$', mtpr_linreg , "Trends of precipitation (1979-2020)", "double", "both" , "coolwarm", mtpr_pval)   

plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year/year}$', mvimd_linreg , "Trends of Moisture convergence (1979-2020)", "double", "both" , "coolwarm", mvimd_pval)   

# 0 = E - P - MD  <==> 0 = -(P-E) + MC <==> P-E = MC

plot_map(lon,lat, -5, 3, 0.5 ,'$\mathregular{Pa/year}$', mslp_linreg , "Trends of Pressure (1979-2020)", "double", "both" , "coolwarm", mslp_pval)   
#plot_map(lon,lat, -0.01, 0.012, 0.002 ,'$\mathregular{K/year}$', sst_linreg , "Trends of SST (1979-2020)", "double", "both" , "coolwarm", sst_pval)   

plot_map(lon,lat, -0.001, 0.0011, 0.0001 ,'fraction of ice in the cell/year', sic_linreg, "Trends of sea ice cover (1979-2020)", "double", "both" , "coolwarm", sic_pval)

plot_map(lon,lat, -0.021, 0.024, 0.003 ,'$\mathregular{K/year}$', sat_linreg , "Trends of SAT (1979-2020)", "double", "both" , "coolwarm", sat_pval)   















    

