# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 10:59:13 2021

@author: noepi
"""

import numpy as np
import netCDF4 as nc


import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.colors as colors

#from seasonal_cycle import year_to_month, year_array, sum_area


#Total water column [kg/m^2]
access_pr_file = 'day_01_31_hour_00_month_january_tcw.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
tcw_day_jan = np.array(dset['tcw'][0:84,:,:]) 
dset.close()


access_pr_file = 'day_01_28_hour_00_month_feb_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_feb = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_march_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_mar = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_30_hour_00_month_april_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_apr = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_may_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_may = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_30_hour_00_month_june_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_jun = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_july_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_jul = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_aout_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_aug = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_30_hour_00_month_sept_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_sep = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_octb_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_oct = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_30_hour_00_month_nov_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_nov = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

access_pr_file = 'day_01_31_hour_00_month_dec_tcw.nc'
dset = nc.Dataset(access_pr_file)
tcw_day_dec = np.array(dset['tcw'][0:84,:,:]) 
dset.close()

dt = 2592000    #secondes d'un mois

def difference(tcw_vector):
    tcw_01 = np.zeros((42,121,1440))
    tcw_31 = np.zeros((42,121,1440))
    dtcw = np.zeros((42,121,1440))
    for i in range(0,42):
        tcw_01[i,:,:] = tcw_vector[i*2,:,:]
        tcw_31[i,:,:] = tcw_vector[i*2+1,:,:]
        dtcw[i,:,:] = (tcw_31[i,:,:]-tcw_01[i,:,:])/dt
    return dtcw

dtcw_jan = difference(tcw_day_jan)
dtcw_feb = difference(tcw_day_feb)
dtcw_mar = difference(tcw_day_mar)
dtcw_apr = difference(tcw_day_apr)
dtcw_may = difference(tcw_day_may)
dtcw_jun = difference(tcw_day_jun)
dtcw_jul = difference(tcw_day_jul)
dtcw_aug = difference(tcw_day_aug)
dtcw_sep = difference(tcw_day_sep)
dtcw_oct = difference(tcw_day_oct)
dtcw_nov = difference(tcw_day_nov)
dtcw_dec = difference(tcw_day_dec)


#Weighted grid----------------------------------------------------------------------


def area_grid(lat, lon):
    """
    Calculate the area of each grid cell
    Area is in square meters
    
    Input
    -----------
    lat: vector of latitude in degrees
    lon: vector of longitude in degrees
    
    Output
    -----------
    area: grid-cell area in square-meters with dimensions, [lat,lon]
    
    Notes
    -----------
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """
    from numpy import meshgrid, deg2rad, gradient, cos
    from xarray import DataArray

    xlon, ylat = meshgrid(lon, lat)
    R = 6371000

    dlat = deg2rad(gradient(ylat, axis=0))
    dlon = deg2rad(gradient(xlon, axis=1))

    dy = dlat * R
    dx = dlon * R * cos(deg2rad(ylat))

    area = dy * dx

    xda = DataArray(
        area,
        dims=["latitude", "longitude"],
        coords={"latitude": lat, "longitude": lon},
        attrs={
            "long_name": "area_per_pixel",
            "description": "area per pixel",
            "units": "m^2",
        },
    )
    return xda

grid_cell_area = area_grid(lat, lon)
total_area_of_earth = grid_cell_area.sum(['latitude','longitude'])



def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(len(array))
    for i in range(0,len(array)):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])     #[kg/year]
    return weighted_mean    

dtcw_jan_weighted = weight(dtcw_jan, grid_cell_area, total_area_of_earth)
dtcw_feb_weighted = weight(dtcw_feb, grid_cell_area, total_area_of_earth)
dtcw_mar_weighted = weight(dtcw_mar, grid_cell_area, total_area_of_earth)
dtcw_apr_weighted = weight(dtcw_apr, grid_cell_area, total_area_of_earth)
dtcw_may_weighted = weight(dtcw_may, grid_cell_area, total_area_of_earth)
dtcw_jun_weighted = weight(dtcw_jun, grid_cell_area, total_area_of_earth)
dtcw_jul_weighted = weight(dtcw_jul, grid_cell_area, total_area_of_earth)
dtcw_aug_weighted = weight(dtcw_aug, grid_cell_area, total_area_of_earth)
dtcw_sep_weighted = weight(dtcw_sep, grid_cell_area, total_area_of_earth)
dtcw_oct_weighted = weight(dtcw_oct, grid_cell_area, total_area_of_earth)
dtcw_nov_weighted = weight(dtcw_nov, grid_cell_area, total_area_of_earth)
dtcw_dec_weighted = weight(dtcw_dec, grid_cell_area, total_area_of_earth)

dtcw_all = [dtcw_jan_weighted, dtcw_feb_weighted ,dtcw_mar_weighted ,dtcw_apr_weighted, dtcw_may_weighted,\
            dtcw_jun_weighted ,dtcw_jul_weighted ,dtcw_aug_weighted ,dtcw_sep_weighted, dtcw_oct_weighted,\
            dtcw_nov_weighted, dtcw_dec_weighted ]

dtcw_mean = np.zeros(42)
for i in range(42):
    for j in range(12):
        dtcw_mean[i] += dtcw_all[j][i]
    dtcw_mean[i]/12
    
    
years = np.zeros(42)
for i in range(42):
    years[i] = 1979 + i

months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


#Figure-------------------------------------------------------------------------



fig = plt.figure(dpi = 100)
ax = plt.subplot(111)

for i in range (len(dtcw_all)):
    ax.plot(years, dtcw_all[i] ,linestyle = "dashed",  label = "{}".format(months[i])) 
ax.plot(years, dtcw_mean, label = "mean", color= 'black')
ax.grid()
ax.set_title("dtcw/dt (1979-2020)")
ax.set_ylabel("[kg/month]")
ax.set_xlabel("Years")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show() 

#Map------------------------------------------------------------------------------
# Determine hemisphere(s) to plot
if np.max(lat) < 0:
  hemi = ("s")
elif np.min(lat) > 0:
  hemi = ("n")
else:
  hemi = ("n", "s")
  
# Define a figure name
colmap = "BrBG"
extmethod = "both"#both
varin='PE'                     #PE , mer or mtpr
hemisphere='s' #!!! 
units='$\mathregular{kg/(m^2*month)*10^{(-6)}}$'
lb=-5
ub=6
filein = '1979'               #moi c'est 1979-2021
filein1 = 'diff between tcw - rigth '#!!!1979-2019  
fname = filein.split("/")[-1] 
# -----------------------
#  Create the colorbar |
# -----------------------
# If lower and upper bounds aren't defined, determine them from the field
if lb is None:
  lb = np.max([1e9, np.min(field)])
if ub is None:
  ub = np.max(field)

nsteps=0.5#!!!
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
def plot(args, name) :
      hemisphere = hemi
      fname = filein.split("/")[-1] + "_" + varin + "_" + hemisphere + "_t"  
      #title = filein1.split("/")[-1]
      title = name  
      
      field = np.squeeze(args[:, :])
      # Open figure
      fig=plt.figure("fig", figsize = (6, 6))
    
      # Draw coastlines, country boundaries, fill continents, meridians & parallels
      map.drawcoastlines(linewidth = 0.25)
      map.drawmeridians(np.arange(0, 360, 45), linewidth=0.5,latmax=90)#,labels=[True,True,True,True],fontsize=8)
      map.drawparallels(np.arange(-90, 90, 10), linewidth=0.5)
      circle = map.drawmapboundary(linewidth=1, color='k')
      circle.set_clip_on(False)
      # Create a contourf object called "cs"
      norm = colors.DivergingNorm(vmin=lb, vcenter=0, vmax=ub)#!!!
      
      cs = map.contourf(x, y, field, clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
      cs.cmap.set_under(col_under)
      cs.cmap.set_over(col_over)
    
      plt.title(title,loc='center')
      #plt.text(cs, 2, 0.5, '(b)', fontsize=13, fontweight = 'bold')
      plt.annotate('', xy=(0.03, 0.92), xycoords='axes fraction',fontsize=13, fontweight = 'bold')
      cax = fig.add_axes([0.91, 0.14, 0.04, 0.7])
      cbar = fig.colorbar(cs, cax=cax)
      cbar.ax.tick_params(labelsize=12)
      cbar.set_label(units,fontsize=13, fontweight = 'bold')
      # Save figure
      filename    = fname  
      imageformat = "png"  
      dpi         = 500     
      #plt.savefig(filename + "." + imageformat, bbox_inches = "tight", dpi = dpi)
      #print('Figure ' + filename + "." + imageformat + " printed")
      
     



plot(left_side_budget, 'Change of the the precipitable water (dtcw/dt) for Jan 1979 ')    #You can plot PE,mer, mptr, mvimd, ...


