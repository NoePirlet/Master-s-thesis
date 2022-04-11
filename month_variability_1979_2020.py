# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 17:05:06 2021

@author: noepi
"""


import numpy as np
import netCDF4 as nc


import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap, addcyclic
import matplotlib.colors as colors
from scipy.stats import linregress


#Import data------------------------------------------------------------------------------------
access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:512,0,:,:])*2592000
dset.close()

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:512,0,:,:])*2592000
dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:512,0,:,:])*2592000
dset.close()

access_pr_file = 'Total_column_water.nc'
dset = nc.Dataset(access_pr_file)
tcw = np.array(dset['tcw'][0:512,0,:,:])
dset.close()

#negative values indicate evaporation and positive values indicate deposition.
access_pr_file = 'Mean_snow_evaporation_rate.nc'    
dset = nc.Dataset(access_pr_file)
mser = np.array(dset['mser'][0:512,0,:,:])*2592000
dset.close()

#------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------
   
        

def repartion_by_months (array,start_month,end_month,unite):
    jan_lst = []
    feb_lst = []
    mar_lst = []
    apr_lst = []
    may_lst = []
    jun_lst = []
    jul_lst = []
    aug_lst = []
    sep_lst = []
    octb_lst = []
    nov_lst = []
    dec_lst = []
  
    for i in range(start_month,end_month):
        if (i+1)%12 == 1:
            jan_lst.append(array[i]*unite)
        if (i+1)%12 == 2:
            feb_lst.append(array[i]*unite)
        if (i+1)%12 == 3:
            mar_lst.append(array[i]*unite)
        if (i+1)%12 == 4:
            apr_lst.append(array[i]*unite)
        if (i+1)%12 == 5:
            may_lst.append(array[i]*unite)
        if (i+1)%12 == 6:
            jun_lst.append(array[i]*unite)
        if (i+1)%12 == 7:
            jul_lst.append(array[i]*unite)
        if (i+1)%12 == 8:
            aug_lst.append(array[i]*unite)
        if (i+1)%12 == 9:
            sep_lst.append(array[i]*unite)
        if (i+1)%12 == 10:
            octb_lst.append(array[i]*unite)
        if (i+1)%12 == 11:
            nov_lst.append(array[i]*unite)
        if (i+1)%12 == 0:
            dec_lst.append(array[i]*unite)
            
        
    return [jan_lst,feb_lst,mar_lst,apr_lst,may_lst,jun_lst,jul_lst,aug_lst,sep_lst,octb_lst,nov_lst,dec_lst] 


km3 = 1/10e12   #conversion en km3

def convert(lst):
    return [ -i for i in lst ]


mer_lst = repartion_by_months(weight(mer, grid_cell_area, total_area_of_earth),0,504,1) 
mtpr_lst = repartion_by_months(weight(mtpr, grid_cell_area, total_area_of_earth),0,504,1)
mvimd_lst = repartion_by_months(weight(mvimd, grid_cell_area, total_area_of_earth),0,504,1)
mser_lst = repartion_by_months(weight(mser, grid_cell_area, total_area_of_earth),0,504,1)


def absolute(lst):
    return [abs(ele) for ele in lst]






#Figures====================================================================================================
years = np.zeros(42)
for i in range(42):
    years[i] = 1979+i


#12 graph 4 curves----------------------------------------------------------------------------------------------------------------
month = ["January", "February" , "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
fig = plt.figure(figsize=(24, 24), dpi = 500)
rows = 6
columns = 2
grid = plt.GridSpec(rows, columns, wspace = .25, hspace = .5)
for i in range(rows*columns):
    exec (f"plt.subplot(grid{[i]})")
    plt.plot(years, absolute(mer_lst[i]), linestyle = "dashed", label = "E")
    plt.plot(years, (mtpr_lst[i]), linestyle = "dashed", label = "P")
    plt.plot(years, absolute(mvimd_lst[i]), linestyle = "dashed", label = "MC")
    plt.plot(years, absolute(mser_lst[i]), linestyle = "dashed", label = "SE")
    plt.xticks(np.arange(1980,2021,5))
    plt.ylabel("[$kg/month$] (abs)")   #km3/year
    plt.xlabel("Years")
    plt.grid()
    plt.title("Terms in budget eq / {0} 1979-2020".format(month[i]))
plt.legend()    
    

#3 graph 12 curves----------------------------------------------------------------------------------------------------------------------


fig = plt.figure(figsize = (18,18), dpi = 100)
ax = plt.subplot(311)
for i in range(len(mer_lst)):
    ax.plot(years, absolute(mer_lst[i]), linestyle = "dashed", label = month[i])
ax.grid()
ax.set_title("Evaporation", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=15)


ax = plt.subplot(312)
for i in range(len(mer_lst)):
    ax.plot(years, absolute(mtpr_lst[i]), linestyle = "dashed", label = month[i])
ax.grid()
ax.set_title("Precipitation", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=15)


ax = plt.subplot(313)
for i in range(len(mer_lst)):
    ax.plot(years, absolute(mvimd_lst[i]), linestyle = "dashed", label = month[i])
ax.grid()
ax.set_title("Moisture convergence", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=15)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

fig.tight_layout()
plt.show()   



#faire graph par mois---------------------------------------------------------------------------------------------------------------------
def saison(array):
    DJF = np.zeros(42)
    MAM = np.zeros(42)
    JJA = np.zeros(42)
    SON = np.zeros(42)
    
    for i in range(len(array[0])):
        DJF[i] = (array[11][i] + array[0][i] + array[1][i])/3
        MAM[i] = (array[2][i] + array[3][i] + array[4][i])/3
        JJA[i] = (array[5][i] + array[6][i] + array[7][i])/3
        SON[i] = (array[8][i] + array[9][i] + array[10][i])/3
    
    return [DJF, MAM, JJA, SON]


def lin_reg (x, y):
    linreg = linregress(x, y)
    m = linreg.slope
    b = linreg.intercept
    p = linreg.pvalue
    return(m,b,p)





saisons = ["DJF \n (Austral summer)", "MAM \n (Austral Autumn)", "JJA \n (Austral Winter)","SON \n (Austral Spring)"]
color = ["darkorange", "red", "royalblue", "green"]


fig = plt.figure(figsize = (18,18), dpi = 100)
ax = plt.subplot(311)
for i in range(4):
    ax.plot(years, lin_reg(years, absolute(saison(mer_lst)[i]))[0]*years +  lin_reg(years, absolute(saison(mer_lst)[i]))[1], color = color[i], label = "a = {:.2f}, p = {:.2f}".format(lin_reg(years, absolute(saison(mer_lst)[i]))[0], lin_reg(years, absolute(saison(mer_lst)[i]))[2]))
    ax.plot(years, absolute(saison(mer_lst)[i]), linestyle = "dashed", color = color[i])
    
ax.grid()
ax.set_title("Evaporation", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
ax.set_xticklabels
ax.tick_params(axis='both', which='major', labelsize=15)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

ax = plt.subplot(312)
for i in range(4):
    ax.plot(years, lin_reg(years, absolute(saison(mtpr_lst)[i]))[0]*years +  lin_reg(years, absolute(saison(mtpr_lst)[i]))[1], color = color[i], label = "a = {:.2f}, p = {:.2f}".format(lin_reg(years, absolute(saison(mtpr_lst)[i]))[0], lin_reg(years, absolute(saison(mtpr_lst)[i]))[2]))
    ax.plot(years, absolute(saison(mtpr_lst)[i]), linestyle = "dashed", color = color[i])
    
ax.grid()
ax.set_title("Precipitation", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
#ax.set_xticklabels(np.arange(1980,2021,5), fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

ax = plt.subplot(313)
for i in range(4):
    ax.plot(years, lin_reg(years, absolute(saison(mvimd_lst)[i]))[0]*years +  lin_reg(years, absolute(saison(mvimd_lst)[i]))[1], color = color[i], label = "a = {:.2f}, p = {:.2f}".format(lin_reg(years, absolute(saison(mvimd_lst)[i]))[0], lin_reg(years, absolute(saison(mvimd_lst)[i]))[2]))
    ax.plot(years, absolute(saison(mvimd_lst)[i]), linestyle = "dashed", label = saisons[i], color = color[i])
    
ax.grid()
ax.set_title("Moisture convergence", fontsize=30)
ax.set_ylabel("[$kg/month$] (abs)", fontsize=30)
ax.set_xlabel("Years", fontsize=30)
#ax.set_xticklabels(np.arange(1980,2021,5), fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

fig.tight_layout()
plt.show()   


 
    













