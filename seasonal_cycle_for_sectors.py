# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 10:35:59 2021

@author: noepi
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct  1 15:33:24 2021

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


"""
Weddell Sea                         480:768
Indian Ocean                        768:1056
Pacific Ocean                       1056:1344
Ross Sea                            1344:1632%1440
Amundsen Bellingshausen Seas        192:480
"""


access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon_RS=np.concatenate((np.array(dset['longitude'][1344:1440]),np.array(dset['longitude'][0:192])))
lon_ABS=np.array(dset['longitude'][192:480])
lon_WS=np.array(dset['longitude'][480:768])
lon_IO=np.array(dset['longitude'][768:1056])
lon_PO=np.array(dset['longitude'][1056:1344])
lon_all= np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr_RS = np.concatenate((np.array(dset['mtpr'][0:512,0,:,1344:1440]), np.array(dset['mtpr'][0:512,0,:,0:192])), axis = 2)*sec_year
mtpr_ABS = np.array(dset['mtpr'][0:512,0,:,192:480])*sec_year
mtpr_WS = np.array(dset['mtpr'][0:512,0,:,480:768])*sec_year
mtpr_IO = np.array(dset['mtpr'][0:512,0,:,768:1056])*sec_year
mtpr_PO = np.array(dset['mtpr'][0:512,0,:,1056:1344])*sec_year
mtpr_all = np.array(dset['mtpr'][0:512,0,:,:])*sec_year             
dset.close()                                                                  

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd_RS = np.concatenate((np.array(dset['mvimd'][0:512,0,:,1344:1440]), np.array(dset['mvimd'][0:512,0,:,0:192])), axis = 2)*sec_year                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
mvimd_ABS = np.array(dset['mvimd'][0:512,0,:,192:480])*sec_year
mvimd_WS = np.array(dset['mvimd'][0:512,0,:,480:768])*sec_year
mvimd_IO = np.array(dset['mvimd'][0:512,0,:,768:1056])*sec_year
mvimd_PO = np.array(dset['mvimd'][0:512,0,:,1056:1344])*sec_year
mvimd_all = np.array(dset['mvimd'][0:512,0,:,:])*sec_year

dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer_RS = np.concatenate((np.array(dset['mer'][0:512,0,:,1344:1440]), np.array(dset['mer'][0:512,0,:,0:192])), axis = 2)*sec_year 
mer_ABS = np.array(dset['mer'][0:512,0,:,192:480])*sec_year
mer_WS = np.array(dset['mer'][0:512,0,:,480:768])*sec_year
mer_IO = np.array(dset['mer'][0:512,0,:,768:1056])*sec_year
mer_PO = np.array(dset['mer'][0:512,0,:,1056:1344])*sec_year
mer_all = np.array(dset['mer'][0:512,0,:,:])*sec_year           
dset.close()


#negative values indicate evaporation and positive values indicate deposition.
#access_pr_file = 'Mean_snow_evaporation_rate.nc'    
#dset = nc.Dataset(access_pr_file)
#mser = np.concatenate((np.array(dset['mser'][0:512,0,:,1344:1440]), np.array(dset['mser'][0:512,0,:,0:192])), axis = 2)*sec_year                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
#dset.close()


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

grid_cell_area_RS = area_grid(lat, lon_RS)
grid_cell_area_WS = area_grid(lat, lon_WS)
grid_cell_area_ABS = area_grid(lat, lon_ABS)
grid_cell_area_IO = area_grid(lat, lon_IO)
grid_cell_area_PO = area_grid(lat, lon_PO)
grid_cell_area_all = area_grid(lat, lon_all)


total_area_of_earth_RS = grid_cell_area_RS.sum(['latitude','longitude'])
total_area_of_earth_WS = grid_cell_area_WS.sum(['latitude','longitude'])
total_area_of_earth_ABS = grid_cell_area_ABS.sum(['latitude','longitude'])
total_area_of_earth_IO = grid_cell_area_IO.sum(['latitude','longitude'])
total_area_of_earth_PO = grid_cell_area_PO.sum(['latitude','longitude'])
total_area_of_earth_all = grid_cell_area_all.sum(['latitude','longitude'])


def sum_year(array):
    somme_year = np.zeros((42,121,array.shape[2]))
    j = -1
    for i in range(0,504):            #somme 1979-fin 2020
        if (i)%12 == 0:
            j += 1
        somme_year[j,:,:] += (array[i,:,:])/12
    
    return somme_year



mer_RS_by_year = sum_year(mer_RS)
mer_WS_by_year = sum_year(mer_WS)  
mer_ABS_by_year = sum_year(mer_ABS)  
mer_IO_by_year = sum_year(mer_IO)  
mer_PO_by_year = sum_year(mer_PO)  
mer_all_by_year = sum_year(mer_all)  
   
mtpr_RS_by_year = sum_year(mtpr_RS)
mtpr_WS_by_year = sum_year(mtpr_WS)
mtpr_ABS_by_year = sum_year(mtpr_ABS)
mtpr_IO_by_year = sum_year(mtpr_IO)
mtpr_PO_by_year = sum_year(mtpr_PO)
mtpr_all_by_year = sum_year(mtpr_all)

mvimd_RS_by_year = sum_year(mvimd_RS)
mvimd_WS_by_year = sum_year(mvimd_WS)
mvimd_ABS_by_year = sum_year(mvimd_ABS)
mvimd_IO_by_year = sum_year(mvimd_IO)
mvimd_PO_by_year = sum_year(mvimd_PO)
mvimd_all_by_year = sum_year(mvimd_all)




def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(len(array))
    for i in range(0,len(array)):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])     #[kg/year]
    return weighted_mean    
    
    
mer_RS_weighted_mean = weight(mer_RS_by_year, grid_cell_area_RS, total_area_of_earth_RS)
mer_WS_weighted_mean = weight(mer_WS_by_year, grid_cell_area_WS, total_area_of_earth_WS)
mer_ABS_weighted_mean = weight(mer_ABS_by_year, grid_cell_area_ABS, total_area_of_earth_ABS)
mer_IO_weighted_mean = weight(mer_IO_by_year, grid_cell_area_IO, total_area_of_earth_IO)
mer_PO_weighted_mean = weight(mer_PO_by_year, grid_cell_area_PO, total_area_of_earth_PO)
mer_all_weighted_mean = weight(mer_all_by_year, grid_cell_area_all, total_area_of_earth_all)
mer = [-mer_RS_weighted_mean,-mer_WS_weighted_mean,-mer_ABS_weighted_mean,-mer_IO_weighted_mean,-mer_PO_weighted_mean,-mer_all_weighted_mean]

mtpr_RS_weighted_mean = weight(mtpr_RS_by_year, grid_cell_area_RS, total_area_of_earth_RS)
mtpr_WS_weighted_mean = weight(mtpr_WS_by_year, grid_cell_area_WS, total_area_of_earth_WS)
mtpr_ABS_weighted_mean = weight(mtpr_ABS_by_year, grid_cell_area_ABS, total_area_of_earth_ABS)
mtpr_IO_weighted_mean = weight(mtpr_IO_by_year, grid_cell_area_IO, total_area_of_earth_IO)
mtpr_PO_weighted_mean = weight(mtpr_PO_by_year, grid_cell_area_PO, total_area_of_earth_PO)
mtpr_all_weighted_mean = weight(mtpr_all_by_year, grid_cell_area_all, total_area_of_earth_all)
mtpr = [mtpr_RS_weighted_mean, mtpr_WS_weighted_mean, mtpr_ABS_weighted_mean, mtpr_IO_weighted_mean, mtpr_PO_weighted_mean, mtpr_all_weighted_mean ]

mvimd_RS_weighted_mean = weight(mvimd_RS_by_year, grid_cell_area_RS, total_area_of_earth_RS)
mvimd_WS_weighted_mean = weight(mvimd_WS_by_year, grid_cell_area_WS, total_area_of_earth_WS)
mvimd_ABS_weighted_mean = weight(mvimd_ABS_by_year, grid_cell_area_ABS, total_area_of_earth_ABS)
mvimd_IO_weighted_mean = weight(mvimd_IO_by_year, grid_cell_area_IO, total_area_of_earth_IO)
mvimd_PO_weighted_mean = weight(mvimd_PO_by_year, grid_cell_area_PO, total_area_of_earth_PO)
mvimd_all_weighted_mean = weight(mvimd_all_by_year, grid_cell_area_all, total_area_of_earth_all)
mvimd = [-mvimd_RS_weighted_mean, -mvimd_WS_weighted_mean, -mvimd_ABS_weighted_mean, -mvimd_IO_weighted_mean, -mvimd_PO_weighted_mean, -mvimd_all_weighted_mean]

PE_RS_weighted_mean = mtpr_RS_weighted_mean - (-mer_RS_weighted_mean)
PE_WS_weighted_mean = mtpr_WS_weighted_mean - (-mer_WS_weighted_mean)
PE_ABS_weighted_mean = mtpr_ABS_weighted_mean - (-mer_ABS_weighted_mean)
PE_IO_weighted_mean = mtpr_IO_weighted_mean - (-mer_IO_weighted_mean)
PE_PO_weighted_mean = mtpr_PO_weighted_mean - (-mer_PO_weighted_mean)
PE_all_weighted_mean = mtpr_all_weighted_mean - (-mer_all_weighted_mean)
PE = [PE_RS_weighted_mean, PE_WS_weighted_mean, PE_ABS_weighted_mean, PE_IO_weighted_mean, PE_PO_weighted_mean, PE_all_weighted_mean]
#Figures====================================================================================================
#Figure avec multi-year trend-----------------------------------------------------------------------------

def years(nbr, start_year):
    years = np.zeros(nbr)
    for i in range(nbr):
        years[i] = start_year+i
    return years
Years = years(42,1979)    
    
def lin_reg (x, y):
    linreg = linregress(x, y)
    m = linreg.slope
    b = linreg.intercept
    p = linreg.pvalue
    return(m,b,p)


mer_linreg = np.zeros((6,3))
mtpr_linreg = np.zeros((6,3))
mvimd_linreg = np.zeros((6,3))
PE_linreg = np.zeros((6,3))
for i in range(len(mer)):
    mer_linreg[i] = lin_reg(Years, mer[i])
    mtpr_linreg[i] = lin_reg(Years, mtpr[i])
    mvimd_linreg[i] =  lin_reg(Years, mvimd[i])
    PE_linreg[i] = lin_reg(Years, PE[i])


label_sectors = ["RS", "WS", "ABS", "IO", "PO", "all"]
color = ["darkorange", "red", "royalblue", "green", "darkviolet", "cyan" ]


fig = plt.figure(figsize = (16,24), dpi = 500)
ax = plt.subplot(411)
for i in range(len(mer)):
    ax.plot(Years, mer[i] , linestyle = "dashed", label = label_sectors[i], color = color[i])  
    ax.plot(Years, mer_linreg[i][0]*Years + mer_linreg[i][1], color = color[i])#, label= "a = {:.2f}, p = {:.3f}".format(mer_linreg[0], mer_linreg[2]))
ax.set_ylabel("[$mm/year$]", fontsize = 30)   #km3/year
ax.set_xlabel("Years", fontsize = 30)
ax.set_title("Evaporation", fontsize = 30)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=15)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)

 
ax = plt.subplot(412)
for i in range(len(mtpr)):
    ax.plot(Years, mtpr[i] , linestyle = "dashed", label = label_sectors[i], color = color[i])
    ax.plot(Years, mtpr_linreg[i][0]*Years + mtpr_linreg[i][1], color = color[i])#, label= "a = {:.2f}, p = {:.3f}".format(mtpr_linreg[0], mtpr_linreg[2]))
ax.set_ylabel("[$mm/year$]", fontsize = 30)   #km3/year
ax.set_xlabel("Years", fontsize = 30) 
ax.set_title("Precipitation", fontsize = 30)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=15)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)



ax = plt.subplot(413)
for i in range(len(mvimd)):
    ax.plot(Years, mvimd[i] , linestyle = "dashed", label = label_sectors[i], color = color[i])
    ax.plot(Years, mvimd_linreg[i][0]*Years + mvimd_linreg[i][1], color = color[i])#, label= "a = {:.2f}, p = {:.3f}".format(mvimd_linreg[0], mvimd_linreg[2]))
ax.set_ylabel("[$mm/year$]", fontsize = 30)   #km3/year 
ax.set_xlabel("Years", fontsize = 30)
ax.set_title("Moisture Convergence", fontsize = 30)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=15)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)



ax = plt.subplot(414)
for i in range(len(PE)):
    ax.plot(Years, PE[i] , linestyle = "dashed", label = label_sectors[i], color = color[i])
    ax.plot(Years, PE_linreg[i][0]*Years + PE_linreg[i][1], color = color[i])#, label= "a = {:.2f}, p = {:.3f}".format(PE_linreg[0], PE_linreg[2]))
ax.set_ylabel("[$mm/year$]", fontsize = 30)   #km3/year
ax.set_xlabel("Years", fontsize = 30)
ax.set_title("P - E", fontsize = 30)
ax.grid()
ax.tick_params(axis='both', which='major', labelsize=15)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)



fig.tight_layout()
plt.show()    

