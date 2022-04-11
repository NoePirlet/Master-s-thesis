# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 20:59:43 2021

@author: noepi
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 17:04:21 2021

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
mtpr_RS = np.concatenate((np.array(dset['mtpr'][0:512,0,:,1344:1440]), np.array(dset['mtpr'][0:512,0,:,0:192])), axis = 2)*sec_month
mtpr_ABS = np.array(dset['mtpr'][0:512,0,:,192:480])*sec_month
mtpr_WS = np.array(dset['mtpr'][0:512,0,:,480:768])*sec_month
mtpr_IO = np.array(dset['mtpr'][0:512,0,:,768:1056])*sec_month
mtpr_PO = np.array(dset['mtpr'][0:512,0,:,1056:1344])*sec_month
mtpr_all = np.array(dset['mtpr'][0:512,0,:,:])*sec_month             
dset.close()                                                                  

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd_RS = np.concatenate((np.array(dset['mvimd'][0:512,0,:,1344:1440]), np.array(dset['mvimd'][0:512,0,:,0:192])), axis = 2)*sec_month                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
mvimd_ABS = np.array(dset['mvimd'][0:512,0,:,192:480])*sec_month
mvimd_WS = np.array(dset['mvimd'][0:512,0,:,480:768])*sec_month
mvimd_IO = np.array(dset['mvimd'][0:512,0,:,768:1056])*sec_month
mvimd_PO = np.array(dset['mvimd'][0:512,0,:,1056:1344])*sec_month
mvimd_all = np.array(dset['mvimd'][0:512,0,:,:])*sec_month

dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer_RS = np.concatenate((np.array(dset['mer'][0:512,0,:,1344:1440]), np.array(dset['mer'][0:512,0,:,0:192])), axis = 2)*sec_month
mer_ABS = np.array(dset['mer'][0:512,0,:,192:480])*sec_month
mer_WS = np.array(dset['mer'][0:512,0,:,480:768])*sec_month
mer_IO = np.array(dset['mer'][0:512,0,:,768:1056])*sec_month
mer_PO = np.array(dset['mer'][0:512,0,:,1056:1344])*sec_month
mer_all = np.array(dset['mer'][0:512,0,:,:])*sec_month           
dset.close()



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



def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(len(array))
    for i in range(0,len(array)):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])     #[kg/year]
    return weighted_mean    
    
    
mer_RS_weighted_mean = weight(mer_RS, grid_cell_area_RS, total_area_of_earth_RS)
mer_WS_weighted_mean = weight(mer_WS, grid_cell_area_WS, total_area_of_earth_WS)
mer_ABS_weighted_mean = weight(mer_ABS, grid_cell_area_ABS, total_area_of_earth_ABS)
mer_IO_weighted_mean = weight(mer_IO, grid_cell_area_IO, total_area_of_earth_IO)
mer_PO_weighted_mean = weight(mer_PO, grid_cell_area_PO, total_area_of_earth_PO)
mer_all_weighted_mean = weight(mer_all, grid_cell_area_all, total_area_of_earth_all)
#mer = [-mer_RS_weighted_mean,-mer_WS_weighted_mean,-mer_ABS_weighted_mean,-mer_IO_weighted_mean,-mer_PO_weighted_mean,-mer_all_weighted_mean]
mer = [-mer_all_weighted_mean,-mer_WS_weighted_mean, -mer_IO_weighted_mean, -mer_PO_weighted_mean, -mer_RS_weighted_mean, -mer_ABS_weighted_mean]

mtpr_RS_weighted_mean = weight(mtpr_RS, grid_cell_area_RS, total_area_of_earth_RS)
mtpr_WS_weighted_mean = weight(mtpr_WS, grid_cell_area_WS, total_area_of_earth_WS)
mtpr_ABS_weighted_mean = weight(mtpr_ABS, grid_cell_area_ABS, total_area_of_earth_ABS)
mtpr_IO_weighted_mean = weight(mtpr_IO, grid_cell_area_IO, total_area_of_earth_IO)
mtpr_PO_weighted_mean = weight(mtpr_PO, grid_cell_area_PO, total_area_of_earth_PO)
mtpr_all_weighted_mean = weight(mtpr_all, grid_cell_area_all, total_area_of_earth_all)
#mtpr = [mtpr_RS_weighted_mean, mtpr_WS_weighted_mean, mtpr_ABS_weighted_mean, mtpr_IO_weighted_mean, mtpr_PO_weighted_mean, mtpr_all_weighted_mean ]
mtpr = [mtpr_all_weighted_mean,mtpr_WS_weighted_mean, mtpr_IO_weighted_mean, mtpr_PO_weighted_mean, mtpr_RS_weighted_mean, mtpr_ABS_weighted_mean]

mvimd_RS_weighted_mean = weight(mvimd_RS, grid_cell_area_RS, total_area_of_earth_RS)
mvimd_WS_weighted_mean = weight(mvimd_WS, grid_cell_area_WS, total_area_of_earth_WS)
mvimd_ABS_weighted_mean = weight(mvimd_ABS, grid_cell_area_ABS, total_area_of_earth_ABS)
mvimd_IO_weighted_mean = weight(mvimd_IO, grid_cell_area_IO, total_area_of_earth_IO)
mvimd_PO_weighted_mean = weight(mvimd_PO, grid_cell_area_PO, total_area_of_earth_PO)
mvimd_all_weighted_mean = weight(mvimd_all, grid_cell_area_all, total_area_of_earth_all)
#mvimd = [-mvimd_RS_weighted_mean, -mvimd_WS_weighted_mean, -mvimd_ABS_weighted_mean, -mvimd_IO_weighted_mean, -mvimd_PO_weighted_mean, -mvimd_all_weighted_mean]
mvimd = [-mvimd_all_weighted_mean,-mvimd_WS_weighted_mean, -mvimd_IO_weighted_mean, -mvimd_PO_weighted_mean, -mvimd_RS_weighted_mean, -mvimd_ABS_weighted_mean]


PE_RS_weighted_mean = mtpr_RS_weighted_mean - (-mer_RS_weighted_mean)
PE_WS_weighted_mean = mtpr_WS_weighted_mean - (-mer_WS_weighted_mean)
PE_ABS_weighted_mean = mtpr_ABS_weighted_mean - (-mer_ABS_weighted_mean)
PE_IO_weighted_mean = mtpr_IO_weighted_mean - (-mer_IO_weighted_mean)
PE_PO_weighted_mean = mtpr_PO_weighted_mean - (-mer_PO_weighted_mean)
PE_all_weighted_mean = mtpr_all_weighted_mean - (-mer_all_weighted_mean)
#PE = [PE_RS_weighted_mean, PE_WS_weighted_mean, PE_ABS_weighted_mean, PE_IO_weighted_mean, PE_PO_weighted_mean, PE_all_weighted_mean]
PE = [PE_all_weighted_mean, PE_WS_weighted_mean, PE_IO_weighted_mean, PE_PO_weighted_mean, PE_RS_weighted_mean, PE_ABS_weighted_mean]

 
#--------------------------------------------------------------------------------------------------------------


def repartion_by_months (array,start_month,end_month,unite):
    jan = 0
    feb = 0
    mar = 0
    apr = 0
    may = 0
    jun = 0
    jul = 0
    aug = 0
    sep = 0
    octb = 0
    nov = 0
    dec = 0
    std = np.zeros(12)
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
            jan_lst.append(array[i])
        if (i+1)%12 == 2:
            feb_lst.append(array[i])
        if (i+1)%12 == 3:
            mar_lst.append(array[i])
        if (i+1)%12 == 4:
            apr_lst.append(array[i])
        if (i+1)%12 == 5:
            may_lst.append(array[i])
        if (i+1)%12 == 6:
            jun_lst.append(array[i])
        if (i+1)%12 == 7:
            jul_lst.append(array[i])
        if (i+1)%12 == 8:
            aug_lst.append(array[i])
        if (i+1)%12 == 9:
            sep_lst.append(array[i])
        if (i+1)%12 == 10:
            octb_lst.append(array[i])
        if (i+1)%12 == 11:
            nov_lst.append(array[i])
        if (i+1)%12 == 0:
            dec_lst.append(array[i])
            
    jan = np.sum(jan_lst) 
    std[0] = np.std(jan_lst)*unite
    feb = np.sum(feb_lst)
    std[1] = np.std(feb_lst)*unite
    mar = np.sum(mar_lst)
    std[2] = np.std(mar_lst)*unite
    apr = np.sum(apr_lst)
    std[3] = np.std(apr_lst)*unite
    may = np.sum(may_lst)
    std[4] = np.std(may_lst)*unite
    jun = np.sum(jun_lst)
    std[5] = np.std(jun_lst)*unite
    jul = np.sum(jul_lst)
    std[6] = np.std(jul_lst)*unite
    aug = np.sum(aug_lst)
    std[7] = np.std(aug_lst)*unite
    sep = np.sum(sep_lst)
    std[8] = np.std(sep_lst)*unite
    octb = np.sum(octb_lst)
    std[9] = np.std(octb_lst)*unite
    nov = np.sum(nov_lst)
    std[10] = np.std(nov_lst)*unite
    dec = np.sum(dec_lst)
    std[11] = np.std(dec_lst)*unite

        
    jan /= ((end_month-start_month)//12)/unite
    feb /= ((end_month-start_month)//12)/unite
    mar /= ((end_month-start_month)//12)/unite
    apr /= ((end_month-start_month)//12)/unite
    may /= ((end_month-start_month)//12)/unite
    jun /= ((end_month-start_month)//12)/unite
    jul /= ((end_month-start_month)//12)/unite
    aug /= ((end_month-start_month)//12)/unite
    sep /= ((end_month-start_month)//12)/unite
    octb /= ((end_month-start_month)//12)/unite
    nov /= ((end_month-start_month)//12)/unite
    dec /= ((end_month-start_month)//12)/unite
        
    return ([jan,feb,mar,apr,may,jun,jul,aug,sep,octb,nov,dec], std) 


mer_sort_month = np.zeros((6, 12))
mtpr_sort_month = np.zeros((6, 12))
mvimd_sort_month = np.zeros((6, 12))
PE_sort_month = np.zeros((6, 12))
mer_std = np.zeros((6,12))
mtpr_std = np.zeros((6,12))
mvimd_std = np.zeros((6,12))
PE_std = np.zeros((6,12))
for i in range(len(mer)):
    mer_sort_month[i], mer_std[i] = repartion_by_months(mer[i],0,504,1) 
    mtpr_sort_month[i], mtpr_std[i] = repartion_by_months(mtpr[i],0,504,1)
    mvimd_sort_month[i], mvimd_std[i] = repartion_by_months(mvimd[i],0,504,1)
    PE_sort_month[i] = (np.abs(mtpr_sort_month[i]) - np.abs(mer_sort_month[i]))
    PE_std[i] = np.sqrt((mtpr_std[i])**2 + (mer_std[i])**2)  #pas sur que ca soit bon ca devrait etre la somme ou autre chose 

#Figures====================================================================================================
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
#label_sectors = ["RS", "WS", "ABS", "IO", "PO", "all"]
label_sectors = ["all", "WS", "IO", "PO", "RS", "ABS"]
#Figure avec courbes saisonnières--------------------------------------------------------------------------

fig = plt.figure(figsize = (16,10) , dpi = 500)
for i in range(len(mer)):
    ax = plt.subplot(231+i)
    ax.errorbar(months, mer_sort_month[i], yerr = mer_std[i],  label = "E", linestyle = "dashed")     #positif
    ax.errorbar(months, mtpr_sort_month[i], yerr = mtpr_std[i],  label = "P", linestyle = "dashed")    #neg
    ax.errorbar(months, PE_sort_month[i], yerr = PE_std[i],  label = "P-E", linestyle = "dashed")
    ax.errorbar(months, mvimd_sort_month[i], yerr = mvimd_std[i],  label = "MC", linestyle = "dashed")          #(positif)
    ax.grid()
    ax.yaxis.set_ticks([0,10,20,30,40,50,60,70,80,90])
    #ax.yaxis.set_major_locator(plt.LinearLocator())
    ax.set_title(label_sectors[i])
    ax.set_ylabel("[$mm/month$]")
    ax.set_xlabel("Month")
    print("mer = {0}, mtpr = {1}, mvimd = {2}, PE = {3}".format(np.mean(mer_sort_month[i]), np.mean(mtpr_sort_month[i]), np.mean(mvimd_sort_month[i]), np.mean(PE_sort_month[i])))
    if i == 2:
        ax.legend()

 
