# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 11:06:22 2022

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


access_pr_file = 'pr_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
lat_histo=np.array(dset['lat'][0:25])
lon_all_histo=np.array(dset['lon'][:])
lon_all_histo=np.append(lon_all_histo, 360)
lon_WS_histo=np.concatenate((np.array(lon_all_histo[160:193]),np.array(lon_all_histo[0:7])))
lon_IO_histo=np.array(lon_all_histo[7:46])
lon_PO_histo=np.array(lon_all_histo[46:84])
lon_RS_histo=np.array(lon_all_histo[84:122])
lon_ABS_histo=np.array(lon_all_histo[122:160])
time=np.array(dset['time'][:])
mtpr_all_histo = np.array(dset['pr'][1548:1980,0:25,:])*sec_month 
mtpr_all_histo_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mtpr_all_histo_v2[i,j] = np.append(mtpr_all_histo[i,j],(mtpr_all_histo[i,j,0]+mtpr_all_histo[i,j,-1])/2)
mtpr_WS_histo = np.concatenate((np.array(mtpr_all_histo_v2[:,:,160:193]), np.array(mtpr_all_histo_v2[:,:,0:7])), axis = 2)
mtpr_IO_histo = np.array(mtpr_all_histo_v2[:,:,7:46])
mtpr_PO_histo = np.array(mtpr_all_histo_v2[:,:,46:84])
mtpr_RS_histo = np.array(mtpr_all_histo_v2[:,:,84:122])
mtpr_ABS_histo = np.array(mtpr_all_histo_v2[:,:,122:160])
dset.close()                                                                  

access_pr_file = 'evspsbl_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
mer_all_histo = np.array(dset['evspsbl'][1548:1980,0:25,:])*sec_month 
mer_all_histo_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mer_all_histo_v2[i,j] = np.append(mer_all_histo[i,j],(mer_all_histo[i,j,0]+mer_all_histo[i,j,-1])/2)
mer_WS_histo = np.concatenate((np.array(mer_all_histo_v2[:,:,160:193]), np.array(mer_all_histo_v2[:,:,0:7])), axis = 2)
mer_IO_histo = np.array(mer_all_histo_v2[:,:,7:46])
mer_PO_histo = np.array(mer_all_histo_v2[:,:,46:84])
mer_RS_histo = np.array(mer_all_histo_v2[:,:,84:122])
mer_ABS_histo = np.array(mer_all_histo_v2[:,:,122:160])
dset.close()


access_pr_file = 'pr_Amon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc'
dset = nc.Dataset(access_pr_file)
mtpr_all_v0 = np.array(dset['pr'][600:1032,0:25,:])*sec_month 
mtpr_all = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mtpr_all[i,j] = np.append(mtpr_all_v0[i,j],(mtpr_all_v0[i,j,0]+mtpr_all_v0[i,j,-1])/2)
mtpr_WS = np.concatenate((np.array(mtpr_all[:,:,160:193]), np.array(mtpr_all[:,:,0:7])), axis = 2)
mtpr_IO = np.array(mtpr_all[:,:,7:46])
mtpr_PO = np.array(mtpr_all[:,:,46:84])
mtpr_RS = np.array(mtpr_all[:,:,84:122])
mtpr_ABS = np.array(mtpr_all[:,:,122:160])
dset.close()                                                 

access_pr_file = 'evspsbl_Amon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc'
dset = nc.Dataset(access_pr_file)
mer_all_v0 = np.array(dset['evspsbl'][600:1032,0:25,:])*sec_month 
mer_all = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mer_all[i,j] = np.append(mer_all_v0[i,j],(mer_all_v0[i,j,0]+mer_all_v0[i,j,-1])/2)
mer_WS = np.concatenate((np.array(mer_all[:,:,160:193]), np.array(mer_all[:,:,0:7])), axis = 2)
mer_IO = np.array(mer_all[:,:,7:46])
mer_PO = np.array(mer_all[:,:,46:84])
mer_RS = np.array(mer_all[:,:,84:122])
mer_ABS = np.array(mer_all[:,:,122:160])
dset.close()


mvimd_all_histo = np.load("mvimd_ACCESS_histo.npy")
for i in range(432):
    for j in range(193):
        mvimd_all_histo[i,0,j] = mvimd_all_histo[i,1,j]
mvimd_WS_histo = np.concatenate((np.array(mvimd_all_histo[:,:,160:193]), np.array(mvimd_all_histo[:,:,0:7])), axis = 2)
mvimd_IO_histo = np.array(mvimd_all_histo[:,:,7:46])
mvimd_PO_histo = np.array(mvimd_all_histo[:,:,46:84])
mvimd_RS_histo = np.array(mvimd_all_histo[:,:,84:122])
mvimd_ABS_histo = np.array(mvimd_all_histo[:,:,122:160])

mvimd_all_ssp = np.load("mvimd_ACCESS_ssp.npy")
for i in range(432):
    for j in range(193):
        mvimd_all_ssp[i,0,j] = mvimd_all_ssp[i,1,j]
mvimd_WS_ssp = np.concatenate((np.array(mvimd_all_ssp[:,:,160:193]), np.array(mvimd_all_ssp[:,:,0:7])), axis = 2)
mvimd_IO_ssp = np.array(mvimd_all_ssp[:,:,7:46])
mvimd_PO_ssp = np.array(mvimd_all_ssp[:,:,46:84])
mvimd_RS_ssp = np.array(mvimd_all_ssp[:,:,84:122])
mvimd_ABS_ssp = np.array(mvimd_all_ssp[:,:,122:160])


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

grid_cell_area_RS_histo = area_grid(lat_histo, lon_RS_histo)
grid_cell_area_WS_histo = area_grid(lat_histo, lon_WS_histo)
grid_cell_area_ABS_histo = area_grid(lat_histo, lon_ABS_histo)
grid_cell_area_IO_histo = area_grid(lat_histo, lon_IO_histo)
grid_cell_area_PO_histo = area_grid(lat_histo, lon_PO_histo)
grid_cell_area_all_histo = area_grid(lat_histo, lon_all_histo)


total_area_of_earth_RS_histo = grid_cell_area_RS_histo.sum(['latitude','longitude'])
total_area_of_earth_WS_histo = grid_cell_area_WS_histo.sum(['latitude','longitude'])
total_area_of_earth_ABS_histo = grid_cell_area_ABS_histo.sum(['latitude','longitude'])
total_area_of_earth_IO_histo = grid_cell_area_IO_histo.sum(['latitude','longitude'])
total_area_of_earth_PO_histo = grid_cell_area_PO_histo.sum(['latitude','longitude'])
total_area_of_earth_all_histo = grid_cell_area_all_histo.sum(['latitude','longitude'])




def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(len(array))
    for i in range(0,len(array)):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])     #[kg/year]
    return weighted_mean    
    
    
mer_RS_weighted_mean_histo = weight(mer_RS_histo, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mer_WS_weighted_mean_histo = weight(mer_WS_histo, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mer_ABS_weighted_mean_histo = weight(mer_ABS_histo, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mer_IO_weighted_mean_histo = weight(mer_IO_histo, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mer_PO_weighted_mean_histo = weight(mer_PO_histo, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mer_all_weighted_mean_histo = weight(mer_all_histo_v2, grid_cell_area_all_histo, total_area_of_earth_all_histo)
mer_histo = [mer_all_weighted_mean_histo, mer_WS_weighted_mean_histo, mer_IO_weighted_mean_histo, mer_PO_weighted_mean_histo, mer_RS_weighted_mean_histo, mer_ABS_weighted_mean_histo]

mtpr_RS_weighted_mean_histo = weight(mtpr_RS_histo, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mtpr_WS_weighted_mean_histo = weight(mtpr_WS_histo, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mtpr_ABS_weighted_mean_histo = weight(mtpr_ABS_histo, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mtpr_IO_weighted_mean_histo = weight(mtpr_IO_histo, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mtpr_PO_weighted_mean_histo = weight(mtpr_PO_histo, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mtpr_all_weighted_mean_histo = weight(mtpr_all_histo_v2, grid_cell_area_all_histo, total_area_of_earth_all_histo)
mtpr_histo = [mtpr_all_weighted_mean_histo, mtpr_WS_weighted_mean_histo, mtpr_IO_weighted_mean_histo, mtpr_PO_weighted_mean_histo, mtpr_RS_weighted_mean_histo, mtpr_ABS_weighted_mean_histo]


mvimd_RS_weighted_mean_histo = weight(mvimd_RS_histo, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mvimd_WS_weighted_mean_histo = weight(mvimd_WS_histo, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mvimd_ABS_weighted_mean_histo = weight(mvimd_ABS_histo, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mvimd_IO_weighted_mean_histo = weight(mvimd_IO_histo, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mvimd_PO_weighted_mean_histo = weight(mvimd_PO_histo, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mvimd_all_weighted_mean_histo = weight(mvimd_all_histo, grid_cell_area_all_histo, total_area_of_earth_all_histo)
mvimd_histo = [mvimd_all_weighted_mean_histo, mvimd_WS_weighted_mean_histo, mvimd_IO_weighted_mean_histo, mvimd_PO_weighted_mean_histo, mvimd_RS_weighted_mean_histo, mvimd_ABS_weighted_mean_histo]


PE_RS_weighted_mean_histo = mtpr_RS_weighted_mean_histo - (mer_RS_weighted_mean_histo)
PE_WS_weighted_mean_histo = mtpr_WS_weighted_mean_histo - (mer_WS_weighted_mean_histo)
PE_ABS_weighted_mean_histo = mtpr_ABS_weighted_mean_histo - (mer_ABS_weighted_mean_histo)
PE_IO_weighted_mean_histo = mtpr_IO_weighted_mean_histo - (mer_IO_weighted_mean_histo)
PE_PO_weighted_mean_histo = mtpr_PO_weighted_mean_histo - (mer_PO_weighted_mean_histo)
PE_all_weighted_mean_histo = mtpr_all_weighted_mean_histo - (mer_all_weighted_mean_histo)
PE_histo = [PE_all_weighted_mean_histo, PE_WS_weighted_mean_histo, PE_IO_weighted_mean_histo, PE_PO_weighted_mean_histo, PE_RS_weighted_mean_histo, PE_ABS_weighted_mean_histo]



mer_RS_weighted_mean = weight(mer_RS, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mer_WS_weighted_mean = weight(mer_WS, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mer_ABS_weighted_mean = weight(mer_ABS, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mer_IO_weighted_mean = weight(mer_IO, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mer_PO_weighted_mean = weight(mer_PO, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mer_all_weighted_mean = weight(mer_all, grid_cell_area_all_histo, total_area_of_earth_all_histo)
#mer = [-mer_RS_weighted_mean,-mer_WS_weighted_mean,-mer_ABS_weighted_mean,-mer_IO_weighted_mean,-mer_PO_weighted_mean,-mer_all_weighted_mean]
mer = [mer_all_weighted_mean,mer_WS_weighted_mean, mer_IO_weighted_mean, mer_PO_weighted_mean, mer_RS_weighted_mean, mer_ABS_weighted_mean]

mtpr_RS_weighted_mean = weight(mtpr_RS, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mtpr_WS_weighted_mean = weight(mtpr_WS, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mtpr_ABS_weighted_mean = weight(mtpr_ABS, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mtpr_IO_weighted_mean = weight(mtpr_IO, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mtpr_PO_weighted_mean = weight(mtpr_PO, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mtpr_all_weighted_mean = weight(mtpr_all, grid_cell_area_all_histo, total_area_of_earth_all_histo)
#mtpr = [mtpr_RS_weighted_mean, mtpr_WS_weighted_mean, mtpr_ABS_weighted_mean, mtpr_IO_weighted_mean, mtpr_PO_weighted_mean, mtpr_all_weighted_mean ]
mtpr = [mtpr_all_weighted_mean,mtpr_WS_weighted_mean, mtpr_IO_weighted_mean, mtpr_PO_weighted_mean, mtpr_RS_weighted_mean, mtpr_ABS_weighted_mean]

mvimd_RS_weighted_mean = weight(mvimd_RS_ssp, grid_cell_area_RS_histo, total_area_of_earth_RS_histo)
mvimd_WS_weighted_mean = weight(mvimd_WS_ssp, grid_cell_area_WS_histo, total_area_of_earth_WS_histo)
mvimd_ABS_weighted_mean = weight(mvimd_ABS_ssp, grid_cell_area_ABS_histo, total_area_of_earth_ABS_histo)
mvimd_IO_weighted_mean = weight(mvimd_IO_ssp, grid_cell_area_IO_histo, total_area_of_earth_IO_histo)
mvimd_PO_weighted_mean = weight(mvimd_PO_ssp, grid_cell_area_PO_histo, total_area_of_earth_PO_histo)
mvimd_all_weighted_mean = weight(mvimd_all_ssp, grid_cell_area_all_histo, total_area_of_earth_all_histo)
mvimd = [mvimd_all_weighted_mean,mvimd_WS_weighted_mean, mvimd_IO_weighted_mean, mvimd_PO_weighted_mean, mvimd_RS_weighted_mean, mvimd_ABS_weighted_mean]


PE_RS_weighted_mean = mtpr_RS_weighted_mean - (mer_RS_weighted_mean)
PE_WS_weighted_mean = mtpr_WS_weighted_mean - (mer_WS_weighted_mean)
PE_ABS_weighted_mean = mtpr_ABS_weighted_mean - (mer_ABS_weighted_mean)
PE_IO_weighted_mean = mtpr_IO_weighted_mean - (mer_IO_weighted_mean)
PE_PO_weighted_mean = mtpr_PO_weighted_mean - (mer_PO_weighted_mean)
PE_all_weighted_mean = mtpr_all_weighted_mean - (mer_all_weighted_mean)
#PE = [PE_RS_weighted_mean, PE_WS_weighted_mean, PE_ABS_weighted_mean, PE_IO_weighted_mean, PE_PO_weighted_mean, PE_all_weighted_mean]
PE = [PE_all_weighted_mean, PE_WS_weighted_mean, PE_IO_weighted_mean, PE_PO_weighted_mean, PE_RS_weighted_mean, PE_ABS_weighted_mean]



def difference(array1, array2):
    array_diff = np.zeros((6,432))
    for i in range(len(array1)):
        array_diff[i] = array1[i] - array2[i]
        """
        for j in range(432):
            if (array1[i][j] >= 0 and array2[i][j] >= 0): 
                array_diff[i,j] = np.abs(array1[i][j]) - np.abs(array2[i][j])
            elif (array1[i][j] <= 0 and array2[i][j] <= 0):
                array_diff[i,j] = np.abs(array1[i][j]) - np.abs(array2[i][j])
            elif array1[i][j] <  0 and array2[i][j] > 0:
                array_diff[i,j] = np.abs(array1[i][j] - array2[i][j])
            elif array1[i][j] > 0 and array2[i][j] < 0:
                array_diff[i,j] = array1[i][j] - array2[i][j]
        """        
    return array_diff

#on veut la différence entre ssp et histo (genre il y a xxxmm/month en moins ou en plus demoisutre convergence)


mer_diff = difference(mer, mer_histo)
mtpr_diff = difference(mtpr, mtpr_histo)
PE_diff = difference(PE, PE_histo)
mvimd_diff = difference(mvimd, mvimd_histo)


#sort by month ===============================================================
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

"""
mer_sort_month_histo = np.zeros((6, 12))
mtpr_sort_month_histo = np.zeros((6, 12))
mvimd_sort_month_histo = np.zeros((6, 12))
PE_sort_month_histo = np.zeros((6, 12))
mer_std_histo = np.zeros((6,12))
mtpr_std_histo = np.zeros((6,12))
mvimd_std_histo = np.zeros((6,12))
PE_std_histo = np.zeros((6,12))
"""
mer_sort_month = np.zeros((6, 12))
mtpr_sort_month = np.zeros((6, 12))
mvimd_sort_month = np.zeros((6, 12))
PE_sort_month = np.zeros((6, 12))
mer_std = np.zeros((6,12))
mtpr_std = np.zeros((6,12))
mvimd_std = np.zeros((6,12))
PE_std = np.zeros((6,12))

for i in range(len(mer_histo)):
    """
    mer_sort_month_histo[i], mer_std_histo[i] = repartion_by_months(mer_histo[i],0,432,1) 
    mtpr_sort_month_histo[i], mtpr_std_histo[i] = repartion_by_months(mtpr_histo[i],0,432,1)
    #mvimd_sort_month_histo[i], mvimd_std_histo[i] = repartion_by_months(mvimd_histo[i],0,432,1)
    PE_sort_month_histo[i] = (np.abs(mtpr_sort_month_histo[i]) - np.abs(mer_sort_month_histo[i]))
    PE_std_histo[i] = np.sqrt((mtpr_std_histo[i])**2 + (mer_std_histo[i])**2)  #pas sur que ca soit bon ca devrait etre la somme ou autre chose 
    """
    mer_sort_month[i], mer_std[i] = repartion_by_months(mer_diff[i],0,432,1) 
    mtpr_sort_month[i], mtpr_std[i] = repartion_by_months(mtpr_diff[i],0,432,1)
    mvimd_sort_month[i], mvimd_std[i] = repartion_by_months(mvimd_diff[i],0,432,1)
    PE_sort_month[i], PE_std[i] = repartion_by_months(PE_diff[i], 0, 432, 1)
    
    #PE_sort_month[i] = (np.abs(mtpr_sort_month[i]) - np.abs(mer_sort_month[i]))
    #PE_std[i] = np.sqrt((mtpr_std[i])**2 + (mer_std[i])**2)  #pas sur que ca soit bon ca devrait etre la somme ou autre chose 


lst= []
lst.append("mer")
for i in range(6):
    lst.append(np.mean(mer_sort_month[i]))
lst.append("mtpr")
for i in range(6):
    lst.append(np.mean(mtpr_sort_month[i]))
lst.append("PE")
for i in range(6):
    lst.append(np.mean(PE_sort_month[i]))
lst.append("mvimd")
for i in range(6):
    lst.append(np.mean(mvimd_sort_month[i]))    
print(lst)



#Figures====================================================================================================
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
#label_sectors = ["RS", "WS", "ABS", "IO", "PO", "all"]
label_sectors = ["all", "WS", "IO", "PO", "RS", "ABS"]
#Figure avec courbes saisonnières--------------------------------------------------------------------------


fig = plt.figure(figsize = (16,10) , dpi = 150)
for i in range(len(mer_histo)):
    ax = plt.subplot(231+i)
    ax.plot(months, mer_sort_month[i],  label = "E", linestyle = "dashed")     #positif
    ax.plot(months, mtpr_sort_month[i],  label = "P", linestyle = "dashed")    #neg
    ax.plot(months, PE_sort_month[i],  label = "P-E", linestyle = "dashed")
    ax.plot(months, mvimd_sort_month[i],  label = "MC", linestyle = "dashed")
    #ax.errorbar(months, mvimd_sort_month_histo[i], yerr = mvimd_std_histo[i],  label = "MC", linestyle = "dashed")          #(positif)
    ax.grid()
    """
    if i == 1:
        ax.yaxis.set_ticks([-10,0,10,20,30,40,50,60])
    elif i == 4:
        ax.yaxis.set_ticks([0,10,20,30,40,50,60])
    else:
        ax.yaxis.set_ticks([-9,-6,-3,0,3,6,9,12,15,18,21,24])
    """    
        
    #ax.yaxis.set_major_locator(plt.LinearLocator())
    ax.set_title(label_sectors[i])
    ax.set_ylabel("[$mm/month$]")
    ax.set_xlabel("Month")
    if i == 2:
        ax.legend()
        
        
        
        

fig = plt.figure(figsize = (8,8) , dpi = 150)
plt.plot(months, mer_sort_month[0],  label = "E", linestyle = "dashed")     #positif
plt.plot(months, mtpr_sort_month[0],  label = "P", linestyle = "dashed")    #neg
plt.plot(months, PE_sort_month[0],  label = "P-E", linestyle = "dashed")
plt.plot(months, mvimd_sort_month[0],  label = "MC", linestyle = "dashed")
#ax.errorbar(months, mvimd_sort_month_histo[i], yerr = mvimd_std_histo[i],  label = "MC", linestyle = "dashed")          #(positif)
plt.grid()
 
plt.ylabel("[$mm/month$]")
plt.xlabel("Month")
plt.legend()        