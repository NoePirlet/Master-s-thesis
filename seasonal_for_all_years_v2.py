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


mer_sort_month, mer_std = repartion_by_months(weight(mer, grid_cell_area, total_area_of_earth),0,504,1) 
mtpr_sort_month, mtpr_std = repartion_by_months(weight(mtpr, grid_cell_area, total_area_of_earth),0,504,1)
mvimd_sort_month, mvimd_std = repartion_by_months(weight(mvimd, grid_cell_area, total_area_of_earth),0,504,1)
mser_sort_month, mser_std = repartion_by_months(weight(mser, grid_cell_area, total_area_of_earth),0,504,1)
tcw_sort_month, tcw_std = repartion_by_months(weight(tcw, grid_cell_area, total_area_of_earth),0,504,1)

PE_sort_month = (np.abs(mtpr_sort_month) - np.abs(mer_sort_month))
PE_std = np.sqrt((mtpr_std)**2 + (mer_std)**2)
def convert(lst):
    return [ -i for i in lst ]
#Figures====================================================================================================
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

#Figure avec courbes saisonnières--------------------------------------------------------------------------
"""
fig = plt.figure(dpi = 100)
plt.plot(months, convert(mer_sort_month), label = "Evaporation (E)", linestyle = "dashed")     #positif
plt.plot(months, convert(mtpr_sort_month), label = "Precipitation (P)", linestyle = "dashed")    #neg
plt.plot(months, convert(mvimd_sort_month), label = "Moisture Convergence (MC)", linestyle = "dashed")          #(positif)
plt.plot(months, convert(mser_sort_month), label = "Snow evaporation (SE)", linestyle = "dashed")     #on met le mois car de base neg = evap donc pour avoir un terme qui ajoute de l'eau dans le budget
plt.plot(months, PE_sort_month, label = "P - E", linestyle = "dashed")
plt.grid()
plt.ylabel("Volume [${km^3}/month$]")  
plt.xlabel("Months")
plt.title("Seasonal mean for each month (1979-2020)")
plt.legend()
plt.show()
"""

fig = plt.figure(dpi = 100)
ax = plt.subplot(111)

ax.errorbar(months, convert(mer_sort_month), yerr = mer_std,  label = "Evaporation (E)", linestyle = "dashed")     #positif
ax.errorbar(months, (mtpr_sort_month), yerr = mtpr_std,  label = "Precipitation (P)", linestyle = "dashed")    #neg
ax.errorbar(months, convert(mvimd_sort_month), yerr = mvimd_std,  label = "Moisture Convergence (MC)", linestyle = "dashed")          #(positif)
ax.errorbar(months, convert(mser_sort_month), yerr = mser_std, label = "Snow evaporation (SE)", linestyle = "dashed")     #on met le mois car de base neg = evap donc pour avoir un terme qui ajoute de l'eau dans le budget
ax.errorbar(months, PE_sort_month, yerr = PE_std,  label = "P - E", linestyle = "dashed")
ax.grid()
ax.set_title("Seasonal mean for each month (1979-2020)")
ax.set_ylabel("[${kg}/month$]")
ax.set_xlabel("Months")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()  

#figure avec répartion budget---------------------------------------------------------------------------
somme_right = np.zeros(12)
somme_right_std = np.zeros(12)
for i in range(12):
    somme_right[i] = ((-mer_sort_month[i]) - mtpr_sort_month[i] ) -(mvimd_sort_month[i])
    somme_right_std[i] = np.sqrt(mer_std[i]**2 + mtpr_std[i]**2 + mvimd_std[i]**2)    #en considérant les termes indépendants

fig = plt.figure(dpi = 500)
ax = plt.subplot(111)

ax.bar(months, convert(mer_sort_month) , label = "Evaporation (E)")     #Dans ma colonne j'ai un apport d'eau par evaporation ==> terme positif
ax.bar(months, convert(mtpr_sort_month) , label = "Precipitation (P)")  #Dans ma colonne j'ai une perte d'eau par précipitation ==> terme négatif
ax.bar(months, convert(mvimd_sort_month) , bottom = convert(mer_sort_month) , label = "Moisture \nConvergence (MC)")   #Dans ma colonne j'ai soit convergence soit divergence d'humidité donc respectivement terme positif et negatif
#ax.plot(months, somme_right, label = "E-P-MC (Total flux = {:.2f}mm/month)".format(np.sum(somme_right)), color = "Red" )
ax.errorbar(months, somme_right, yerr = somme_right_std,  label = "E-P-MC and standard \ndeviation (95%)", color = "red")
ax.set_ylabel("[$mm/month$]")
ax.set_xlabel("Months")
#ax.set_title("Relative importance of each process in the atmospheric freshwater budget (1979-2020)", loc = 'center')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()    

#calcul percentage----------------------------------------------------------------------------------------
abs_sum = np.zeros(12)
final_lst = []
for i in range(12):
    abs_sum = abs(mer_sort_month[i]) + abs(mtpr_sort_month[i]) + abs(mvimd_sort_month[i])
    percentage = (abs(mer_sort_month[i])/abs_sum, abs(mtpr_sort_month[i])/abs_sum, abs(mvimd_sort_month[i])/abs_sum)
    final_lst.append(percentage)
print(final_lst)
















