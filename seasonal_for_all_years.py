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


def sum_area (array):
    #surface below 60°S
    r = 6371000  #[m]
    h = r - r*np.cos(60)
    surface = np.pi*h*(2*r-h)  #[m^2]
    somme = np.zeros(512)
    for i in range(512):
        somme[i] = np.sum(array[i,:,:])*surface
    return somme   
        

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
            jan_lst.append(sum_area(array)[i])
        if (i+1)%12 == 2:
            feb_lst.append(sum_area(array)[i])
        if (i+1)%12 == 3:
            mar_lst.append(sum_area(array)[i])
        if (i+1)%12 == 4:
            apr_lst.append(sum_area(array)[i])
        if (i+1)%12 == 5:
            may_lst.append(sum_area(array)[i])
        if (i+1)%12 == 6:
            jun_lst.append(sum_area(array)[i])
        if (i+1)%12 == 7:
            jul_lst.append(sum_area(array)[i])
        if (i+1)%12 == 8:
            aug_lst.append(sum_area(array)[i])
        if (i+1)%12 == 9:
            sep_lst.append(sum_area(array)[i])
        if (i+1)%12 == 10:
            octb_lst.append(sum_area(array)[i])
        if (i+1)%12 == 11:
            nov_lst.append(sum_area(array)[i])
        if (i+1)%12 == 0:
            dec_lst.append(sum_area(array)[i])
            
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


km3 = 1/10e12   #conversion en km3

mer_sort_month, mer_std = repartion_by_months(mer,0,504,km3) 
mtpr_sort_month, mtpr_std = repartion_by_months(mtpr,0,504,km3)
mvimd_sort_month, mvimd_std = repartion_by_months(mvimd,0,504,km3)
mser_sort_month, mser_std = repartion_by_months(mser,0,504,km3)


PE_sort_month = (np.abs(mtpr_sort_month) - np.abs(mer_sort_month))
PE_std = mtpr_std - mer_std
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

#ax.errorbar(months, convert(mer_sort_month), yerr = mer_std,  label = "Evaporation (E)", linestyle = "dashed")     #positif
#ax.errorbar(months, convert(mtpr_sort_month), yerr = mtpr_std,  label = "Precipitation (P)", linestyle = "dashed")    #neg
#ax.errorbar(months, convert(mvimd_sort_month), yerr = mvimd_std,  label = "Moisture Convergence (MC)", linestyle = "dashed")          #(positif)
ax.errorbar(months, convert(mser_sort_month), yerr = mser_std, label = "Snow evaporation (SE)", linestyle = "dashed")     #on met le mois car de base neg = evap donc pour avoir un terme qui ajoute de l'eau dans le budget
#ax.errorbar(months, PE_sort_month, yerr = PE_std,  label = "P - E", linestyle = "dashed")
ax.grid()
ax.set_title("Seasonal mean for each month (1979-2020)")
ax.set_ylabel("Volume [${km^3}/month$]")
ax.set_xlabel("Months")
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
