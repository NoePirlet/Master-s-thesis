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
from matplotlib.offsetbox import AnchoredText

sec_year = 31536000
sec_month = 2592000


access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:512,0,:,:])*sec_year                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
dset.close()                                                                   #2592000 si on est en [kg/((m^2)*month)]

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:512,0,:,:])*sec_year
dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:512,0,:,:])*sec_year
dset.close()

access_pr_file = 'Total_column_water.nc'
dset = nc.Dataset(access_pr_file)
tcw = np.array(dset['tcw'][0:512,0,:,:])
dset.close()

#negative values indicate evaporation and positive values indicate deposition.
access_pr_file = 'Mean_snow_evaporation_rate.nc'    
dset = nc.Dataset(access_pr_file)
mser = np.array(dset['mser'][0:512,0,:,:])*sec_year      
dset.close()

year = 1997 #année avec laquelle on veut travailler

def year_to_month (year):
    temp = year-1978
    nbr_of_month = 12*temp
    start_month = nbr_of_month-12
    end_month = nbr_of_month
    year_vec = [start_month,end_month]
    return year_vec

year_to_months = year_to_month(year)


def year_array (array, year_to_months):
        new_array = array[year_to_months[0]:year_to_months[1],:,:]
        return new_array

short_mer = year_array(mer, year_to_months)
short_mtpr = year_array(mtpr, year_to_months)
short_mvimd = year_array(mvimd, year_to_months)
short_tcw = year_array(tcw, year_to_months)
short_mser = year_array(mser, year_to_months)


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



def sum_year(array):
    somme_year = np.zeros((42,121,1440))
    j = -1
    for i in range(0,504):            #somme 1979-fin 2020
        if (i)%12 == 0:
            j += 1
        somme_year[j,:,:] += (array[i,:,:])/12
    
    return somme_year

"""
Weddell Sea                         48,336
Indian Ocean                        336,624
Pacific Ocean                       624,912
Ross Sea                            912,1200
Amundsen Bellingshausen Seas        1200,1488
"""

mer_by_year = sum_year(mer)     
mtpr_by_year = sum_year(mtpr)
mvimd_by_year = sum_year(mvimd)
mser_by_year = sum_year(mser)


def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(len(array))
    for i in range(0,len(array)):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])     #[kg/year]
    return weighted_mean    
    
    
mer_weighted_mean = weight(mer_by_year, grid_cell_area, total_area_of_earth)
mtpr_weighted_mean = weight(mtpr_by_year, grid_cell_area, total_area_of_earth)
mvimd_weighted_mean = weight(mvimd_by_year, grid_cell_area, total_area_of_earth)
mser_weighted_mean = weight(mser_by_year, grid_cell_area, total_area_of_earth)
PE_weighted_mean = mtpr_weighted_mean - (-mer_weighted_mean)


short_mer_weighted_mean = weight(short_mer, grid_cell_area, total_area_of_earth)
short_mtpr_weighted_mean = weight(short_mtpr, grid_cell_area, total_area_of_earth)
short_mvimd_weighted_mean = weight(short_mvimd, grid_cell_area, total_area_of_earth)
short_mser_weighted_mean = weight(short_mser, grid_cell_area, total_area_of_earth)
short_PE_weighted_mean = short_mtpr_weighted_mean - (-short_mer_weighted_mean)

    
#Figures====================================================================================================
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

#Figure avec courbes saisonnières--------------------------------------------------------------------------
"""
fig = plt.figure(dpi = 100)
plt.plot(months, -short_mer_weighted_mean , label = "Evaporation (E)", linestyle = "dashed")     
plt.plot(months, short_mtpr_weighted_mean , label = "Precipitation (P)", linestyle = "dashed")    
plt.plot(months, -short_mvimd_weighted_mean , label = "Moisture Convergence (MC)", linestyle = "dashed")          
plt.plot(months, -short_mser_weighted_mean , label = "Snow evaporation (SE)", linestyle = "dashed")     #on met le mois car de base neg = evap donc pour avoir un terme qui ajoute de l'eau dans le budget
plt.plot(months, short_PE_weighted_mean , label = "P - E", linestyle = "dashed")
plt.grid()
plt.ylabel("[$kg/month$]")  
plt.xlabel("Months")
plt.title("Seasonal curve for each process in the budget (1997)")
plt.legend()
plt.show()

"""

#figure avec répartion budget---------------------------------------------------------------------------

somme_right = ((-short_mer_weighted_mean) - short_mtpr_weighted_mean) -(short_mvimd_weighted_mean)

somme_right_std = np.std(somme_right)
"""
fig = plt.figure(dpi = 100)
ax = plt.subplot(111)

ax.bar(months, -short_mer_weighted_mean , label = "Evaporation (E)")     #Dans ma colonne j'ai un apport d'eau par evaporation ==> terme positif
ax.bar(months, -short_mtpr_weighted_mean , label = "Precipitation (P)")  #Dans ma colonne j'ai une perte d'eau par précipitation ==> terme négatif
ax.bar(months, -short_mvimd_weighted_mean , bottom = -short_mer_weighted_mean , label = "Moisture Convergence (MC)")   #Dans ma colonne j'ai soit convergence soit divergence d'humidité donc respectivement terme positif et negatif
#ax.plot(months, somme_right, label = "E-P-MC", color = "Red" )
ax.errorbar(months, somme_right, yerr = somme_right_std,  label = "E-P-MC", color = "red")
ax.set_ylabel("[$mm/month$]")
ax.set_xlabel("Months")
ax.set_title("Relative importance of each process in the atmospheric freshwater budget (1997)", loc = 'center')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.show()  

"""

#Figure avec multi-year trend-----------------------------------------------------------------------------

def years(nbr, start_year):
    years = np.zeros(nbr)
    for i in range(nbr):
        years[i] = start_year+i
    return years
    
    
def lin_reg (x, y):
    linreg = linregress(x, y)
    m = linreg.slope
    b = linreg.intercept
    p = linreg.pvalue
    return(m,b,p)


mer_linreg = lin_reg(years(42,1979), -mer_weighted_mean)
mtpr_linreg = lin_reg(years(42,1979), mtpr_weighted_mean)
mvimd_linreg =  lin_reg(years(42,1979), -mvimd_weighted_mean)
mser_linreg = lin_reg(years(42,1979), -mser_weighted_mean)
PE_linreg = lin_reg(years(42,1979), PE_weighted_mean)


fig = plt.figure(figsize = (8,10), dpi = 250)
#fig.suptitle("All sectors")
ax1 = plt.subplot(411)
ax1.plot(years(42,1979), -mer_weighted_mean , linestyle = "dashed")  #neg term
#plt.plot(years(42,1979), mer_linreg[0]*years(42,1979) + mer_linreg[1], label= "m = {:.2f}, p = {:.3f}".format(mer_linreg[0], mer_linreg[2]))
ax1.set_ylabel("[$mm/year$]")   #km3/year
ax1.set_xlabel("Years")
ax1.set_xticks([1980,1985,1990,1995,2000,2005,2010,2015,2020])
ax1.set_title("Evaporation")
ax1.grid()

 
ax2 = plt.subplot(412)
ax2.plot(years(42,1979), mtpr_weighted_mean , linestyle = "dashed")
ax2.plot(years(42,1979), mtpr_linreg[0]*years(42,1979) + mtpr_linreg[1])
ax2.set_ylabel("[$mm/year$]")   #km3/year
ax2.set_xlabel("Years") 
ax2.set_xticks([1980,1985,1990,1995,2000,2005,2010,2015,2020])
ax2.set_title("Precipitation")
at = AnchoredText(
    "Linear trend : {:.2f}".format(mtpr_linreg[0]), prop=dict(size=12, color = "darkorange"), frameon=True, loc='lower right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax2.add_artist(at)
#ax2.text(1,0,"m = {:.2f}".format(mtpr_linreg[0]), fontsize=18, color='k')
ax2.grid()



ax3 = plt.subplot(413)
ax3.plot(years(42,1979), -mvimd_weighted_mean , linestyle = "dashed")
ax3.plot(years(42,1979), mvimd_linreg[0]*years(42,1979) + mvimd_linreg[1])#, label= "m = {:.2f}, p = {:.3f}".format(mvimd_linreg[0], mvimd_linreg[2]))
ax3.set_ylabel("[$mm/year$]")   #km3/year 
ax3.set_xticks([1980,1985,1990,1995,2000,2005,2010,2015,2020])
ax3.set_xlabel("Years")
ax3.set_title("Moisture Convergence")
at = AnchoredText(
    "Linear trend : {:.2f}".format(mvimd_linreg[0]), prop=dict(size=12, color = "darkorange"), frameon=True, loc='lower right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax3.add_artist(at)
#ax3.text(1,0,"m = {:.2f}".format(mvimd_linreg[0]), fontsize=18, color='k')
ax3.grid()


#plt.subplot(514)
#plt.xticks([])
#plt.plot(years(42,1979), -mser_weighted_mean , linestyle = "dashed")
#plt.plot(years(42,1979), mser_linreg[0]*years(42,1979) + mser_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(mser_linreg[0],mser_linreg[1], mser_linreg[2]))
#plt.ylabel("[$mm/year$]")   #km3/year 
#plt.xlabel("Years")
#plt.title("Snow Evaporation")
#plt.legend()

ax4 = plt.subplot(414)
ax4.plot(years(42,1979), PE_weighted_mean , linestyle = "dashed")
ax4.plot(years(42,1979), PE_linreg[0]*years(42,1979) + PE_linreg[1])#, label= "m = {:.2f}, p = {:.3f}".format(PE_linreg[0], PE_linreg[2]))
ax4.set_ylabel("[$mm/year$]")   #km3/year
ax4.set_xlabel("Years")
ax4.set_xticks([1980,1985,1990,1995,2000,2005,2010,2015,2020])
ax4.set_title("P - E")
at = AnchoredText(
    "Linear trend : {:.2f}".format(PE_linreg[0]), prop=dict(size=12, color = "darkorange"), frameon=True, loc='lower right')
at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax4.add_artist(at)
#ax4.text(1,0,"m = {:.2f}".format(PE_linreg[0]), fontsize=18, color='k')
ax4.grid()


fig.tight_layout()
plt.show()    

##All in one--------------------------------------------------------------------
"""
fig = plt.figure(figsize = (6,6), dpi = 500)


plt.plot(years(42,1979), -mer_weighted_mean , linestyle = "dashed", label = "Evaporation")  #neg term
plt.plot(years(42,1979), mer_linreg[0]*years(42,1979) + mer_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(mer_linreg[0],mer_linreg[1], mer_linreg[2]))
plt.plot(years(42,1979), mtpr_weighted_mean , linestyle = "dashed", label = "Precipitation")
plt.plot(years(42,1979), mtpr_linreg[0]*years(42,1979) + mtpr_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(mtpr_linreg[0],mtpr_linreg[1], mtpr_linreg[2]))
plt.plot(years(42,1979), -mvimd_weighted_mean , linestyle = "dashed", label = "Moisture Convergence")
plt.plot(years(42,1979), mvimd_linreg[0]*years(42,1979) + mvimd_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(mvimd_linreg[0],mvimd_linreg[1], mvimd_linreg[2]))
plt.plot(years(42,1979), PE_weighted_mean , linestyle = "dashed", label = "P - E")
plt.plot(years(42,1979), PE_linreg[0]*years(42,1979) + PE_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(PE_linreg[0],PE_linreg[1], PE_linreg[2]))
plt.ylabel("[$mm/year$]")   
plt.xlabel("Years")
plt.grid()
plt.legend()

#plt.subplot(514)
#plt.xticks([])
#plt.plot(years(42,1979), -mser_weighted_mean , linestyle = "dashed")
#plt.plot(years(42,1979), mser_linreg[0]*years(42,1979) + mser_linreg[1], label= "a = {:.2f}, b = {:.2f}, p = {:.3f}".format(mser_linreg[0],mser_linreg[1], mser_linreg[2]))
#plt.ylabel("[$mm/year$]")   #km3/year 
#plt.xlabel("Years")
#plt.title("Snow Evaporation")
#plt.legend()

fig.tight_layout()
plt.show()  
"""