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

year = 1979 #année avec laquelle on veut travailler

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
"""
def sum_area (new_array):
    somme = np.zeros(12)
    for i in range(12):
        somme[i] = np.sum(new_array[i,:,:])
    return somme
        
#surface below 60°S
r = 6371000  #[m]
surface = 2*np.pi*(r**2)*(1-np.cos(30))   #theta = 90-60 = 30
#h = r - r*np.cos(60)
#surface = np.pi*h*(2*r-h)  #[m^2]
"""
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
    somme_year = np.zeros(42)
    j = -1
    for i in range(0,504):            #somme 1979-fin 2020
        if (i)%12 == 0:
            j += 1
            if j == 43:
                break
        somme_year[j] += np.sum(array[i,:,:]) 
    return somme_year

mer_by_year = sum_year(mer)

def weight(array, grid_cell_area, total_area_of_earth):
    weighted_mean = np.zeros(42)
    for i in range(0,42):
        weighted_mean[i] = ((array[i,:,:] * grid_cell_area) / total_area_of_earth).sum(['latitude','longitude'])
    return weighted_mean    
    
    
mer_weighted_mean = weight(mer, grid_cell_area, total_area_of_earth)
mtpr_weighted_mean = weigth(mtpr, grid_cell_area, total_area_of_earth)
mvimd_m


#---------------------------------------------------------------------------------------
somme_mer = sum_area(short_mer)*surface   #[kg/month] 
somme_mtpr = sum_area(short_mtpr)*surface
somme_mvimd = sum_area(short_mvimd)*surface 
somme_mser = sum_area(short_mser)*surface    
somme_tcw = sum_area(short_tcw)*surface   #[kg]
somme_tcw_time = somme_tcw/2592000        #[kg/month]

  
km3 = 1/10e12   #conversion en km3


somme_PE = (np.abs(somme_mtpr) - np.abs(somme_mer))*km3
somme_right = (-somme_PE - (somme_mvimd)*km3)

#-----------------------------------------------------------------------------------------------------
def sum_year(array):
    somme_year = np.zeros(42)
    j = -1
    for i in range(0,504):            #somme 1979-fin 2020
        if (i)%12 == 0:
            j += 1
            if j == 43:
                break
        somme_year[j] += np.sum(array[i,:,:]) 
    return somme_year

somme_year_mer = sum_year(mer)*surface
somme_year_mtpr = sum_year(mtpr)*surface
somme_year_mvimd = sum_year(mvimd)*surface

somme_year_PE = (np.abs(somme_year_mtpr) - np.abs(somme_year_mer))*km3
        
#Figures====================================================================================================
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

#Figure avec courbes saisonnières--------------------------------------------------------------------------
fig = plt.figure(dpi = 100)
plt.plot(months, -(somme_mer)*km3, label = "Evaporation (E)", linestyle = "dashed")     #positif
plt.plot(months, -(somme_mtpr)*km3, label = "Precipitation (P)", linestyle = "dashed")    #neg
plt.plot(months, -(somme_mvimd)*km3, label = "Moisture Convergence (MC)", linestyle = "dashed")          #(positif)
plt.plot(months, -somme_mser*km3, label = "Snow evaporation (SE)", linestyle = "dashed")     #on met le mois car de base neg = evap donc pour avoir un terme qui ajoute de l'eau dans le budget
plt.plot(months, somme_PE, label = "P - E", linestyle = "dashed")
plt.grid()
plt.ylabel("Volume [${km^3}/month$]")  
plt.xlabel("Months")
plt.title("Seasonal curve for each process in the budget")
plt.legend()
plt.show()

#figure avec le résidu du budget------------------------------------------------------------------------
"""
fig = plt.figure(dpi = 100)
plt.plot(months, (somme_tcw_time)*km3, label = "Water column", linestyle = "dashed") 
#plt.plot(months, somme_right, label = "Budget", linestyle = "dashed")
#plt.plot(months, (somme_tcw_time*km3)-(somme_right), label = "diff both sides" )
plt.grid()
plt.ylabel("Volume [km3]")   #km3/month !!!!! and km3!!!!!
plt.legend()
plt.show()
"""
#figure avec répartion budget---------------------------------------------------------------------------

fig = plt.figure(dpi = 100)
plt.bar(months, -(somme_mer)*km3, label = "Evaporation (E)")     #Dans ma colonne j'ai un apport d'eau par evaporation ==> terme positif
plt.bar(months, -(somme_mtpr)*km3, label = "Precipitation (P)")  #Dans ma colonne j'ai une perte d'eau par précipitation ==> terme négatif
plt.bar(months, -(somme_mvimd)*km3, bottom = -(somme_mer)*km3, label = "Moisture Convergence (MC)")   #Dans ma colonne j'ai soit convergence soit divergence d'humidité donc respectivement terme positif et negatif
plt.plot(months, somme_right, label = "E-P-MC", color = "Red" )
plt.ylabel("Volume [${km^3}/month$]")
plt.xlabel("Months")
plt.title("Relative importance of each process in the freshwater budget")
plt.legend()
plt.show()

#Figure avec multi-year trend-----------------------------------------------------------------------------

def years(nbr, start_year):
    years = np.zeros(nbr)
    for i in range(nbr):
        years[i] = start_year+i
    return years
    
    
def lin_reg (x, y):
    m, b = np.polyfit(x, y, 1)
    return(m,b)

fig = plt.figure(figsize = (8,10), dpi = 500)

plt.subplot(411)
plt.xticks([])
plt.plot(years(42,1979), np.abs(somme_year_mer)*km3, linestyle = "dashed")
plt.plot(years(42,1979), lin_reg(years(42,1979), np.abs(somme_year_mer)*km3)[0]*years(42,1979) + lin_reg(years(42,1979), np.abs(somme_year_mer)*km3)[1], label= "a = {0}, b = {1}".format(int(lin_reg(years(42,1979), np.abs(somme_year_mer)*km3)[0]),int(lin_reg(years(42,1979), np.abs(somme_year_mer)*km3)[1]) ))
plt.legend()
plt.ylabel("Volume [${km^3}/year$]")   #km3/year
plt.xlabel("Years")
plt.title("Evaporation")

 
plt.subplot(412)
plt.xticks([])
plt.plot(years(42,1979), np.abs(somme_year_mtpr)*km3, linestyle = "dashed")
plt.plot(years(42,1979), lin_reg(years(42,1979), np.abs(somme_year_mtpr)*km3)[0]*years(42,1979) + lin_reg(years(42,1979), np.abs(somme_year_mtpr)*km3)[1], label= "a = {0}, b = {1}".format(int(lin_reg(years(42,1979), np.abs(somme_year_mtpr)*km3)[0]),int(lin_reg(years(42,1979), np.abs(somme_year_mtpr)*km3)[1]) ))
plt.ylabel("Volume [${km^3}/year$]")   #km3/year
plt.xlabel("Years") 
plt.title("Precipitation")
plt.legend()


plt.subplot(413)
plt.xticks([])
plt.plot(years(42,1979), np.abs(somme_year_mvimd)*km3, linestyle = "dashed")
plt.plot(years(42,1979), lin_reg(years(42,1979), np.abs(somme_year_mvimd)*km3)[0]*years(42,1979) + lin_reg(years(42,1979), np.abs(somme_year_mvimd)*km3)[1], label= "a = {0}, b = {1}".format(int(lin_reg(years(42,1979), np.abs(somme_year_mvimd)*km3)[0]),int(lin_reg(years(42,1979), np.abs(somme_year_mvimd)*km3)[1]) ))
plt.ylabel("Volume [${km^3}/year$]")   #km3/year 
plt.xlabel("Years")
plt.title("Moisture Convergence")
plt.legend()

plt.subplot(414)
plt.plot(years(42,1979), somme_year_PE, linestyle = "dashed")
plt.plot(years(42,1979), lin_reg(years(42,1979), np.abs(somme_year_PE))[0]*years(42,1979) + lin_reg(years(42,1979), np.abs(somme_year_PE))[1], label= "a = {0}, b = {1}".format(int(lin_reg(years(42,1979), np.abs(somme_year_PE))[0]),int(lin_reg(years(42,1979), np.abs(somme_year_PE))[1]) ))
plt.ylabel("Volume [${km^3}/year$]")   #km3/year
plt.xlabel("Years")
plt.title("P - E")
plt.legend()



fig.tight_layout()
plt.show()    
    

    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
    
  
        
        
