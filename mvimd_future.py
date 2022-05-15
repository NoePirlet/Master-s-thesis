# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 17:21:09 2022

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


access_pr_file = 'intuaw_Emon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
lat_histo=np.array(dset['lat'][0:26])
lon_all_histo=np.array(dset['lon'][:])
"""
lon_all_histo=np.append(lon_all_histo, 360)
"""
umvimd_all_histo_v2 = np.array(dset['intuaw'][1548:1980,0:26,:])*sec_month
""" 
umvimd_all_histo_v2 = np.zeros((432,26,193))
for i in range(432):
    for j in range(25):
            umvimd_all_histo_v2[i,j] = np.append(umvimd_all_histo[i,j],(umvimd_all_histo[i,j,0]+umvimd_all_histo[i,j,-1])/2)
"""            
dset.close()                                                                  

access_pr_file = 'intvaw_Emon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
vmvimd_all_histo_v2 = np.array(dset['intvaw'][1548:1980,0:26,:])*sec_month
""" 
vmvimd_all_histo_v2 = np.zeros((432,26,193))
for i in range(432):
    for j in range(26):
            vmvimd_all_histo_v2[i,j] = np.append(vmvimd_all_histo[i,j],(vmvimd_all_histo[i,j,0]+vmvimd_all_histo[i,j,-1])/2)
"""            
dset.close()



access_pr_file = 'intuaw_Emon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc'
dset = nc.Dataset(access_pr_file)
umvimd_all_ssp_v2 = np.array(dset['intuaw'][600:1032,0:26,:])*sec_month
""" 
umvimd_all_ssp_v2 = np.zeros((432,26,193))
for i in range(432):
    for j in range(26):
            umvimd_all_ssp_v2[i,j] = np.append(umvimd_all_ssp[i,j],(umvimd_all_ssp[i,j,0]+umvimd_all_ssp[i,j,-1])/2)
"""        
dset.close()                                                                  

access_pr_file = 'intvaw_Emon_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_201501-210012.nc'
dset = nc.Dataset(access_pr_file)
vmvimd_all_ssp_v2 = np.array(dset['intvaw'][600:1032,0:26,:])*sec_month 
"""
vmvimd_all_ssp_v2 = np.zeros((432,26,193))
for i in range(432):
    for j in range(26):
            vmvimd_all_ssp_v2[i,j] = np.append(vmvimd_all_ssp[i,j],(vmvimd_all_ssp[i,j,0]+vmvimd_all_ssp[i,j,-1])/2)
"""            
dset.close()

"""
access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
lat_era=np.array(dset['latitude'][:])
lon_era=np.array(dset['longitude'][:])
mvimd_era = np.array(dset['mvimd'][0:432,0,:,:])*sec_month*(-1)               #*sec_month*(-1)       #*sec_year
dset.close()
"""
"""
access_pr_file = 'Vertical integral of eastward water vapour flux.nc'
dset = nc.Dataset(access_pr_file)
lat_era=np.array(dset['latitude'][:])
lon_era=np.array(dset['longitude'][:])
umvimd_era = np.array(dset['p71.162'][0:432,0,:,:])*sec_month              #*sec_month*(-1)       #*sec_year
dset.close()

access_pr_file = 'Vertical integral of northward water vapour flux.nc'
dset = nc.Dataset(access_pr_file)
vmvimd_era = np.array(dset['p72.162'][0:432,0,:,:])*sec_month              #*sec_month*(-1)       #*sec_year
dset.close()
"""



#def dx et dy ===============================================================================


mvimd_v0 = np.zeros((432,25,192))
mvimd_ssp_v0 = np.zeros((432,25,192)) 

from numpy import meshgrid, deg2rad, gradient, cos
R = 6371000   
"""
for i in range(len(lat_histo)):
    for j in range(len(lon_all_histo)):
        for k in range(432):  
                mvimd[k,i,j] = (96/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))*(umvimd_all_histo_v2[k,i,(j+1)%193]-umvimd_all_histo_v2[k,i,j])\
                    + (144/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))\
                    *(vmvimd_all_histo_v2[k,(i-1),j]*np.cos(deg2rad(lat_histo[i-1]))-vmvimd_all_histo_v2[k,i,j]*np.cos(deg2rad(lat_histo[i])))
                        
                mvimd_ssp[k,i,j] = (96/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))*(umvimd_all_ssp_v2[k,i,(j+1)%193]-umvimd_all_ssp_v2[k,i,j])\
                    + (144/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))\
                    *(vmvimd_all_ssp_v2[k,(i-1),j]*np.cos(deg2rad(lat_histo[i-1]))-vmvimd_all_ssp_v2[k,i,j]*np.cos(deg2rad(lat_histo[i])))
"""

for i in range(len(lat_histo)-1):
    for j in range(len(lon_all_histo)):
        for k in range(432):  
                mvimd_v0[k,i,j] = (96/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))*(umvimd_all_histo_v2[k,i,(j+1)%192]-umvimd_all_histo_v2[k,i,j])\
                    + (144/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))\
                    *(vmvimd_all_histo_v2[k,(i+1),j]*np.cos(deg2rad(lat_histo[i+1]))-vmvimd_all_histo_v2[k,i,j]*np.cos(deg2rad(lat_histo[i])))
                        
                mvimd_ssp_v0[k,i,j] = (96/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))*(umvimd_all_ssp_v2[k,i,(j+1)%192]-umvimd_all_ssp_v2[k,i,j])\
                    + (144/(np.pi*R*np.cos(deg2rad(lat_histo[i]))))\
                    *(vmvimd_all_ssp_v2[k,(i+1),j]*np.cos(deg2rad(lat_histo[i+1]))-vmvimd_all_ssp_v2[k,i,j]*np.cos(deg2rad(lat_histo[i])))


lon_all_histo=np.append(lon_all_histo, 360)

mvimd = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):    
            mvimd[i,j] = np.append(mvimd_v0[i,j],(mvimd_v0[i,j,0]+mvimd_v0[i,j,-1])/2)

mvimd_ssp = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mvimd_ssp[i,j] = np.append(mvimd_ssp_v0[i,j],(mvimd_ssp_v0[i,j,0]+mvimd_ssp_v0[i,j,-1])/2)

for i in range(432):
        mvimd[i,0,:] = 0
        mvimd_ssp[i,0,:] = 0

"""
mvimd_era = np.zeros((432,121,1440))
for i in range(len(lat_era)):
    for j in range(len(lon_era)):
        for k in range(432):  
            mvimd_era[k,i,j] = (720/(np.pi*R*np.cos(deg2rad(lat_era[i]))))*(umvimd_era[k,i,(j+1)%1440]-umvimd_era[k,i,j] + vmvimd_era[k,(i-1),j]*np.cos(deg2rad(lat_era[i-1]))-vmvimd_era[k,i,j]*np.cos(deg2rad(lat_era[i])))
"""        
        
    
#np.save("mvimd_era5", mvimd_era)

#=====================================================================================================

#faire par mois---------------------------------------------------------------------------------------------------------------------

def repartion_by_months (array, start, end):
    
    lenght = int((len(array))/12)
    jan = np.zeros((lenght,len(array[0]),len(array[0,0])))
    feb = np.zeros((lenght,len(array[0]),len(array[0,0])))
    mar = np.zeros((lenght,len(array[0]),len(array[0,0])))
    apr = np.zeros((lenght,len(array[0]),len(array[0,0])))
    may = np.zeros((lenght,len(array[0]),len(array[0,0])))
    jun = np.zeros((lenght,len(array[0]),len(array[0,0])))
    jul = np.zeros((lenght,len(array[0]),len(array[0,0])))
    aug = np.zeros((lenght,len(array[0]),len(array[0,0])))
    sep = np.zeros((lenght,len(array[0]),len(array[0,0])))
    octb = np.zeros((lenght,len(array[0]),len(array[0,0])))
    nov = np.zeros((lenght,len(array[0]),len(array[0,0])))
    dec = np.zeros((lenght,len(array[0]),len(array[0,0])))
  
    for i in range(start, end):
        if (i+1)%12 == 1:
            idx = (i//12)
            jan[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 2:
            idx = (i//12)
            feb[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 3:
            idx = (i//12)
            mar[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 4:
            idx = (i//12)
            apr[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 5:
            idx = (i//12)
            may[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 6:
            idx = (i//12)
            jun[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 7:
            idx = (i//12)
            jul[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 8:
            idx = (i//12)
            aug[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 9:
            idx = (i//12)
            sep[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 10:
            idx = (i//12)
            octb[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 11:
            idx = (i//12)
            nov[idx,:,:] =  (array[i,:,:])
        if (i+1)%12 == 0:
            idx = (i//12)
            dec[idx,:,:] =  (array[i,:,:])
            
        
    return [jan,feb,mar,apr,may,jun,jul,aug,sep,octb,nov,dec] 


def saison(array):
    DJF = np.zeros((len(array[0]),len(array[0][0]),len(array[0][0,0])))
    MAM = np.zeros((len(array[0]),len(array[0][0]),len(array[0][0,0])))
    JJA = np.zeros((len(array[0]),len(array[0][0]),len(array[0][0,0])))
    SON = np.zeros((len(array[0]),len(array[0][0]),len(array[0][0,0])))
    
    for i in range(len(array[0])):
        DJF[i,:,:] = (array[11][i,:,:] + array[0][i,:,:] + array[1][i,:,:])/3
        MAM[i,:,:] = (array[2][i,:,:] + array[3][i,:,:] + array[4][i,:,:])/3
        JJA[i,:,:] = (array[5][i,:,:] + array[6][i,:,:] + array[7][i,:,:])/3
        SON[i,:,:] = (array[8][i,:,:] + array[9][i,:,:] + array[10][i,:,:])/3
    
    return [DJF, MAM, JJA, SON]



mvimd_by_season_future = saison(repartion_by_months(mvimd, 0,432))
#mvimd_by_saison = saison(repartion_by_months(mvimd_era, 0,432))
mvimd_by_season_future_ssp = saison(repartion_by_months(mvimd_ssp, 0,432))


mean_mvimd_by_season_future = np.zeros((4,25,193))
mean_mvimd_by_season_future_ssp = np.zeros((4,25,193))
#mean_mvimd_era_by_season_future = np.zeros((4,121,1440))
for i in range(4):
    for j in range(36):
        mean_mvimd_by_season_future[i,:,:] += mvimd_by_season_future[i][j,:,:]/36
        mean_mvimd_by_season_future_ssp[i,:,:] += mvimd_by_season_future_ssp[i][j,:,:]/36
        #mean_mvimd_era_by_season_future[i,:,:] += mvimd_by_saison[i][j,:,:]/36
        
       

#figure--------------------------------------------------------------------------------------
#Plots =======================================================================
def plot_map(lon,lat, lb, ub, step, units ,args, name, which, arrow, my_color, sic_ = None):
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
    
    
    # ----------------
    # Plot the figure |
    # ----------------
    #(args, name, number) :
    #for jt in range(1):
    #jt = number                                  #la fin y a pas car on a pas les données jusqu'a fin 2021 :) 
    hemisphere = hemi
    #fname = filein.split("/")[-1] + "_" + varin + "_" + hemisphere + "_t"  
    #title = filein1.split("/")[-1] 
    
    
    field = np.squeeze(args[:, :])
    #sic = np.squeeze(sic_[:,:])
    
    seasons = ["summer", "autumn", "winter", "spring"]
    
    
    # Open figure
    #fig=plt.figure("fig", figsize = (6, 6), dpi = 500)
    fig = plt.figure(figsize = (12,12) , dpi = 150)
    for i in range(4):
        
        ax = plt.subplot(221+i)
        title = "{0} {1}".format(name, seasons[i]) 
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
            #si = map.contour(x,y, sic[i], [15,100], colors = "black" , linestyles = "dashed")
            cs = map.contourf(x, y, field[i], clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
            #ds = map.contourf(x, y, array_pval[i], [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
            cs.cmap.set_under(col_under)
            cs.cmap.set_over(col_over)
        elif which == "simple":
            # Create a contourf object called "cs"
            norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
            
            #si = map.contour(x,y, sic[i], [15,100], colors = "black" , linestyles = "dashed")
            cs = map.contourf(x, y, field[i], clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
            #ds = map.contourf(x, y, array_pval[i], [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
            cs.cmap.set_under(col_under)
            cs.cmap.set_over(col_over)
        
        
        #pt = map.scatter(my_x, my_y, marker = 'o', color='g', s = 10)
        #map.scatter(200, 100, marker = 'o', color='g', zorder=5, s = 10)
        #map.scatter(1300, 90, marker = 'o', color='g', zorder=5, s = 10)
        
        #Arrow======
        
        #u_norm = u / np.sqrt(u ** 2.0 + v ** 2.0)
        #v_norm = v / np.sqrt(u ** 2.0 + v ** 2.0)
        #map.quiver(x[::6,::20], y[::6,::20], u_norm[i,::6,::20], v_norm[i,::6,::20])
        
        
    
        plt.title(title,loc='center', pad=20)
        #plt.text(cs, 2, 0.5, '(b)', fontsize=13, fontweight = 'bold')
        #plt.annotate('(a)', xy=(0.03, 0.92), xycoords='axes fraction',fontsize=13, fontweight = 'bold')
        cax = fig.add_axes([0.98, 0.14, 0.04, 0.7])
        cbar = fig.colorbar(cs, cax=cax)
        cbar.ax.tick_params(labelsize=12)
        cbar.set_label(units,fontsize=13, fontweight = 'bold')

        
    plt.show()




#plot_map(lon,lat, 0, 55, 5 ,'$\mathregular{mm/month}$', mean_mer_by_season , "Evaporation (1979-2014) Austral", "double", "both" , "Reds", sic_ = mean_sic_by_season)   
#plot_map(lon,lat, 0, 110, 10 ,'$\mathregular{mm/month}$', mean_mtpr_by_season, "Precipitation (1979-2014) Autral", "simple", "max" , "Blues")   


#plot_map(lon_future,lat_future, 0, 55, 5 ,'$\mathregular{mm/month}$', mean_mer_by_season_future , "Evaporation (1979-2014) historical Austral", "double", "both" , "Reds", sic_ = mean_sic_by_season_future)   
#plot_map(lon_future,lat_future, 0, 110, 10 ,'$\mathregular{mm/month}$', mean_mtpr_by_season_future, "Precipitation (1979-2014) historical Autral", "simple", "max" , "Blues")   
plot_map(lon_all_histo ,lat_histo[0:25] , -20, 100, 10 ,'$\mathregular{mm/month}$', (-1)*mean_mvimd_by_season_future , "Moisture convergence (1979-2014) historical Austral", "double", "both" , "BrBG")   
plot_map(lon_all_histo ,lat_histo[0:25] , -20, 100, 10 ,'$\mathregular{mm/month}$', (-1)*mean_mvimd_by_season_future_ssp , "Moisture convergence (1979-2014) ssp585 Austral", "double", "both" , "BrBG")   


#plot_map(lon_era ,lat_era , -20, 100, 10 ,'$\mathregular{mm/month}$', mean_mvimd_era_by_season_future*(-1) , "Moisture convergence (1979-2014) historical Austral", "double", "both" , "BrBG")   



#np.save("mvimd_ACCESS_histo_corrected", mvimd)
#np.save("mvimd_ACCESS_ssp_corrected", mvimd_ssp)

