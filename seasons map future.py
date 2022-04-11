# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 11:02:50 2022

@author: noepi
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 21:04:46 2021

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


access_pr_file = 'pr_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
lat_future=np.array(dset['lat'][0:25])
lon_future=np.array(dset['lon'][:])
lon_future=np.append(lon_future, 360)
time=np.array(dset['time'][:])
mtpr_future = np.array(dset['pr'][1548:1980,0:25,:])*sec_month                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
mtpr_future_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mtpr_future_v2[i,j] = np.append(mtpr_future[i,j],(mtpr_future[i,j,0]+mtpr_future[i,j,-1])/2)
dset.close()                                                                   #2592000 si on est en [kg/((m^2)*month)]

access_pr_file = 'evspsbl_Amon_ACCESS-ESM1-5_historical_r1i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
mer_future = np.array(dset['evspsbl'][1548:1980,0:25,:])*sec_month                      #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
mer_future_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mer_future_v2[i,j] = np.append(mer_future[i,j],(mer_future[i,j,0]+mer_future[i,j,-1])/2)
dset.close()

access_pr_file = 'siconca_SImon_ACCESS-ESM1-5_historical_r11i1p1f1_gn_185001-201412.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
#lat_sic=np.array(dset['j'][0:25])
#lon_sic=np.array(dset['i'][:])
#lon_sic=np.append(lon_sic, 360)
sic_future = np.array(dset['siconca'][1548:1980,0:25,:])  
sic_future_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            sic_future_v2[i,j] = np.append(sic_future[i,j],(sic_future[i,j,0]+sic_future[i,j,-1])/2)      
dset.close()

access_pr_file = 'intuaw_Emon_ACCESS-ESM1-5_historical_r11i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
mvimd_future_u = np.array(dset['intuaw'][1548:1980,0:25,:])*(-1)               #*sec_month*(-1)       #*sec_year
mvimd_future_u_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mvimd_future_u_v2[i,j] = np.append(mvimd_future_u[i,j],(mvimd_future_u[i,j,0]+mvimd_future_u[i,j,-1])/2)
dset.close()

access_pr_file = 'intvaw_Emon_ACCESS-ESM1-5_historical_r11i1p1f1_gn_185001-201412.nc'
dset = nc.Dataset(access_pr_file)
mvimd_future_v = np.array(dset['intvaw'][1548:1980,0:25,:])*(-1)               #*sec_month*(-1)       #*sec_year
mvimd_future_v_v2 = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
            mvimd_future_v_v2[i,j] = np.append(mvimd_future_v[i,j],(mvimd_future_v[i,j,0]+mvimd_future_v[i,j,-1])/2)
dset.close()

mvimd_future = np.zeros((432,25,193))
for i in range(432):
    for j in range(25):
        for k in range(193):
            mvimd_future[i,j,k] = mvimd_future_u_v2[i,j,k] + mvimd_future_v_v2[i,j,k]



access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:432,0,:,:])*sec_month
dset.close()


access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:432,0,:,:])*(-sec_month)          #2592000 nécessaire car on a [kg/((m^2)*s)] et qu'on veut représenter sur notre map en [kg/((m^2)*month)]
dset.close()

access_pr_file = 'Sea_ice_cover.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
sic = np.array(dset['siconc'][0:432,0,:,:])        
dset.close()

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:432,0,:,:])*sec_month*(-1)                #*sec_month*(-1)       #*sec_year
dset.close()

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


mer_by_season_future = saison(repartion_by_months(mer_future_v2, 0, 432))
mtpr_by_season_future = saison(repartion_by_months(mtpr_future_v2, 0,432))
sic_by_season_future = saison(repartion_by_months(sic_future_v2, 0,432))
mvimd_by_season_future = saison(repartion_by_months(mvimd_future, 0,432))

mer_by_season = saison(repartion_by_months(mer, 0, 432))
mtpr_by_season = saison(repartion_by_months(mtpr, 0,432))
sic_by_season = saison(repartion_by_months(sic, 0, 432))
mvimd_by_saison = saison(repartion_by_months(mvimd, 0,432))


mean_mer_by_season_future = np.zeros((4,25,193))
mean_mtpr_by_season_future = np.zeros((4,25,193))
mean_mer_by_season = np.zeros((4,121,1440))
mean_mtpr_by_season = np.zeros((4,121,1440))
mean_sic_by_season = np.zeros((4,121,1440))
mean_sic_by_season_future = np.zeros((4,25,193))
mean_mvimd_by_season_future = np.zeros((4,25,193))
for i in range(4):
    for j in range(36):
        mean_mvimd_by_season_future[i,:,:] += mvimd_by_season_future[i][j,:,:]/36
        mean_sic_by_season_future[i,:,:] += sic_by_season_future[i][j,:,:]/36
        mean_sic_by_season[i,:,:] += sic_by_season[i][j,:,:]/36
        mean_mer_by_season_future[i,:,:] += mer_by_season_future[i][j,:,:]/36
        mean_mtpr_by_season_future[i,:,:] += mtpr_by_season_future[i][j,:,:]/36
        mean_mer_by_season[i,:,:] += mer_by_season[i][j,:,:]/36
        mean_mtpr_by_season[i,:,:] += mtpr_by_season[i][j,:,:]/36
        
       

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
    sic = np.squeeze(sic_[:,:])
    
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
            norm = colors.DivergingNorm(vmin=lb, vcenter= 5 , vmax=ub)#!!!
            si = map.contour(x,y, sic[i], [15,100], colors = "black" , linestyles = "dashed")
            cs = map.contourf(x, y, field[i], clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
            #ds = map.contourf(x, y, array_pval[i], [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
            cs.cmap.set_under(col_under)
            cs.cmap.set_over(col_over)
        elif which == "simple":
            # Create a contourf object called "cs"
            norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
            
            si = map.contour(x,y, sic[i], [15,100], colors = "black" , linestyles = "dashed")
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


plot_map(lon_future,lat_future, 0, 55, 5 ,'$\mathregular{mm/month}$', mean_mer_by_season_future , "Evaporation (1979-2014) historical Austral", "double", "both" , "Reds", sic_ = mean_sic_by_season_future)   
#plot_map(lon_future,lat_future, 0, 110, 10 ,'$\mathregular{mm/month}$', mean_mtpr_by_season_future, "Precipitation (1979-2014) historical Autral", "simple", "max" , "Blues")   
#plot_map(lon_future,lat_future, -20, 100, 10 ,'$\mathregular{mm/month}$', mean_mvimd_by_season_future , "Moisture convergence (1979-2014) historical Austral", "double", "both" , "BrBG")   


    