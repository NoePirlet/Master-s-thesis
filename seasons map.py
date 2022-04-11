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


access_pr_file = 'Mean_total_precipitation_rate.nc'
dset = nc.Dataset(access_pr_file)
lat=np.array(dset['latitude'][:])
lon=np.array(dset['longitude'][:])
time=np.array(dset['time'][:])
mtpr = np.array(dset['mtpr'][0:512,0,:,:])*sec_year          #*sec_month        #*sec_year                       #365*24*3600 si on est en [kg/((m^2)*year)]!!!!!!
dset.close()                                                                   #2592000 si on est en [kg/((m^2)*month)]


access_pr_file = '2m_temperature.nc'      
dset = nc.Dataset(access_pr_file)
sat = np.array(dset['t2m'][0:512,:,:])-273.15  #[C]
dset.close()

access_pr_file = 'Mean_vertically_integrated_moisture_divergence.nc'
dset = nc.Dataset(access_pr_file)
mvimd = np.array(dset['mvimd'][0:512,0,:,:])*sec_month*(-1)                #*sec_month*(-1)       #*sec_year
dset.close()

access_pr_file = 'Mean_evaporation_rate.nc'
dset = nc.Dataset(access_pr_file)
mer = np.array(dset['mer'][0:512,0,:,:])*sec_month*(-1)          #*sec_month*(-1)           #*sec_year
dset.close()

access_pr_file = 'Mean_sea_level_pressure.nc'       #[P]
dset = nc.Dataset(access_pr_file)
mslp = np.array(dset['msl'][0:512,0,:,:])/100             #/100       #nothing
dset.close()

access_pr_file = 'SST.nc'                           #[K]
dset = nc.Dataset(access_pr_file)
sst = np.array(dset['sst'][0:512,0,:,:])
dset.close()

access_pr_file = 'Sea_ice_cover.nc'                 #[No dim]
dset = nc.Dataset(access_pr_file)
sic = np.array(dset['siconc'][0:512,0,:,:])        
dset.close()


access_pr_file = '100m u-component of wind era5.nc'                 
dset = nc.Dataset(access_pr_file)
u100 = np.array(dset['u100'][0:504,0,:,:])
dset.close()

access_pr_file = '100m v-component of wind era5.nc'                 
dset = nc.Dataset(access_pr_file)
v100 = np.array(dset['v100'][0:504,0,:,:])
dset.close() 

access_pr_file = '10m u-component of wind.nc'                 
dset = nc.Dataset(access_pr_file)
u10 = np.array(dset['u10'][0:504,0,:,:])
dset.close() 

access_pr_file = '10m v-component of wind.nc'                 
dset = nc.Dataset(access_pr_file)
v10 = np.array(dset['v10'][0:504,0,:,:])
dset.close() 
   

#access_pr_file = 'Total_column_water.nc'
#dset = nc.Dataset(access_pr_file)
#tcw = np.array(dset['tcw'][0:512,0,:,:])
#dset.close()

#negative values indicate evaporation and positive values indicate deposition.
#access_pr_file = 'Mean_snow_evaporation_rate.nc'    
#dset = nc.Dataset(access_pr_file)
#mser = np.array(dset['mser'][0:512,0,:,:])*sec_year      
#dset.close()


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
short_pressure = year_array(mslp, year_to_months)
short_sic = year_array(sic, year_to_months)  

#faire par mois---------------------------------------------------------------------------------------------------------------------

def repartion_by_months (array, start, end):
    
    lenght = int((len(array))/12)
    jan = np.zeros((lenght,121,1440))
    feb = np.zeros((lenght,121,1440))
    mar = np.zeros((lenght,121,1440))
    apr = np.zeros((lenght,121,1440))
    may = np.zeros((lenght,121,1440))
    jun = np.zeros((lenght,121,1440))
    jul = np.zeros((lenght,121,1440))
    aug = np.zeros((lenght,121,1440))
    sep = np.zeros((lenght,121,1440))
    octb = np.zeros((lenght,121,1440))
    nov = np.zeros((lenght,121,1440))
    dec = np.zeros((lenght,121,1440))
  
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
    DJF = np.zeros((len(array[0]),121,1440))
    MAM = np.zeros((len(array[0]),121,1440))
    JJA = np.zeros((len(array[0]),121,1440))
    SON = np.zeros((len(array[0]),121,1440))
    
    for i in range(len(array[0])):
        DJF[i,:,:] = (array[11][i,:,:] + array[0][i,:,:] + array[1][i,:,:])/3
        MAM[i,:,:] = (array[2][i,:,:] + array[3][i,:,:] + array[4][i,:,:])/3
        JJA[i,:,:] = (array[5][i,:,:] + array[6][i,:,:] + array[7][i,:,:])/3
        SON[i,:,:] = (array[8][i,:,:] + array[9][i,:,:] + array[10][i,:,:])/3
    
    return [DJF, MAM, JJA, SON]



mer_by_season = saison(repartion_by_months(mer, 0 ,504)) 
mtpr_by_season = saison(repartion_by_months(mtpr, 0 ,504)) 
mvimd_by_season = saison(repartion_by_months(mvimd, 0 ,504))
PE = mtpr -(-mer)
PE_by_season = saison(repartion_by_months(PE, 0, 504)) 
pressure_by_season = saison(repartion_by_months(mslp, 0 ,504)) 
sic_by_season = saison(repartion_by_months(sic, 0 ,504))
sst_by_season = saison(repartion_by_months(sst, 0, 504))
u100_by_season = saison(repartion_by_months(u100, 0, 504))
v100_by_season = saison(repartion_by_months(v100, 0, 504))
u10_by_season = saison(repartion_by_months(u10, 0, 504))
v10_by_season = saison(repartion_by_months(v10, 0, 504))
sat_by_season  = saison(repartion_by_months(sat, 0, 504))


mean_sic_by_season = np.zeros((4,121,1440))
mean_mer_by_season = np.zeros((4,121,1440))
mean_mtpr_by_season = np.zeros((4,121,1440))
mean_mvimd_by_season = np.zeros((4,121,1440))
mean_pressure_by_season = np.zeros((4,121,1440))
mean_sst_by_season = np.zeros((4,121,1440))
mean_u100_by_season = np.zeros((4,121,1440))
mean_v100_by_season = np.zeros((4,121,1440))
mean_u10_by_season = np.zeros((4,121,1440))
mean_v10_by_season = np.zeros((4,121,1440))
#mean_PE_by_season = np.zeros((4,121,1440))
mean_sat_by_season = np.zeros((4,121,1440))
for i in range(4):
    for j in range(42):
        mean_sat_by_season[i,:,:] += (sat_by_season[i][j,:,:])/42
        mean_mvimd_by_season[i,:,:] += mvimd_by_season[i][j,:,:]/42
        mean_pressure_by_season[i,:,:] += pressure_by_season[i][j,:,:]/42
        mean_sic_by_season[i,:,:] += sic_by_season[i][j,:,:]/42
        mean_mer_by_season[i,:,:] += mer_by_season[i][j,:,:]/42
        mean_u10_by_season[i,:,:] += u10_by_season[i][j,:,:]/42
        mean_v10_by_season[i,:,:] += v10_by_season[i][j,:,:]/42
        mean_u100_by_season[i,:,:] += u100_by_season[i][j,:,:]/42
        mean_v100_by_season[i,:,:] += v100_by_season[i][j,:,:]/42
        
        mean_sst_by_season[i,:,:] += sst_by_season[i][j,:,:]/42
        
        mean_mtpr_by_season[i,:,:] += mtpr_by_season[i][j,:,:]/42
        
        
        #mean_PE_by_season[i,:,:] = PE_by_season[i][j,:,:]
        

mer_by_season_1997 = saison(repartion_by_months(short_mer, 0 ,12))
mtpr_by_season_1997 = saison(repartion_by_months(short_mtpr, 0 ,12))
mvimd_by_season_1997 = saison(repartion_by_months(short_mvimd, 0 ,12))
pressure_by_season_1997 = saison(repartion_by_months(short_pressure, 0 ,12))   
sic_by_season_1997 = saison(repartion_by_months(short_sic, 0 ,12))

#pour mettre sous la forme désiré
mer_by_season_1997_3d = np.zeros((4,121,1440))
mtpr_by_season_1997_3d = np.zeros((4,121,1440))
mvimd_by_season_1997_3d = np.zeros((4,121,1440))
pressure_by_season_1997_3d = np.zeros((4,121,1440))
sic_by_season_1997_3d = np.zeros((4,121,1440))
sic_by_season_3d = np.zeros((4,121,1440))
for i in range(4):
    mer_by_season_1997_3d[i,:,:] = mer_by_season_1997[i]
    mtpr_by_season_1997_3d[i,:,:] = mtpr_by_season_1997[i]
    mvimd_by_season_1997_3d[i,:,:] = mvimd_by_season_1997[i]
    pressure_by_season_1997_3d[i,:,:] = pressure_by_season_1997[i]   
    sic_by_season_1997_3d[i,:,:] = sic_by_season_1997[i]
    


   
#Regression linéaire ------------------------------------------------------------------------

def years(nbr):
    years = np.zeros(nbr)
    for i in range(nbr):
        years[i] = 1979 + i
    return years


def lin_reg (x, y):
    linreg = linregress(x, y)
    m = linreg.slope
    b = linreg.intercept
    p = linreg.pvalue
    return(m,b,p)


def linreg_grid(years, array):
    linreg_g =  np.full((121, 1440), None)
    p_val =  np.full((121, 1440), None)
    for i in range(121):
        for j in range(1440):
            lin_val = lin_reg(years, array[:,i,j])
            linreg_g[i,j] = lin_val[0]
            p_val[i,j] = lin_val[2]
    return linreg_g,p_val
    

lin_mer_season = np.zeros((4,121,1440))
mer_pval = np.zeros((4,121,1440))
lin_mtpr_season = np.zeros((4,121,1440))
mtpr_pval = np.zeros((4,121,1440))
lin_mvimd_season = np.zeros((4,121,1440))
mvimd_pval = np.zeros((4,121,1440))
lin_pressure_season = np.zeros((4,121,1440))
pressure_pval = np.zeros((4,121,1440))
lin_PE_season = np.zeros((4,121,1440))
PE_pval = np.zeros((4,121,1440))

for i in range(4):
    lin_mer_season[i], mer_pval[i] = linreg_grid(years(42), -mer_by_season[i])
    lin_mtpr_season[i], mtpr_pval[i] = linreg_grid(years(42), mtpr_by_season[i])
    lin_mvimd_season[i], mvimd_pval[i] = linreg_grid(years(42), -mvimd_by_season[i])
    lin_pressure_season[i], pressure_pval[i] = linreg_grid(years(42), pressure_by_season[i])
    lin_PE_season[i], PE_pval[i] = linreg_grid(years(42), PE_by_season[i])

#figure--------------------------------------------------------------------------------------
#Plots =======================================================================
def plot_map(lon,lat, lb, ub, step, units ,args, name, which, arrow, my_color, array_pval = None, sic_= None, u =None ,v =None):
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
    
    #zm = np.ma.masked_less(p_val, 0.05)
    
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
    
    #my_lon = [-64,145,145]
    #my_lat = [-67,-68,-64]
    
    #création d'une zone d'étude 
    my_lat = []
    my_lon = []
    
    #ABS box
    for j in range(20):
        for i in range(28):
            latitude = -65 + -0.25*j
            my_lat.append(latitude)
        
            longitude = -65 - 0.25*i
            my_lon.append(longitude)  
    """
    #PO Sea box        
    for j in range(16):
        for i in range(64):
            latitude = -61 + -0.25*j
            my_lat.append(latitude)
        
            longitude = 127 + 0.25*i
            my_lon.append(longitude)
    """
    #PO Coastal box        
    for j in range(16):
        for i in range(64):
            latitude = -63 + -0.25*j
            my_lat.append(latitude)
        
            longitude = 127 + 0.25*i
            my_lon.append(longitude)
       
        
        
    
    #my_lon = [-66,135,140]
    #my_lat = [-67,-65,-63]
    my_x, my_y = map(my_lon, my_lat)
    
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
            #si = map.contour(x,y, sic[i], [0.15,1], colors = "black" , linestyles = "dashed")
            cs = map.contourf(x, y, field[i], clevs, cmap = cmap, vmin=lb, vmax=ub, norm=norm, latlon = False, extend = extmethod)
            #ds = map.contourf(x, y, array_pval[i], [0,0.05], colors = "black", alpha = 0.1,  hatches=['///'])
            cs.cmap.set_under(col_under)
            cs.cmap.set_over(col_over)
        elif which == "simple":
            # Create a contourf object called "cs"
            norm = colors.Normalize(vmin=lb, vmax=ub)#!!!
            
            #si = map.contour(x,y, sic[i], [0.15,1], colors = "black" , linestyles = "dashed")
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

        #ax.plot(my_x, my_y,'bo', markersize = 1, alpha = 0.1)
        
        
        
    plt.show()




#plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year}$', lin_mer_season , "Trends of evaporation (1979-2020) Austral", "double", "both" , "coolwarm", mer_pval, mean_sic_by_season)   
#plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year}$', lin_mtpr_season , "Trends of precipitation (1979-2020) Austral", "double", "both" , "coolwarm_r", mtpr_pval)   
#plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{mm/year}$', lin_mvimd_season , "Trends of Moisture convergence (1979-2020) Austral", "double", "both" , "coolwarm", mvimd_pval)   
#plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{Pa/year}$', lin_pressure_season , "Trends of Pressure (1979-2020) Austral", "double", "both" , "coolwarm_r", pressure_pval)   
#plot_map(lon,lat, -1.2, 1.4, 0.2 ,'$\mathregular{Pa/year}$', lin_PE_season , "Trends of PE (1979-2020) Austral", "double", "both" , "coolwarm_r", PE_pval)   


#plot_map(lon,lat, 0, 55, 5 ,'$\mathregular{mm/month}$', mer_by_season_1997_3d , "Evaporation (1997) Austral", "simple", "max" , "Reds", mer_pval, sic_by_season_1997_3d)   
#plot_map(lon,lat, 0, 110, 10 ,'$\mathregular{mm/month}$', mtpr_by_season_1997_3d , "Precipitation (1997) Autral", "simple", "max" , "Blues", mtpr_pval)   
#plot_map(lon,lat, -20, 100, 10 ,'$\mathregular{mm/month}$', mvimd_by_season_1997_3d , "Moisture convergence (1997) Austral", "double", "both" , "BrBG", mvimd_pval)   
#plot_map(lon,lat, 980, 1018, 2 ,'$\mathregular{hPa}$', pressure_by_season_1997_3d , "Pressure (1997) Austral", "simple", "both" , "coolwarm_r", pressure_pval)   


#plot_map(lon,lat, 0, 55, 5 ,'$\mathregular{mm/month}$', mean_mer_by_season , "Evaporation (1979-2020) Austral", "double", "both" , "Reds",sic_ = mean_sic_by_season)   
#plot_map(lon,lat, 0, 110, 10 ,'$\mathregular{mm/month}$', mean_mtpr_by_season, "Precipitation (1979-2020) Autral", "simple", "max" , "Blues")   
#plot_map(lon,lat, -20, 100, 10 ,'$\mathregular{mm/month}$', mean_mvimd_by_season , "Moisture convergence (1979-2020) Austral", "double", "both" , "BrBG")   
#plot_map(lon,lat, 980, 1018, 2 ,'$\mathregular{hPa}$', mean_pressure_by_season , "Pressure (1979-2020) Austral", "simple", "both" , "coolwarm_r")   
#plot_map(lon,lat, 271, 275.5, 0.5 ,'$\mathregular{K}$', mean_sst_by_season , "SST (1979-2020) Austral", "simple", "both" , "coolwarm_r", sic_ = mean_sic_by_season)   
  
plot_map(lon,lat, -30, 9, 3 ,'$\mathregular{Celsius}$', mean_sat_by_season , "SAT (1979-2020) Austral", "double", "both" , "coolwarm")   
         
"""
!!!!!!!!!!!!! pour la map de l'évaporatation plot la sic
              pour la map de la précipitation plot la topographie 
"""
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
