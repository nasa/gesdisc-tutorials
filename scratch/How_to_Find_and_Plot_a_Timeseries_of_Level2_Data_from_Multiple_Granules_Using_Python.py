#!/usr/bin/env python
# coding: utf-8

# # How to Find and Plot a timeseries Level 2 Data from Multiple Granules using Python

# ## Overview
# 
# Level 2 data from Low Earth Orbit (LEO) satellites are often organized into \"granules\" that contain spatial sets of data that may not match up with the region of interest.  If the region is smaller than the granule, the user can create a spatial subset.  If the region spans mulitple granules, each granule may be plotted individually, or the data can be aggregated and plotted at once.
# 
# This How-To assumes that the user has searched for the granules from the [JPSS-1 CrIS Level 2 CLIMCAPS: Atmosphere cloud and surface geophysical state V2](https://doi.org/10.5067/LESQUBLWS18H) data over the **Gulf of Mexico** on **2021-09-23** within a bounding box defined by the following pairs of longitude and latitude: **(-100, 15)** and **(-75, 30)**.  This data search may be achieved either by using the GES DISC web interface or by performing a granule search as demostrated below using the [Common Metadata Repository (CMR)](https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html).
# 
# The data plotted in this example is the **Surface Air Temperature** from the CLIMCAPS product, but can be modified to apply to other other Level 2 products.
# 
# 
# --- 
# 
# *Dataset Citation:*
# 
# Chris Barnet (2019), Sounder SIPS: JPSS-1 CrIS Level 2 CLIMCAPS: Atmosphere cloud and surface geophysical state V2, Greenbelt, MD, USA, Goddard Earth Sciences Data and Information Services Center (GES DISC), Accessed: [Data Access Date], 10.5067/LESQUBLWS18H
# 
# ## Prerequisites
# 
# This example code is written in Python 3.9 and requires these libraries and files: 
# - numpy
# - requests
# - netCDF 
# - matplotlib
# - cartopy  
# - netrc file with valid Earthdata Login credentials ([How to Generate Earthdata Prerequisite Files](https://disc.gsfc.nasa.gov/information/howto?keywords=prerequisite&title=How%20to%20Generate%20Earthdata%20Prerequisite%20Files))
# - Approval to access the GES DISC archives with your Earthdata credentials (https://disc.gsfc.nasa.gov/earthdata-login)
# - wget (https://www.gnu.org/software/wget/)
# 
# *Note: If the files to be aggregated are not in netCDF format, a different library (e.g. Xarray) will have to be used.*

# ### 1. Import Libraries

# In[1]:


import glob
import numpy as np
import netCDF4
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import os
import requests
from geopy.distance import lonlat, distance
import tai
from cmr import CollectionQuery, GranuleQuery, ToolQuery, ServiceQuery, VariableQuery
import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
import csv
import pdb
import pandas as pd  
import myio
from myio import list2txt, txt2list
import h5py

# ### 2. Extract URLs from the CMR and Download Granules
# 
# To search for URLs of the granules of interest, we perform a CMR query using the ShortName (SNDRJ1IML2CCPRET) and Version ID of the dataset as well as the time range and bounding box of the time and region of interest. These URLs will be stored in a list for downloading.

# In[2]:


cmr_url = 'https://cmr.earthdata.nasa.gov/search/granules'


#Manaus
# There is Manaus TCCON data from 20140930_20150727
# There is ObsPack data from "2017-04-29T20:23:59Z" to "2020-07-12T17:30:54Z"
###title_string = "Manaus"
###pointlon = -60.209
###pointlat = -2.595
###pointradius = 100 # km
###time_start = datetime.datetime(2014,9,6)
###time_end = datetime.datetime(2023,7,1)

#OCO-2 data
ShortName = "OCO2_L2_Lite_FP"
VersionID = '11.1r'



#Mauna Loa
title_string = "Mauna Loa"
pointlon = -155.57
pointlat = 19.536
pointradius = 100 # km
time_start = datetime.datetime(2022,3,1)
time_end = datetime.datetime(2023,3,1)


'''
#Caltech
pointlon = -118.13
pointlat = 34.14
pointradius = 100 # km
'''

'''
#Plumas County
title_string = "Plumas County"
pointlon = -120.83
pointlat = 40.01
pointradius = 100 # km
time_start = datetime.datetime(2020,1,1)
time_end = datetime.datetime(2023,1,1)
'''


datacircle_string = "{:8.3f},{:8.3f},{:8f}".format(pointlon,pointlat,pointradius).replace(' ','')



# This is the aggregated file of data to be plotted
aggfile = title_string.replace(' ','')+'_'+ShortName+'.'+VersionID+'_'+time_start.strftime('%Y%m%d')+'to'+time_end.strftime('%Y%m%d')+'_Circle'+datacircle_string+'.h5'






# python-cmr can be used to create a list of URLs.

api = GranuleQuery()
gran_met = api.short_name(ShortName).version(VersionID).temporal(time_start.strftime('%Y-%m-%dT%H:%M:%SZ'),time_end.strftime('%Y-%m-%dT%H:%M:%SZ')).get(1000000) 




# First, get the list of all entry URL links returned by granule metadata query
all_links = [entry['links'] for entry in gran_met]

# Next, get just the granule URL links, exclude inherited links (i.e. from collection metadata).
gran_links = [link for group in all_links for link in group if 'inherited' not in link.keys()]

# Get just the URLs (identified by 'via S3' in the link title)
s3_urls = [link['href'] for link in gran_links if 'via S3' in link['title']]

# Get just the URLs (identified by 'OPENDAP location for the granule' in the link title)
opendap_urls = [link['href'] for link in gran_links if 'OPENDAP location for the granule' in link['title']]

# Finally, get just the URLs (identified by 'Download' in the link title)
data_urls = [link['href'] for link in gran_links if 'Download' in link['title']]






# Another alternative would be to use the GES DISC subsetter to only download the data in the region of interest as described in the following How-To: https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Use%20the%20Web%20Services%20API%20for%20Subsetting
    
# The following block of code uses wget to download the Level 2 files in the region of interest to a directory called the ShortName.VersionID of the product.

# In[ ]:


data_rootdir = './data/' 
datadir = data_rootdir +ShortName + '.' + VersionID # this is where I will download the files
if os.path.exists(datadir) is False:
    os.mkdir(datadir)


list2txt(data_urls,'dataurls.txt')


# change the data urls to local path for the local files.
localdatafiles = []
for  url in data_urls:
    localdatafiles.append(datadir + '/' + os.path.basename(url))


os.system('wget -nc -i dataurls.txt -P '+datadir)  # skip for now since I already did this
'oco2L2local.txt'

list2txt(localdatafiles,'oco2L2local.txt')




# ### 3. Aggregate and Plot

# The following block of code reads the names of the local netCDF files and aggregates the data in the region of interest.


# In[4]:

# coordinate variables
lon_var = 'longitude'
lat_var = 'latitude'
time_var = 'time'

coordinate_vars = [lon_var,lat_var,time_var]

# science variables
science_vars = ['xco2','xco2_quality_flag','Sounding/land_fraction','solar_zenith_angle']









if os.path.exists(aggfile) is True:
    # Don't create the file if it already exists
    print(aggfile+' already exists.') 
     
else:
    # I'll create it.        
    N=0

    # put all coordinate and science variables into a dictionary
    plot_dic = {}
    for var in coordinate_vars+science_vars:
        plot_dic[var] = []

    for ifile,file in enumerate(localdatafiles):
        print('Reading '+file)
        fid = netCDF4.Dataset(file ,mode='r',format='NETCDF4')
        longitude = fid.variables[lon_var][:]
        latitude = fid.variables[lat_var][:]

        for i in range(len(longitude)):
            dist = distance(lonlat(pointlon,pointlat),lonlat(longitude[i],latitude[i])).km
            if dist <= pointradius:
                print('I found one')
                N = N+1
                # Now I'll fill the plot_dic with the relevant coordinate and science data
                for var in plot_dic.keys():
                    plot_dic[var].append(fid[var][:].data[i])


    # now lets make it a dictionary of arrays rather than lists
    for var in plot_dic.keys():
        plot_dic[var] = np.array(plot_dic[var])

    # now write out all of the variables so I don't have to aggregate it again

    myio.dict2h5(plot_dic,aggfile,fill_value=-999999.0)


The previous stuff will be skipped because I already have an aggfile.    


    

# Now lets read the file that we just created

saved_dic = h5py.File(aggfile, "r")
plot_dic = {}
for key1 in saved_dic.keys():
    if type(saved_dic[key1]) == h5py._hl.dataset.Dataset:
        plot_dic[key1] = saved_dic[key1][:]
    else: # I'll assume it's a group
        for key2 in saved_dic[key1].keys():
            plot_dic[key1+'/'+key2] = saved_dic[key1][key2][:] # I'm not going any deeper than this for now.

    
    
# Now let's plot the variables on a map 

# this sets up the plot


ax = plt.axes(projection=ccrs.PlateCarree(),label='allmap')
ax.coastlines()
gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='black')
lowerlimit = np.min(plot_dic['xco2'])
upperlimit = np.max(plot_dic['xco2'])
ncolors=10
cmap=plt.cm.rainbow
bounds = np.linspace(lowerlimit,upperlimit,ncolors+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
c=plt.scatter(plot_dic['longitude'],plot_dic['latitude'],c=plot_dic['xco2'], cmap=cmap, norm=norm, edgecolor='black', linewidth=0)
plt.colorbar(ax=ax, orientation="vertical", pad=0.110, spacing='proportional', ticks=bounds, boundaries=bounds, format='%7.2f', extend = "both")
ax.set_title(title_string)
filename = title_string.replace(' ','')+'_'+time_start.strftime('%Y%m%d')+'to'+time_end.strftime('%Y%m%d')+'_allmap.png'
plt.savefig(filename)
ax.cla()

plt.show()

# Now apply the quality flags 

ocogood = plot_dic['xco2_quality_flag'] == 0
ocogoodland = np.logical_and(ocogood,plot_dic['Sounding/land_fraction'] > 0)
ocogoodwater = np.logical_and(ocogood,plot_dic['Sounding/land_fraction'] ==0)


# this sets up the plot
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='black')
lowerlimit = np.min(plot_dic['xco2'])
upperlimit = np.max(plot_dic['xco2'])
ncolors=10
cmap=plt.cm.rainbow
bounds = np.linspace(lowerlimit,upperlimit,ncolors+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
c=plt.scatter(plot_dic['longitude'][ocogood],plot_dic['latitude'][ocogood],c=plot_dic['xco2'][ocogood], cmap=cmap, norm=norm, edgecolor='black', linewidth=0)
plt.colorbar(ax=ax, orientation="vertical", pad=0.110, spacing='proportional', ticks=bounds, boundaries=bounds, format='%7.2f', extend = "both")
ax.set_title(title_string)





ax = plt.axes(projection=ccrs.PlateCarree(),label='goodmap')
ax.coastlines()
gl=ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=0.75, color='black')
lowerlimit = np.min(plot_dic['xco2'][ocogood])
upperlimit = np.max(plot_dic['xco2'][ocogood])
ncolors=10
cmap=plt.cm.rainbow
bounds = np.linspace(lowerlimit,upperlimit,ncolors+1)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
c=plt.scatter(plot_dic['longitude'][ocogood],plot_dic['latitude'][ocogood],c=plot_dic['xco2'][ocogood], cmap=cmap, norm=norm, edgecolor='black', linewidth=0)
plt.colorbar(ax=ax, orientation="vertical", pad=0.110, spacing='proportional', ticks=bounds, boundaries=bounds, format='%7.2f', extend = "both")
ax.set_title(title_string)
filename = title_string.replace(' ','')+'_'+time_start.strftime('%Y%m%d')+'to'+time_end.strftime('%Y%m%d')+'_goodmap.png'
plt.savefig(filename)
ax.cla()







# now plot a time series of the data over Mauna Loa

# we must first convert the TAI time to a datetime array
oco2l2_utc = tai.tai_to_datetime(plot_dic['time'][:],1970)


# OCO-2 data
fig, ax = plt.subplots()
ax.plot_date(oco2l2_utc[ocogood], plot_dic['xco2'][ocogood],label='OCO-2 L2',color='k',marker='o')

ax.set_title(title_string)

ax.set_ylim([np.floor(lowerlimit),np.ceil(upperlimit)])
ax.set_ylabel('xCO2 [ppm]')
ax.fmt_xdata = DateFormatter('%Y-%m-%d')
ax.grid(True)
fig.autofmt_xdate()
plt.legend()
filename = title_string.replace(' ','')+'_'+time_start.strftime('%Y%m%d')+'to'+time_end.strftime('%Y%m%d')+'_xco2_ts_oco2.png'
plt.savefig(filename)
ax.cla()



if title_string == "Mauna Loa":
    # Lets compare this plot with hourly in-situ measurements from Mauna Loa

    mlo_insitu_daily = pd.read_csv('./data/MaunaLoa/daily_in_situ_co2_mlo.csv',skiprows=34,names=['Yr','Mn','Dy','CO2','NB','scale']) 

    # Make a datetime array for the observations
    mlo_insitu_utc = []
    for i in range(len(mlo_insitu_daily['Dy'])):
        mlo_insitu_utc.append(datetime.datetime(mlo_insitu_daily['Yr'][i],mlo_insitu_daily['Mn'][i],mlo_insitu_daily['Dy'][i]))

    # convert the insitu data to an array
    mlo_insitu_utc = np.array(mlo_insitu_utc) 
    mlo_insitu_xco2 = list(mlo_insitu_daily['CO2'])

    # Now convert the NaN to a fill value that is valid.
    for i, xco2 in enumerate(mlo_insitu_xco2):
        if 'NaN' in xco2:
            mlo_insitu_xco2[i]=-999999.0
        else:
            mlo_insitu_xco2[i] = float(xco2)

    mlo_insitu_xco2 = np.array(mlo_insitu_xco2)

    mlo_valid = mlo_insitu_xco2 != -999999.
    # OCO-2 data with Daily In-situ
    fig, ax = plt.subplots()
    ax.plot_date(mlo_insitu_utc[mlo_valid], mlo_insitu_xco2[mlo_valid], label='In-situ',color='k',marker='x',ms=8)
    ax.plot_date(oco2l2_utc[ocogoodland], plot_dic['xco2'][ocogoodland],label='OCO-2 L2 Land',color='r',marker='o',ms=8)
    ax.plot_date(oco2l2_utc[ocogoodwater], plot_dic['xco2'][ocogoodwater],label='OCO-2 L2 Water',color='b',marker='o',ms=5)

    ax.set_title(title_string)
    ax.set_xlim([time_start,time_end])
    ax.set_ylim([np.floor(lowerlimit)-5,np.ceil(upperlimit)+5])
    ax.set_ylabel('xCO2 [ppm]')
    ax.fmt_xdata = DateFormatter('%Y-%m-%d')
    ax.grid(True)
    fig.autofmt_xdate()
    plt.legend()
    filename = title_string.replace(' ','')+'_'+time_start.strftime('%Y%m%d')+'to'+time_end.strftime('%Y%m%d')+'_xco2_ts_oco2withinsitu.png'
    plt.savefig(filename)

    
    
