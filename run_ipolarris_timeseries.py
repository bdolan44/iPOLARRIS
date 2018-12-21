import numpy as np
import os
import glob
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

import numpy as np

import RadarData
import datetime

import RadarConfig
import plot_driver
from polarris_config import run_exper
from polarris_config import get_data
import warnings
warnings.filterwarnings('ignore')
import GeneralFunctions as GF
from skewPy import SkewT
from collections import OrderedDict

import os
import sys


from matplotlib.dates import DateFormatter,HourLocator
dayFormatter = DateFormatter('%H%M')      # e.g., 12
hourFormatter = DateFormatter('%H')      # e.g., 12


def hasNumbers(inputString):
     return any(char.isdigit() for char in inputString)



configfile = sys.argv[1:]
config = {}
#print sys.argv[1:]

with open(configfile[0]) as f:
    for line in f:
        #print line
        if not line.startswith("#"):
            key, val, comment = line.split('==')
            vval = val.replace(" ","")
            numck = hasNumbers(vval)
            if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'mphys':
                numck = False
            if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'extra' or key.replace(" ", "") == 'ptype':
                vval = vval.strip("''")
            #print numck
            #print vval,key
            if key.replace(" ", "") == 'image_dir':
                numck = True

            if numck is True or vval == 'None' or vval == 'True' or vval == 'False':
                try:
                    config[(key.replace(" ", ""))] = eval(vval)
                except:
                    if "datetime" in vval:
                        config[(key.replace(" ", ""))] = vval
            else:
                config[(key.replace(" ", ""))] = vval
            

with open(config['radar_files'], 'r') as f:
    rfiles = f.read().splitlines()
#rfiles= glob.glob('*.nc')
rvar = xr.open_mfdataset(rfiles,concat_dim='d')#,preprocess=lambda ds:ds.drop(['time']),concat_dim='d')

lon_0 = 131.04444
lat_0 = -12.24917

lat_r = -12.24917
lon_r = -131.04444

tm = []
for d in rfiles:
    #print d
    dformat = config['wdate_format']
    base = os.path.basename(d)
    radcdate=np.str(base[config['time_parse'][0]:config['time_parse'][1]])
#    print radcdate
    date=datetime.datetime.strptime(radcdate,dformat)
    tm.append(date)



rdata = RadarData.RadarData(rvar,tm,ddata = None,dz =config['dz_name'],zdr=config['dr_name'],
                                              kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],
                                              u=config['uname'],v=config['vname'],w=config['wname'],x=config['xname'],
                                              rr=config['rr_name'],band = 'C',vr = 'vr',lat_r=lat_r,lon_r=lon_r,
                                              y=config['yname'],z=config['zname'],lat=config['xname'], lon=config['yname'],lat_0=lat_0,lon_0=lon_0,
                                              exper=config['exper'],mphys=config['mphys'],radar_name =config['radarname'],
                                              z_thresh=0,conv_types =  config['conv_types'],
                                               strat_types = config['strat_types'])
                                               
                                               
rdata.calc_cs_shy()
rdata.set_hid()

###Do some quick masking of the data####
mask = np.zeros([rdata.data.dims['d'],rdata.data.dims['z'],rdata.data.dims['y'],rdata.data.dims['x']])
whbad = np.logical_or(np.logical_or(np.logical_or(np.logical_or(rdata.data[rdata.dz_name].values>-20.,rdata.data[rdata.kdp_name].values>-10.),rdata.data[rdata.kdp_name].values<10.),rdata.data[rdata.zdr_name].values<10.),rdata.data[rdata.dz_name].values<70.)

mask[whbad] = 1
if np.nanmin(rdata.data['CSS'].values)<1.:
    mask[rdata.data['CSS'].values<=0] = 0
else:
    mask[np.isnan(rdata.data['CSS'].values)]=0

rdata.data['CSS'] = rdata.data['CSS'].where(mask ==1)
config['image_dir'] ='./'
#########################################

################################################################################
##################Now you can just start plotting!##############################
################################################################################


################################################################################
##First make a timeseries of rain rate, unconditional and conditional. This puts strat, conv, and total on the same plot but you can split the out by putting cs==False.
## The conditional rain rate is achieved by sending threshold = 0.
fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0)
ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=-50,ls='--',typ='uncond')

ax.set_ylabel('Rain Rate (mm/hr)')
ax.set_title('Precipitation Timeseries TWP-ICE')
plt.tight_layout()
plt.savefig('{i}Precip_timeseries_convstrat_{e}_{m}_{x}.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)
############################################################################

################################################################################
##Next let's make quantile (50,90,99) plots of the vertical velocity. This splits it by up and down, but you can turn split_updn == False
fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=True)
ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=False)
ax.set_xlabel('Vertical velocity m/s')
ax.set_title('Vertical velocity profiles TWP-ICE')
plt.tight_layout()
plt.savefig('{i}Quantile_vvel_{e}_{m}_{x}.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)
################################################################################

################################################################################
##Next let's make mean vertical profile of reflectivity
fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_verprof(rdata.data[rdata.dz_name],rdata.data[rdata.z_name],ax,split_updn=False,lab='dz',thresh=-50)
ax.set_title('Vertical profile of reflectivity')
ax.set_xlabel('Reflectivity')
plt.tight_layout()
plt.savefig('{i}MeanProfile_refl_{e}_{m}_{x}.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)


################################################################################
##Next let's make a reflectivity CFAD

cfaddat,vbins = plot_driver.cfad(rdata.data[rdata.dz_name],rdata,rdata.data[rdata.z_name].sel(d=0),nbins=40)

fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_cfad(cfaddat,rdata.data[rdata.z_name].sel(d=0).values,vbins,ax,levels=True)
ax.set_xlabel('Reflectivity')
ax.set_ylabel('Height (km)')
ax.set_title('TWP-ICE CFAD')
plt.tight_layout()
plt.savefig('{i}CFAD_refl_{e}_{m}_{x}.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)