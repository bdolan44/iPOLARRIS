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
from polarris_driver_new import polarris_driver
import os
import sys





configfile = sys.argv[1:]
#print sys.argv[1:]

rdata, config = polarris_driver(configfile)
#config['image_dir'] ='./'
#########################################

if sys.argv[2:]:
    configfile1 = sys.argv[2:]
    rdata2, config2 = polarris_driver(configfile1)

    print('calculating CFAD differences')

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=True,cscfad=False)
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,'Z$_{dr}$',config,config2,bins=np.arange(-2,8,0.2),savefig=True,cscfad=False)
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,'K$_{dp}$',config,config2,bins=np.arange(-2,6,0.2),savefig=True,cscfad=False)
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad=False)


    fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad=None)
    ##Convective

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=True,cscfad='convective')
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,'Z$_{dr}$',config,config2,bins=np.arange(-2,8,0.2),savefig=True,cscfad='convective')
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,'K$_{dp}$',config,config2,bins=np.arange(-2,6,0.2),savefig=True,cscfad='convective')
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad='convective')


    fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad='convective',savefig=True)

################################################################################
##################Now you can just start plotting!##############################
################################################################################

### To see the variables that are available to plot, type:

#rdata.data.keys()


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
plt.clf()
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
plt.clf()
################################################################################

################################################################################
##Next let's make mean vertical profile of reflectivity
fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_verprof(rdata.data[rdata.dz_name],rdata.data[rdata.z_name],ax,split_updn=False,lab='dz',thresh=-50)
ax.set_title('Vertical profile of reflectivity')
ax.set_xlabel('Reflectivity')
plt.tight_layout()
plt.savefig('{i}MeanProfile_refl_{e}_{m}_{x}.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)

plt.clf()
################################################################################
##Next let's make a reflectivity CFAD

cfaddat,vbins = plot_driver.cfad(rdata.data[rdata.dz_name],rdata,rdata.data[rdata.z_name],var=rdata.dz_name,nbins=40)

fig,ax = plt.subplots(1,1,figsize=(10,10))
ax = plot_driver.plot_cfad(cfaddat,rdata.data[rdata.z_name].values,vbins,ax,levels=True,cont=True)
ax.set_xlabel('Reflectivity')
ax.set_ylabel('Height (km)')
ax.set_title('TWP-ICE CFAD')
plt.tight_layout()
plt.savefig('{i}CFAD_refl_{e}_{m}_{x}_new.png'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extra']),dpi=400)
plt.clf()