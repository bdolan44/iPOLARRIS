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
                                              rr=None,band = 'C',vr = 'vr',lat_r=lat_r,lon_r=lon_r,
                                              y=config['yname'],z=config['zname'],lat=config['xname'], lon=config['yname'],lat_0=lat_0,lon_0=lon_0,
                                              exper=config['exper'],mphys=config['mphys'],radar_name =config['radarname'],
                                              z_thresh=0,conv_types =  config['conv_types'],
                                               strat_types = config['strat_types'])
                                               
                                               
rdata.calc_cs_shy()

wup50 = rdata[rdata.w_name].where(rdata[rdata.w_name]>0).quantile(.5,dim=['x','y','d'])
wup90 = rdata[rdata.w_name].where(rdata[rdata.w_name]>0).quantile(.9,dim=['x','y','d'])
wup99 = rdata[rdata.w_name].where(rdata[rdata.w_name]>0).quantile(.99,dim=['x','y','d'])


wdn50 = rdata[rdata.w_name].where(rdata[rdata.w_name]<0).quantile(.5,dim=['x','y','d'])
wdn90 = rdata[rdata.w_name].where(rdata[rdata.w_name]<0).quantile(.1,dim=['x','y','d'])
wdn99 = rdata[rdata.w_name].where(rdata[rdata.w_name]<0).quantile(.01,dim=['x','y','d'])                                               

plt.plot(wup50,rdata[rdata.z_name],color='goldenrod',label='50th')
plt.plot(wup90,rdata[rdata.z_name],color='k',label='90th')
plt.plot(wup99,rdata[rdata.z_name],color='r',label='99th')
plt.legend(loc='best')
plt.plot(wdn50,rdata[rdata.z_name],color='goldenrod')
plt.plot(wdn90,rdata[rdata.z_name],color='k')
plt.plot(wdn99,rdata[rdata.z_name],color='r')

plt.title("TWP-ICE")

plt.xlabel('Vertical Velocity (m/s)')
plt.ylabel('Height (km)')
plt.xlim(-10,15)
plt.savefig('wup_down_wrf.png',dpi=200)