import numpy as n
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
from polarris_config import get_data
import warnings
warnings.filterwarnings('ignore')
import GeneralFunctions as GF


def foo(s1):
    return '{}'.format(s1.rstrip())
#######################################################
#######################################################
def get_time(time_parse,filename,dformat):
    tp = time_parse[0]
    te = time_parse[1]
#    print filename
    base = os.path.basename(filename)
    radar_name = base[:3]
    radcdate=np.str(base[tp:tp+te])
    date=datetime.datetime.strptime(radcdate,dformat)
    
    return date

radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_mc3e.txt'
image_dir = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/'
yp = 'wrf'
exper = 'MC3E'
mphys = 'SBM'
date = '20110523'
time_parse=[11,19]
wdate_format ='%Y-%m-%d_%H:%M:%S'

ptype = 'png'


dum =[]
with open(radar_files) as f: 
#    dum.append(foo(f.readline()))
#    dum.append(foo(f.readline()))

    for line in f:
        dat = (line)
        dum.append(foo(dat))
dzcfad_a=[]
drcfad_a=[]
kdcfad_a=[]
rhcfad_a=[]
hidcfad_a=[]

hts_a=[]
dzbins_a=[]
drbins_a=[]
kdbins_a=[]
rhbins_a=[]


times_a=[]
histzzdr_a = []
edg_a = []

cbins = np.arange(-25,26,0.5)
dzbins = np.arange(-10,60,1)
drbins = np.arange(-2,4,0.1)
kdbins = np.arange(-2,4,0.1)
rrbins = np.logspace(0.01,100.01,30)


for d in dum:
    print d
#        dn = 
    tm=get_time(time_parse,d,wdate_format)
    vvar = xr.open_dataset(d)

    rdat = get_data(exper = exper,type=yp,mphys=mphys,date=date,file=vvar, tm = [tm],pol_on =1)
    
    dzcfad, hts, dzbins = GF.cfad(data = rdat.data[rdat.dz_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = dzbins,
                        ret_z=1,ret_bin = 1,thresh=-10,norm=False)

    dzcfad_a.append(dzcfad)
    dzbins_a.append(dzbins)

    bad = np.where(rdat.data[rdat.dz_name].data<-10.)
    rdat.data[rdat.zdr_name].data[bad]=np.nan
    
    drcfad, hts, drbins = GF.cfad(data = rdat.data[rdat.zdr_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = drbins,
                        ret_z=1,ret_bin = 1,thresh=-10,norm = False)
    drcfad_a.append(drcfad)
    drbins_a.append(drbins)

    kdcfad, hts, kdbins = GF.cfad(data = rdat.data[rdat.kdp_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = kdbins,
                        ret_z=1,ret_bin = 1,thresh=-10, norm = False)
    kdcfad_a.append(kdcfad)
    kdbins_a.append(kdbins)
    
    hts_a.append(hts)
    times_a.append(tm)
    
    hist,edg = GF.hist2d(rdat.data[rdat.zdr_name].data,rdat.data[rdat.dz_name].data,binsx = np.arange(-2,4.1,0.2),binsy=np.arange(-10,60,1))

    histzzdr_a.append(hist)
    edg_a.append(edg)
    vvar.close()
    
    
    rconf = RadarConfig.RadarConfig()
    rconf.date = times_a
    rconf.mphys = rdat.mphys
    print rconf.mphys
    rconf.exper = rdat.exper
    rconf.radar_name = rdat.radar_name
    del rdat

st = rconf.sav_title()

dzcfad_all = np.sum(dzcfad_a,axis=0)
fig, ax = GF.cfad_plot('DZ',cfad = dzcfad_all, hts = hts_a[0],  bins = dzbins_a[0],cfad_on = 0,rconf = rconf,tspan = times_a,maxval=20,cont=True,levels = True)
plt.savefig('{id}CFADDZ{s}.{t}'.format(id=image_dir,s=st,t=ptype),dpi=200)


drcfad_all = np.sum(drcfad_a,axis=0)
fig, ax = GF.cfad_plot('DR',cfad = drcfad_all, hts = hts_a[0],  bins = drbins_a[0],cfad_on = 0,rconf = rconf,tspan = times_a,maxval=20,cont=True,levels = True)
plt.savefig('{id}CFADDR{s}.{t}'.format(id=image_dir,s=st,t=ptype),dpi=200)

kdcfad_all = np.sum(kdcfad_a,axis=0)
fig, ax = GF.cfad_plot('DR',cfad = kdcfad_all, hts = hts_a[0],  bins = kdbins_a[0],cfad_on = 0,rconf = rconf,tspan = times_a,maxval=20,cont=True,levels = True)
plt.savefig('{id}CFADKD{s}.{t}'.format(id=image_dir,s=st,t=ptype),dpi=200)


