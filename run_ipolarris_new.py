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

def run_exper(radar_files,exper,mphys,date,time_parse,wdate_format,yp):

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
    warea=[]
    histzzdr_a = []
    edg_a = []

    cbins = np.arange(-25,26,0.5)
    dzbins = np.arange(-10,60,1)
    drbins = np.arange(-2,6,0.05)
    kdbins = np.arange(-2,6,0.05)
    rrbins = np.logspace(0.01,100.01,30)

    water_vert = []
    graup_vert =[]
    hail_vert = []
    snow_vert = []
    hid_hts = []
    
    wcfad_a =[]
    wbins_a =[]

    hid_cfad = []

    for d in dum:
        print d
    #        dn = 
        tm=get_time(time_parse,d,wdate_format)
        vvar = xr.open_dataset(d)

        rdat = get_data(exper = exper,type=yp,mphys=mphys,date=date,file=vvar, tm = [tm],pol_on =1)
        bad = np.where(rdat.data[rdat.dz_name].data<-10.)
        rdat.data[rdat.dz_name].data[bad]=np.nan
    
        dzcfad, hts, dzbins = GF.cfad(data = rdat.data[rdat.dz_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = dzbins,
                            ret_z=1,ret_bin = 1,thresh=-10)

        dzcfad_a.append(dzcfad)
        dzbins_a.append(dzbins)

        rdat.data[rdat.zdr_name].data[bad]=np.nan
    
        drcfad, hts, drbins = GF.cfad(data = rdat.data[rdat.zdr_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = drbins,
                            ret_z=1,ret_bin = 1,thresh=-10)
        drcfad_a.append(drcfad)
        drbins_a.append(drbins)

        kdcfad, hts, kdbins = GF.cfad(data = rdat.data[rdat.kdp_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = kdbins,
                            ret_z=1,ret_bin = 1,thresh=-10)
        kdcfad_a.append(kdcfad)
        kdbins_a.append(kdbins)
    

        wcfad, hts, wbins = GF.cfad(data = rdat.data[rdat.w_name].data,hts = rdat.data[rdat.z_name][:].data,value_bins = cbins,
                            ret_z=1,ret_bin = 1,thresh=-10)
        wcfad_a.append(wcfad)
        wbins_a.append(wbins)

        hts_a.append(hts)
        times_a.append(tm)
    
        hist,edg = GF.hist2d(rdat.data[rdat.zdr_name].data,rdat.data[rdat.dz_name].data,binsx = np.arange(-2,4.1,0.2),binsy=np.arange(-10,60,1))

        histzzdr_a.append(hist)
        edg_a.append(edg)
        vvar.close()
    
        hidwater = [1,2,10]
        hidgraup = [7,8]
        hidhail = [9]
        hidsnow =[3,4,5,6]
        
        rdat.hid[bad]=-1
        
        
        hidhts, mwrf_water_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,hidwater,rdat.species,z_resolution =1.0,z_ind=1)
        hidhts, mwrf_graup_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,hidgraup,rdat.species,z_resolution =1.0,z_ind=1)
        hidhts, mwrf_hail_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,hidhail,rdat.species,z_resolution =1.0,z_ind=1)
        hidhts, mwrf_snow_vert = GF.hid_vertical_fraction(rdat.hid,rdat.data[rdat.z_name][:].data,hidsnow,rdat.species,z_resolution =1.0,z_ind=1)

        water_vert.append(mwrf_water_vert)
        graup_vert.append(mwrf_graup_vert)
        hail_vert.append(mwrf_graup_vert)
        snow_vert.append(mwrf_graup_vert)
    
        hid_hts.append(hidhts)
        hid_cfad.append(GF.hid_cdf(rdat.hid, rdat.data[rdat.z_name][:].data,rdat.species,z_resolution=1.0, pick=None,z_ind =0))
        print type(rdat.T.data)
#         tmp, m_warea_wrf = GF.updraft_width_profile(rdat.data[rdat.w_name].data,rdat.data[rdat.w_name][:].data,thresh=5.0, temps=np.arange(20,-60,-5),\
#             z_ind=0,tcoord = True,temp = rdat.T.data)
#         warea_wrf = m_warea_wrf*rdat.dx*rad.dy/rdat.ntimes
#         if rdat.data[rdat.x_name].units == "[deg]":
#             ware_wrf=warea_wrf*110.*110.
# 
# 
#         warea.append(warea_wrf)
        
        rconf = RadarConfig.RadarConfig()
        rconf.date = times_a
        rconf.mphys = rdat.mphys
        rconf.exper = rdat.exper
        rconf.radar_name = rdat.radar_name
        del rdat


    dat = {'rconf': rconf,
          'dzcfad':dzcfad_a,
          'drcfad':drcfad_a,
          'kdcfad':kdcfad_a,
          'wcfad':wcfad_a,
          'hts':hts_a,
          'histzzdr':histzzdr_a,
          'edgzzdr':edg_a,
          'hidhts':hid_hts,
          'water_vert':water_vert,
          'graup_vert':graup_vert,
          'hail_vert':hail_vert,
          'snow_vert':snow_vert,
#           'warea':warea,
          'time':times_a,
          'dzbins':dzbins,
          'drbins':drbins,
          'kdbins':kdbins,
          'wbins':cbins,
          'hidcfad':hid_cfad}
    return dat

