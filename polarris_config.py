"""
This is the configruation file for setting up iPOLARRIS. Herein, the speicifcs of the dataset need to be defined such as the experiment, location and type of reading, etc.
Written by Brenda Dolan
May 2017
bdolan@atmos.colostate.edu
"""

import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import RadarData
import datetime
import matplotlib.pyplot as plt
import numpy as np
import GeneralFunctions as GF
import RadarConfig

###############
def get_data(exper = 'TWPICE',tm=0,type='wrf',mphys='4ICE',date='2006123',file=r'wrf_twp_files.txt',pol_on = False):

    ############MC3E radar#######
    if exper == 'MC3E':
        if type == 'obs':
            radarname = 'CSAPR'
            band = 'C'
            read_list = 'False'
#            radar_files = '/Volumes/rawah/data/jib/MC3E_darwin_analysis/csapr/{d}/ncfiles/'.format(d=date)
            dd_files = '/Volumes/rawah/data/jib/MC3E_darwin_analysis/csapr/{d}/cdffiles/'.format(d=date)
            rdate_format = '%Y%m%d_%H%M%S'
            ddate_format = '%Y%m%d_%H%M'
            dz_name = 'DBZCS'
            dr_name = 'ZDRCS'
            kd_name = 'KDPCS'
            t_name= None
            rh_name = 'RHOCS'
            uname = None
            vname = None
            wname = None
            xname = 'x'
            yname = 'y'
            zname = 'z'
#            w_name = 'Wvar'
            doff = 0
            ddoff = 0
            ddadd= 13
            dext = 'cdf'
            ext='nc'
            dd_on = 'True'
            wdate_format=rdate_format
            usedate_range = 'True'
            time_parse=[doff,ddadd]
        elif type == 'wrf':
            radarname='WRF_Cband'
            read_list = 'True'
#            radar_files = r'wrf_mc3e_files.txt'
            usedate_range = 'False'
            type  = 'wrf'
            dd_on = 'False'
            dz_name = 'zhh01'
            dr_name = 'zdr01'
            kd_name='kdp01'
            rh_name='rhohv01'
            t_name='t_air'
            uname = 'u'
            vname = 'v'
            wname = 'w'
            xname = 'longitude'
            yname = 'latitude'
            zname='hgt'
            band = 'C'
            mphys=mphys
            time_parse=[11,19]
            wdate_format ='%Y-%m-%d_%H-%M-%S'


        else:
            print 'I do not recognize the type'

        sfiles = '/Volumes/rawah/data/jib/MC3E_darwin_analysis/'
        sdate_format = '%Y%m%d%H'
        sstat = 'LMN'
        dates1 = date
        ad = ''
        lat = 36.79616
        lon = -97.450546
        alt = 0.327


    ########TWP-ICE############
    elif exper == 'TWPICE':
        if type == 'obs':
            radarname='CPOL'
            band = 'C'
            read_list = 'False'
#            radar_files = '/Volumes/rawah/data/jib/twp_ice/cpol/{d}/gridded/{d}/'.format(d=date)
            dd_files = '/Volumes/rawah/data/jib/twp_ice/dualdoppler_files/{d}/DUALDOP/'.format(d=date)
            sfiles = '/Volumes/rawah/data/jib/MC3E_darwin_analysis/'
            rdate_format = '%Y%m%d_%H%M%S'
            ddate_format = '%Y%m%d_%H%M'
            sdate_format = '%Y%m%d%H'
            doff = 4
            ddoff = 12
            ddadd = 4
            dext = 'nc'
            ext='nc'
            wdate_format=rdate_format
            dates1 = date
            ad = dates1+'_'
            dz_name = 'DZ'
            dr_name = 'CR'
            t_name= None
            kd_name = 'KD'
            rh_name = 'RH'
#            w_name = 'Wvar'
            uname = None
            vname = None
            wname = None
            xname = 'x'
            yname = 'y'
            zname = 'z'
            sstat = 'DWN'
            time_parse=[doff,ddoff]

        elif type == 'wrf':
            radarname='WRF_Cband'
            read_list = 'True'
#            radar_files = r'wrf_twpice_files.txt'
            usedate_range = 'False'
            type  = 'wrf'
            dd_on = 'False'
            dz_name = 'zhh01'
            dr_name = 'zdr01'
            kd_name='kdp01'
            rh_name='rhohv01'
            t_name='t_air'
            uname = 'u'
            vname = 'v'
            wname = 'w'
            xname = 'longitude'
            yname = 'latitude'
            zname='hgt'
            band = 'C'
            time_parse=[11,19]
            wdate_format ='%Y-%m-%d_%H-%M-%S'

        else:
            print 'I do not recognize the type'

        lat = -12.24917
        lon = 131.04444
        dd_on='True'
        usedate_range = 'True'

    else:
        print 'Experiment not in database'

    mphys=mphys
    
#    print mphys

    ###########################Now grab the radar data##############

    def foo(s1):
        return '{}'.format(s1.rstrip())

    rdata = RadarData.RadarData(file,tm,dz = dz_name,zdr=dr_name,
                                              kdp=kd_name,rho=rh_name,temp=t_name,
                                              u=uname,v=vname,w=wname,x=xname,
                                              y=yname,z=zname,lat=lat, lon=lon,band = band,
                                              exper=exper,mphys=mphys,radar_name = radarname)
                                              
    #vvar.close()
    #rdata.radar_name = radarname
    if type == 'wrf':
        rdata.convert_t()
    print 'Calculating polarimetric fields like HID and rain...'
    if pol_on == True:
        rdata.calc_pol_analysis()

    return(rdata)
    
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
    
        hidhts,hiddum =GF.hid_cdf(rdat.hid, rdat.data[rdat.z_name][:].data,rdat.species,z_resolution=1.0, pick=None,z_ind =0)
        print np.shape(hiddum)
        hid_hts.append(hidhts)
        hid_cfad.append(hiddum)
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
