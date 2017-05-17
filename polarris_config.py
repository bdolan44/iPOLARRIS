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
            radarname='WRF Cband'
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
            radarname='WRF Cband'
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
    #rdata.convert_t()
    print 'Calculating polarimetric fields like HID and rain...'
    if pol_on == True:
        rdata.calc_pol_analysis()

    return(rdata)