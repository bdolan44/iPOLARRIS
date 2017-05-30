import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import numpy as np
import RadarData
import datetime
#import matplotlib.pyplot as plt

from polarris_config import get_data
from plot_driver import make_single_pplots
from plot_driver import make_compare_pplots
#######################################################
def get_time(time_parse,filename,dformat):
    tp = time_parse[0]
    te = time_parse[1]
    print time_parse
#    print filename
    base = os.path.basename(filename)
    print base
    radar_name = base[:3]

    radcdate=np.str(base[tp:tp+te])
    print radcdate
    date=datetime.datetime.strptime(radcdate,dformat)
    print date
    return date

def foo(s1):
    return '{}'.format(s1.rstrip())
#######################################################


# mc3e_wrf = get_data(exper='MC3E',radar_files=r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt')
radar_files = r'/Users/bdolan/scratch/POLARRIS_2/cpol_twp_files.txt'
type = 'obs'
exper = 'TWPICE'
mphys = '4ICE'
date = '2006123'


flags = {'cfad_4panel_flag':False,
         'hid_cfad_flag':False,
         'joint_flag':True,
         'cfad_individ_flag':False,
         'hid_prof':False,
         'all_cappi':False,
         'all_xsec':False,
         'up_width':False,
         'qr_cappi':False,
         'qr_rhi':False,}


dum =[]
with open(radar_files) as f: 
    for line in f:
        dat = (line)
        dum.append(foo(dat))
mydat= []
tm = []

if flags['cfad_4panel_flag'] == True:
    dzcfad = []
    drcfad = []
    kdcfad = []
    rhcfad = []
    rrcfad = []
    


for d in dum:
    print d
#        dn = 
    tm=get_time(time_parse,d,wdate_format)
    vvar = xr.open_dataset(dum,concat_dim='t')
    #vvar.xr.Dataset({'time': datetime.datetime(2000, 1, 1)})
#        mydat.append(vvar)            
#    vvar.close()

    if flags['hid_cfad_flag'] == True or flags['hid_prof'] == True or flags['call_cappi'] == True or flags['all_xsec'] == True:
        pol_on == True
    else:
        pol_on == False
    
    
    rdat = get_data(exper = exper,type=type,mphys=mphys,date=date,file=vvar, tm = tm, pol_on = pol_on)


#         
    if flags['cfad_4panel_flag'] == True:
        cbins = np.arange(-25,26,0.5)
        dzbins = np.arange(-10,70,1)
        drbins = np.arange(-2,6,0.1)
        kdbins = np.arange(-2,4,0.05)
        rrbins = np.logspace(0.01,100.01,30)

        dzcfad.append(rdat.cfad(rdat.dz_name,bins=cbins,z_resolution=1.0,levels='levs'))
        drcfad.append(rdat.cfad(rdat.dr_name,bins=cbins,z_resolution=1.0,levels='levs'))
        kdcfad.append(rdat.cfad(rdat.kd_name,bins=cbins,z_resolution=1.0,levels='levs'))
        wcfad.append(rdat.cfad(rdat.w_name,bins=cbins,z_resolution=1.0,levels='levs'))
        rhocfad.append(rdat.cfad(rdat.rho_name,bins=cbins,z_resolution=1.0,levels='levs'))
        rrcfad.append(rdat.cfad(rdat.rr_name,bins=cbins,z_resolution=1.0,levels='levs'))
    rdat.close()
