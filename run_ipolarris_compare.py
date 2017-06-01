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

#mc3e_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_mc3e.txt'
mc3e_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt'
mc3e_yp = 'wrf'
mc3e_exper = 'MC3E'
mc3e_mphys = '4ICE'
mc3e_date = '20110523'
mc3e_time_parse=[11,19]
#mc3e_wdate_format = '%Y%m%d_%H%M%S'
mc3e_wdate_format = '%Y-%m-%d_%H-%M-%S'
#mc3e_wdate_format = '%Y-%m-%d_%H:%M:%S'


#twpice_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_twpice.txt'
twpice_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_twp_files.txt'
twpice_yp = 'wrf'
twpice_exper = 'TWPICE'
twpice_mphys = '4ICE'
twpice_date = '2006123'
twpice_time_parse=[11,19]
twpice_wdate_format = '%Y-%m-%d_%H-%M-%S'
#twpice_wdate_format = '%Y-%m-%d_%H:%M:%S'
#twpice_wdate_format = '%Y%m%d_%H%M%S'


ptype = 'png'
#image_dir = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/'
image_dir = r'/Users/bdolan/scratch/iPOLARRIS_images_test/'

flags = {'cfad_4panel_flag':False,  # 4 panel CFADs of Dbz, Zdr, Kdp and W
         'hid_cfad_flag':False,     # HID CFAD
         'joint_flag':False,        #Z vs. Zdr
         'cfad_individ_flag':False,  #Make separate images for dbz, Zdr, Kdp, W and Rho
         'hid_prof':False,          #Profile of HID species with height
         'all_cappi':False,         # make a cappi cross section for all times. change parameters in plot_driver.py
         'all_xsec':True,          # make RHI cross sections for all times. change parmaters in plot_driver.py
         'up_width':False,          # Make a plot of the vertical velocity profile with temperature.
         'qr_cappi':False,          # Make cappi cross section of mixing ratios. change parameters in plot_driver.py
         'qr_rhi':False}            # make RHI cross sections of mixing ratios. change parameters in plot_driver.py

plot_single = 0         #Set this flag to 1, plus the flags above, for individual plots at each time.
plot_int = 1           #Set this flag for integrated time plots over all files in the list.
plot_diff = 0            #set this flag to compare two different datasets in the integrated plots.

mc3e_dat = run_exper(mc3e_radar_files,mc3e_exper,mc3e_mphys,mc3e_date,mc3e_time_parse,mc3e_wdate_format,mc3e_yp,plot_on=plot_single,flags=flags,image_dir=image_dir)
twpice_dat = run_exper(twpice_radar_files,twpice_exper,twpice_mphys,twpice_date,twpice_time_parse,twpice_wdate_format,twpice_yp,plot_on=1,flags=flags,image_dir=image_dir)


dat1 = mc3e_dat
dat2 = twpice_dat

if plot_int == 1:

    #######Make plots of Integrated CFADS 
    plot_driver.plot_cfad_int(dat1,typ='dz',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_int(dat1,typ='dr',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_int(dat1,typ='kd',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_int(dat1,typ='w',image_dir = image_dir,ptype=ptype)
    

    #######Now HID##############
    plot_driver.plot_hid_int(dat1,typ='hid',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_hid_prof_int(dat1,typ='hid',image_dir=image_dir,ptype=ptype)

    ########Now 2D histograms######
    plot_driver.plot_joint_int(dat1,typ='zzdr',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_int(dat1,typ='zkdp',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_int(dat1,typ='zw',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_int(dat1,typ='wr',image_dir=image_dir,ptype=ptype)
    ########Updraft Width##########
    plot_driver.plot_upwidth_int(dat1,image_dir=image_dir,ptype=ptype)


if plot_diff == 1:

    ############PLOTTING Differences#################

    ########Make plots of CFAD trios (dat1 , dat2 , difference)##############
    plot_driver.plot_cfad_compare(dat1,dat2,typ='dz',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_compare(dat1,dat2,typ='dr',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_compare(dat1,dat2,typ='kd',image_dir = image_dir,ptype=ptype)
    plot_driver.plot_cfad_compare(dat1,dat2,typ='w',image_dir = image_dir,ptype=ptype)

    #########Now HID##############
    plot_driver.plot_hid_2panel(dat1,dat2,typ='hid',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_hid_profile(dat1,dat2,typ='hid',image_dir=image_dir,ptype=ptype)
# 
    ########Now 2D histograms######
    plot_driver.plot_joint_comp(dat1,dat2,typ='zzdr',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_comp(dat1,dat2,typ='zkdp',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_comp(dat1,dat2,typ='zw',image_dir=image_dir,ptype=ptype)
    plot_driver.plot_joint_comp(dat1,dat2,typ='wr',image_dir=image_dir,ptype=ptype)

    ########Updraft Width##########
    plot_driver.plot_upwidth(dat1,dat2,image_dir=image_dir,ptype=ptype)
