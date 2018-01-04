"""
run_ipoloarris_compare.py

Brenda Dolan, CSU, May 2017
bdolan@atmos.colostate.edu

This program runs through two different file lists and generates plots and difference plots between the two. This could be between simulations
and observations, between locations, or between microphysics schemes.

Several other helper functions are used:
- polarris_config.py defines some of the variable names. If you want to change to a different radar, you will need ot modify that program
- plot_driver.py houses many of the plotting functions. If you want to change the cross-section plots, bounds used for a project, etc, that program will need to be modified.

-RadarData.py is the class fuction that houses the data, and is where the plotting for single times is done.
-RadarConfig.py defines many of the variable ranges and other more "global" attributes.
-GeneralFunctions also has some general functions (ha) and houses the plotting code for time-integrated plots.

"""

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

#######################################################################
#######change the parameters below to run different cases #############
#######################################################################

disc =0
ptype = 'png'
ext = 'regdiff' ##Use this variable to label any sensitivity studies, like changing axis ratio or whatever.

if disc==1:
    image_dir = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/'
    d1_wdate_format = '%Y-%m-%d_%H:%M:%S'
    d2_wdate_format = '%Y-%m-%d_%H:%M:%S'
else:
    image_dir = r'/Users/bdolan/scratch/iPOLARRIS_images_test/'
    d1_wdate_format = '%Y-%m-%d_%H-%M-%S'
    d2_wdate_format = '%Y-%m-%d_%H-%M-%S'

#mc3e_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_mc3e.txt'
d1_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt'
d1_yp = 'wrf'
d1_exper = 'MC3E'
d1_mphys = '4ICE'
d1_date = '20110523'
d1_time_parse=[11,19]
#mc3e_wdate_format = '%Y%m%d_%H%M%S'


#twpice_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_twpice.txt'
d2_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_twp_files.txt'
d2_yp = 'wrf'
d2_exper = 'TWPICE'
d2_mphys = '4ICE'
d2_date = '2006123'
d2_time_parse=[11,19]


dat1 = run_exper(d1_radar_files,d1_exper,d1_mphys,d1_date,d1_time_parse,d1_wdate_format,d1_yp,image_dir=image_dir)
dat2 = run_exper(d2_radar_files,d2_exper,d2_mphys,d2_date,d2_time_parse,d2_wdate_format,d2_yp,image_dir=image_dir)


############PLOTTING Differences#################

########Make plots of CFAD trios (dat1 , dat2 , difference)##############
plot_driver.plot_cfad_compare(dat1,dat2,typ='dz',image_dir = image_dir,ptype=ptype,extra=ext)
plot_driver.plot_cfad_compare(dat1,dat2,typ='dr',image_dir = image_dir,ptype=ptype,extra=ext)
plot_driver.plot_cfad_compare(dat1,dat2,typ='kd',image_dir = image_dir,ptype=ptype,extra=ext)
plot_driver.plot_cfad_compare(dat1,dat2,typ='w',image_dir = image_dir,ptype=ptype,extra=ext)

#########Now HID##############
plot_driver.plot_hid_2panel(dat1,dat2,typ='hid',image_dir=image_dir,ptype=ptype,extra=ext)
plot_driver.plot_hid_profile(dat1,dat2,typ='hid',image_dir=image_dir,ptype=ptype,extra=ext)
# 
########Now 2D histograms######
plot_driver.plot_joint_comp(dat1,dat2,typ='zzdr',image_dir=image_dir,ptype=ptype,extra=ext)
plot_driver.plot_joint_comp(dat1,dat2,typ='zkdp',image_dir=image_dir,ptype=ptype,extra=ext)
plot_driver.plot_joint_comp(dat1,dat2,typ='zw',image_dir=image_dir,ptype=ptype,extra=ext)
plot_driver.plot_joint_comp(dat1,dat2,typ='wr',image_dir=image_dir,ptype=ptype,extra=ext)

########Updraft Width##########
plot_driver.plot_upwidth(dat1,dat2,image_dir=image_dir,ptype=ptype,extra=ext)
