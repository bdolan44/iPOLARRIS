"""
run_ipoloarris_single.py

Brenda Dolan, CSU, May 2017
bdolan@atmos.colostate.edu

This program runs a single set of files defined in the radar_files list and runs them through
the iPOLARRIS radar calculations (HID, rain, etc.) and plots them.

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
plot_single = 1         #Set this flag to 1, plus the flags above, for individual plots at each time.
plot_int = 1           #Set this flag for integrated time plots over all files in the list.

ext = 'test' ##Use this variable to label any sensitivity studies, like changing axis ratio or whatever.

if disc==1:
    image_dir = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/'
    wdate_format = '%Y-%m-%d_%H:%M:%S'
else:
    image_dir = r'/Users/bdolan/scratch/iPOLARRIS_images_test/'
    wdate_format = '%Y-%m-%d_%H-%M-%S'


#radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_mc3e.txt'
radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt'
yp = 'wrf'
exper = 'MC3E'
date = '20110523'
mphys = '4ICE'

#radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_twpice.txt'
#radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_twpice_files.txt'
#exper = 'TWPICE'
#date = '20060123'
time_parse=[11,19]

#mphys = 'SBM'


flags = {'cfad_4panel_flag':True,  # 4 panel CFADs of Dbz, Zdr, Kdp and W
         'hid_cfad_flag':True,     # HID CFAD
         'joint_flag':True,        #Z vs. Zdr
         'cfad_individ_flag':True,  #Make separate images for dbz, Zdr, Kdp, W and Rho
         'hid_prof':True,          #Profile of HID species with height
         'all_cappi':True,         # make a cappi cross section for all times. change parameters in plot_driver.py
         'all_xsec':True,          # make RHI cross sections for all times. change parmaters in plot_driver.py
         'up_width':True,          # Make a plot of the vertical velocity profile with temperature.
         'qr_cappi':True,          # Make cappi cross section of mixing ratios. change parameters in plot_driver.py
         'qr_rhi':True}            # make RHI cross sections of mixing ratios. change parameters in plot_driver.py

###########Run through each file in the radar_files list and plot if that is turned on########

dat1 = run_exper(radar_files,exper,mphys,date,time_parse,wdate_format,yp,plot_on=plot_single,flags=flags,image_dir=image_dir)

if plot_int == 1:
    #######Make plots of Integrated CFADS 
    plot_driver.plot_cfad_int(dat1,typ='dz',image_dir = image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_cfad_int(dat1,typ='dr',image_dir = image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_cfad_int(dat1,typ='kd',image_dir = image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_cfad_int(dat1,typ='w',image_dir = image_dir,ptype=ptype,extra=ext)
    

    #######Now HID##############
    plot_driver.plot_hid_int(dat1,typ='hid',image_dir=image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_hid_prof_int(dat1,typ='hid',image_dir=image_dir,ptype=ptype,extra=ext)

    ########Now 2D histograms######
    plot_driver.plot_joint_int(dat1,typ='zzdr',image_dir=image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_joint_int(dat1,typ='zkdp',image_dir=image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_joint_int(dat1,typ='zw',image_dir=image_dir,ptype=ptype,extra=ext)
    plot_driver.plot_joint_int(dat1,typ='wr',image_dir=image_dir,ptype=ptype,extra=ext)
    ########Updraft Width##########
    plot_driver.plot_upwidth_int(dat1,image_dir=image_dir,ptype=ptype,extra=ext)

