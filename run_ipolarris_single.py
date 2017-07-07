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

def match_dd(rdate,ddates):
    dum=abs(rdate-np.array(ddates))

    try:
        mval=np.argwhere(dum == np.min(dum))
        diff=np.min(dum)
        if diff.total_seconds() < 600.:
                return mval
        else:
                return 'no'
    except ValueError: 
        print 'Something DD is not working here!'
        
def match_snd(rdate,sdates):
    dum=abs(rdate-np.array(sdates))

    try:
        mval=np.argwhere(dum == np.min(dum))
        diff=np.min(dum)
        #allow 12 hours between the radar obs and the sounding
        if diff.total_seconds() < 43200:
                return mval
        else:
                return 'no'
    except ValueError: 
        print 'Something sound is not working here!'


def find_dd_match(config):
    rdum =[]
    with open(config.radar_files) as f: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    ddum =[]

    with open(config.dd_files) as g: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in g:
            dat = (line)
            ddum.append(foo(dat))

   ddates = []
   for v,dname in enumerate(ddum):

        base = os.path.basename(dname)
    
        radcdate=np.str(base[config.ddoff:config.ddoff+config.ddadd])
        #print radcdate
        #ddate_format='%Y%m%d_%H%M'
        if config.ad != '':
            dates=datetime.datetime.strptime(config.dates1+'_{r}'.format(a=config.ad,r=radcdate),config.ddate_format)
        else:
            dates=datetime.datetime.strptime('{r}'.format(a=config.ad,r=radcdate),config.ddate_format)
        ddates.append(dates)

    slist = sorted(glob.glob('{p}*{s}_*.txt'.format(p=config.sfiles,s=config.sstat)))


    sdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
#            print base
        radcdate=np.str(base[13:13+9])

        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config.sdate_format)
        sdates.append(dates)

    radlist=[]

mdfiles = {}
msfiels = {}
    for v,cname in enumerate(rdum):
#            print cname
        base = os.path.basename(cname)
        radcdate=np.str(base[config.doff:config.doff+15])
        dates=datetime.datetime.strptime(radcdate,config.rdate_format)
        #print dates, etime,stime
        if (dates <= config.etime) and (dates >= config.stime):
            #print cname
            mval = match_dd(dates,ddates)
            #print dates,ddates
        
            if mval != 'no':
                dfile = dlist[mval]
                print 'Found DD match!', dfile
                

            #now find a sounding match
            mv = match_snd(dates,sdates)
            if mv != 'no':
                print slist[mv]
                snd = SkewT.Sounding(slist[mv])
                radar.add_sounding_object(snd) # this will add the sounding object to the radar object
                            # and then will take the heights and temps
                radar.interp_sounding()
            radar.calc_pol_analysis()
            if radar.date is 'None':
                radar.date = date
            radlist.append(radar)
        
    return radlist



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

config = {}
with open("polarris_config.txt") as f:
    for line in f:
        #print line
        if not line.startswith("#"):
            key, val, comment = line.split('==')
            config[(key.replace(" ", ""))] = val.replace(" ", "")

if config.dd_on == True:
    


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

