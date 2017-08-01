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
import sys

import dill as pickle

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
    with open(config['radar_files']) as f: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    ddum =[]

    with open(config['dd_files']) as g: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in g:
            dat = (line)
            ddum.append(foo(dat))

    ddates = []
    for v,dname in enumerate(ddum):

        base = os.path.basename(dname)

        radcdate=np.str(base[config['ddoff']:config['ddoff']+config['ddadd']])
#        print radcdate, dname
        #ddate_format='%Y%m%d_%H%M'
        if eval(config['ad']) != '':
#            print radcdate
            dates=datetime.datetime.strptime(config['date']+'_{r}'.format(a=config['ad'],r=radcdate),config['ddate_format'])
        else:
            dates=datetime.datetime.strptime('{r}'.format(a=config['ad'],r=radcdate),config['ddate_format'])
        ddates.append(dates)


    radlist=[]

    mdfiles = {}
    for v,cname in enumerate(rdum):
        #print cname
        base = os.path.basename(cname)
        radcdate=np.str(base[config['doff']:config['doff']+15])
#        print radcdate
        dates=datetime.datetime.strptime(radcdate,config['rdate_format'])
        #print dates, etime,stime
        if (dates >= config['etime']) and (dates <= config['stime']):
            #print cname
            mval = match_dd(dates,ddates)
            #print dates,ddates
        
            if mval != 'no':
                dfile = ddum[mval]
                #print 'Found DD match!', dfile
                mdfiles[cname] = dfile
            else:
                mdfiles[cname] = None

        
    return mdfiles

def find_snd_match(config):
    rdum =[]
    with open(config['radar_files']) as f: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in f:
            dat = (line)
            rdum.append(foo(dat))

    slist = sorted(glob.glob('{p}*{s}_*.txt'.format(p=config['sfiles'],s=(config['sstat']))))
    sdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
#            print base
        radcdate=np.str(base[13:13+9])

        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config['sdate_format'])
        sdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
#            print cname
        base = os.path.basename(cname)
        radcdate=np.str(base[config['doff']:config['doff']+15])
        dates=datetime.datetime.strptime(radcdate,config['rdate_format'])
        #print dates, etime,stime
        if (dates >= config['etime']) and (dates <= config['stime']):
            #print cname
            #now find a sounding match
            mv = match_snd(dates,sdates)
            if mv != 'no':
#                print 'sounding match',slist[mv]
                msfiles[cname] = slist[mv]
            else:
                return None
        
    return msfiles


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


def hasNumbers(inputString):
     return any(char.isdigit() for char in inputString)

#######################################################################
#######change the parameters below to run different cases #############
#######################################################################
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
            if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'mphys' or key.replace(" ","") == 'vr_name' or key.replace(" ","") == 'vdop_name' or key.replace(" ","") =='u_name' or key.replace(" ","") =='v_name':
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
                        

#print config['dd_on']
if config['dd_on'] == True:
    dmatch = find_dd_match(config)
#    print 'dmatch',dmatch
else:
    dmatch = None

if config['snd_on'] == True:
    smatch = find_snd_match(config)
else:
    smatch = None

#print 'smatch',smatch

dat1 = run_exper(config, dmatch = dmatch, smatch=smatch,interactive = False)

if not sys.argv[2:]:
    if config['plot_int'] == 1:
        #######Make plots of Integrated CFADS 
        if config['plot_cs'] == 1:
            hold = config['extra']
            config['extra']='{e}_convective'.format(e=hold)
            plot_driver.plot_cfad_int(dat1,config,typ='dzc')
            plot_driver.plot_cfad_int(dat1,config,typ='drc')
            plot_driver.plot_cfad_int(dat1,config,typ='kdc')
            if dat1['runw'] is True:
                plot_driver.plot_cfad_int(dat1,config,typ='wc')
            config['extra']=hold

            hold = config['extra']
            config['extra']='{e}_stratiform'.format(e=hold)

            plot_driver.plot_cfad_int(dat1,config,typ='dzs')
            plot_driver.plot_cfad_int(dat1,config,typ='drs')
            plot_driver.plot_cfad_int(dat1,config,typ='kds')
            if dat1['runw'] is True:
                plot_driver.plot_cfad_int(dat1,config,typ='ws')
            config['extra']=hold

        plot_driver.plot_cfad_int(dat1,config,typ='dz')
        plot_driver.plot_cfad_int(dat1,config,typ='dr')
        plot_driver.plot_cfad_int(dat1,config,typ='kd')
        if dat1['runw'] is True:
            plot_driver.plot_cfad_int(dat1,config,typ='w')
    
        #######Now HID##############
        plot_driver.plot_hid_int(dat1,config,typ='hid')
        if config['plot_cs'] == 1:
            hold = config['extra']
            config['extra']='{e}_convective'.format(e=hold)
            plot_driver.plot_hid_int(dat1,config,typ='hidc')
            config['extra']=hold

            hold = config['extra']
            config['extra']='{e}_stratiform'.format(e=hold)
            plot_driver.plot_hid_int(dat1,config,typ='hids')
            config['extra']=hold

        plot_driver.plot_hid_prof_int(dat1,config,typ='hid')

        ########Now 2D histograms######
        plot_driver.plot_joint_int(dat1,config,typ='zzdr')
        plot_driver.plot_joint_int(dat1,config,typ='zkdp')
        if dat1['runw'] is True:
            plot_driver.plot_joint_int(dat1,config,typ='zw')
            ########Updraft Width##########
            plot_driver.plot_upwidth_int(dat1,config)
            if dat1['runrr'] is True:
                plot_driver.plot_joint_int(dat1,config,typ='wr')

else:
    configfile1 = sys.argv[2:]
    config1 = {}
    #print sys.argv[1:]

    with open(configfile1[0]) as f:
        for line in f:
            #print line
            if not line.startswith("#"):
                key, val, comment = line.split('==')
                vval = val.replace(" ","")
                numck = hasNumbers(vval)
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'mphys' or key.replace(" ","") == 'vr_name' or key.replace(" ","") == 'vdop_name' or key.replace(" ","") =='u_name' or key.replace(" ","") =='v_name':
                    numck = False
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'extra' or key.replace(" ", "") == 'ptype':
                    vval = vval.strip("''")
                #print numck
                #print vval,key
                if key.replace(" ", "") == 'image_dir':
                    numck = True

                if numck is True or vval == 'None' or vval == 'True' or vval == 'False':
                    try:
                        config1[(key.replace(" ", ""))] = eval(vval)
                    except:
                        if "datetime" in vval:
                            config1[(key.replace(" ", ""))] = vval
                else:
                    config1[(key.replace(" ", ""))] = vval
            
    #print config['dd_on']
    if config1['dd_on'] == True:
        dmatch1 = find_dd_match(config1)
    #    print 'dmatch',dmatch
    else:
        dmatch1 = None

    if config['snd_on'] == True:
        smatch1 = find_snd_match(config1)
    else:
        smatch1 = None

    #print 'smatch',smatch

    dat2 = run_exper(config1, dmatch = dmatch1, smatch=smatch1,interactive = False)


    ############PLOTTING Differences#################

    ########Make plots of CFAD trios (dat1 , dat2 , difference)##############
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='dz')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='dr')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='kd')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='w')

    
    ########Make plots of CFAD trios (dat1 , dat2 , difference)##############
    hold = config['extra']
    config['extra']='{e}_convective'.format(e=hold)
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='dzc')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='drc')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='kdc')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='wc')
    config['extra']=hold

    ########Make plots of CFAD trios (dat1 , dat2 , difference)##############
    hold = config['extra']
    config['extra']='{e}_stratiform'.format(e=hold)
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='dzs')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='drs')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='kds')
    plot_driver.plot_cfad_compare(dat1,dat2,config,typ='ws')
    config['extra']=hold


    #########Now HID##############
    plot_driver.plot_hid_2panel(dat1,dat2,config,typ='hid')
    plot_driver.plot_hid_profile(dat1,dat2,config,typ='hid')

    hold = config['extra']
    config['extra']='{e}_convective'.format(e=hold)
    plot_driver.plot_hid_2panel(dat1,dat2,config,typ='hidc')
    plot_driver.plot_hid_profile(dat1,dat2,config,typ='hidc')
    config['extra']=hold

    hold = config['extra']
    config['extra']='{e}_stratiform'.format(e=hold)
    plot_driver.plot_hid_2panel(dat1,dat2,config,typ='hids')
    plot_driver.plot_hid_profile(dat1,dat2,config,typ='hids')
    config['extra']=hold

    # 
    ########Now 2D histograms######
    plot_driver.plot_joint_comp(dat1,dat2,config,typ='zzdr')
    plot_driver.plot_joint_comp(dat1,dat2,config,typ='zkdp')
    plot_driver.plot_joint_comp(dat1,dat2,config,typ='zw')
#    plot_driver.plot_joint_comp(dat1,dat2,config,typ='wr')

    ########Updraft Width##########
    plot_driver.plot_upwidth(dat1,dat2,config)


    pickle.dump(dat1, open( "{e}_{x}.p".format(e=config['exper'],x=config['extra']), "wb" ) )
    pickle.dump(dat2, open( "{e}_{x}.p".format(e=config1['exper'],x=config1['extra']), "wb" ) )

