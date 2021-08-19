"""
This is the configruation file for setting up iPOLARRIS. Herein, the specifics of the dataset need to be defined such as the experiment, location and type of reading, etc.
Written by Brenda Dolan (CSU) and Anthony Di Stefano (UBC)
Released: May 2017
Last Modified: June 2021
bdolan@atmos.colostate.edu
"""
from __future__ import print_function

import datetime
import glob
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import re
import sys
import xarray as xr
import time

import RadarData
import GeneralFunctions as GF
import RadarConfig
import plot_driver
from skewPy import SkewT

def fix_my_data(ds):
    return(ds.drop(['VTZCS','CVECS']))

def find_dd_match(rdum,ddum,rdate,ddates):

    radlist=[]

    mdfiles = {}
    for v,cname in enumerate(rdum):
        #print cname
        base = os.path.basename(cname)
        dates = rdate[v]
        #print dates, etime,stime
            #print cname
        mval = match_dd(dates,ddates)
        #print( dates,ddates)
        
        if mval != 'no':
            dfile = ddum[mval[0]]
            print('Found DD match!', dfile)
            mdfiles[cname] = dfile
        else:
            mdfiles[cname] = None

        
    return mdfiles

def match_dd(rdate,ddates):
    dum=abs(rdate-np.array(ddates))

    try:
        mval=np.argwhere(dum == np.min(dum))[0]
        diff=np.min(dum)
        if diff.total_seconds() < 600.:
                return mval
        else:
                return 'no'
    except ValueError: 
        print ('Something DD is not working here!')


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
        print ('Something sound is not working here!')


def find_snd_match(config):
    rdum =[]
    with open(config['rfiles']) as f: 
        #dum.append(foo(f.readline()))
        #dum.append(foo(f.readline()))

        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    #print('sfiles:',config['sfiles'])
    #rdum = glob.glob(config['rfiles']+'*')
    #slist = sorted(glob.glob('{p}*{s}_*.txt'.format(p=config['sfiles'],s=(config['sstat']))))
    with open(config['sfiles']) as f:
        slist = f.read().splitlines()
    #slist = glob.glob(config['sfiles']+'*')
    sdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
#            print base
        radcdate=np.str(base[config['sdstart']:config['sdend']])
        #print('radcdate',radcdate)
        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config['sdate_format'])
        sdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
#            print cname
        base = os.path.basename(cname)
        radcdate=np.str(base[config['rdstart']:config['rdend']])
        #dates=datetime.datetime.strptime(radcdate,config['rdate_format'])
        dates = datetime.datetime.strptime('{r}'.format(r=radcdate),config['rdate_format'])
        sdt = datetime.datetime.strptime(config['sdatetime'],config['sdatetime_format'])
        edt = datetime.datetime.strptime(config['edatetime'],config['edatetime_format'])
        #if (dates >= config['etime']) and (dates <= config['stime']):
        if (dates <= edt) and (dates >= sdt):
            #print cname
            #now find a sounding match
            mv = match_snd(dates,sdates)
            if mv != 'no':
                #print ('sounding match',mv[0][0])
                msfiles[cname] = np.array(slist)[mv[0][0]]
            else:
                return None
        
    return msfiles

def foo(s1):
    return '{}'.format(s1.rstrip())

def reduce_dim(ds):
    try:
        t1= ds['time'][0].values
    except KeyError as ke:
        #print(f"{ke} skipping preprocessing")
        return(ds)
    for v in ds.data_vars.keys():
        try:
            ds[v]=ds[v].sel(time=t1).drop('time')
        except ValueError as e:
            pass
#            print(e)
#            print(v)
    return(ds)
    
from matplotlib.dates import DateFormatter,HourLocator
dayFormatter = DateFormatter('%H%M')      # e.g., 12
hourFormatter = DateFormatter('%H')      # e.g., 12


def hasNumbers(inputString):
     return any(char.isdigit() for char in inputString)


def polarris_driver(configfile):

    # =====
    # (1) Read in config file line by line.
    # =====

    config = {} # Load variable for config file data
    #print('ready to roll')

    print('\nReading '+str(configfile[0])+'...')
    with open(configfile[0]) as f:
        lines1 = [mm for mm in (line.replace('\t',' ') for line in f) if mm]
        lines2 = [nn for nn in (line.strip() for line in lines1) if nn] # NEW! Allow new lines in config file - can be skipped over!
        for line in lines2: #f:
            #print line
            if not line.startswith("#"):
                #print('line',line)
                key, val, comment = line.split('==')
                vval = val.replace(" ","")
                numck = hasNumbers(vval)
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'drop_vars' or key.replace(" ", "") == 'radarname' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'vr_name' or key.replace(" ", "") == 'mphys':
                    numck = False
                if key.replace(" ", "") == 'exper' or  key.replace(" ", "") == 'ptype':
                    vval = vval.strip("''")
                #print numck
                #print vval,key
                if key.replace(" ", "") == 'image_dir':
                    numck = True
                if key.replace(" ", "") == 'rfiles':
                    numck = True

                if numck is True or vval == 'None' or vval == 'True' or vval == 'False':
                    try:
                        config[(key.replace(" ", ""))] = eval(vval)
                    except:
                        if "datetime" in vval:
                            config[(key.replace(" ", ""))] = vval
                else:
                    config[(key.replace(" ", ""))] = vval
    
    time.sleep(3)
    print('Read-in complete.\n')

    # =====
    # (2) Find input radar files and concatenate the data. Rename x, y, z variables. 
    # =====

    #print('Finding and concatenating radar files in '+config['rfiles']+'...')
    drop_vars=config['drop_vars']
    with open(config['rfiles'], 'r') as f:
        #print(config['rfiles'])
        rfiles = f.read().splitlines()
    
    #print((config['exper']),(config['mphys']))
    print('Station/experiment: '+config['exper'])
    print('Input: '+config['mphys'])
    print('Start: '+config['sdatetime'])
    print('End: '+config['edatetime'])
    time.sleep(3)

    if config['exper'] == 'MC3E' and config['mphys'] == 'obs':
        print("special handling for ",config['exper'])

        file = open(config['rfiles'], "r")
        rf1=[]
        rf2=[]
        for line in file:
            print(line)
            if re.search('vtzms', line):
                rf1.append(line.rstrip('\n'))
            else:
                #print('other')
                rf2.append(line.rstrip('\n'))
        #print(rf1)
        rvar1 = xr.open_mfdataset(rf1,autoclose=True,concat_dim='d',preprocess=fix_my_data)
        rvar2=  xr.open_mfdataset(rf2,autoclose=True,concat_dim='d')
        rvar = xr.concat((rvar1,rvar2),dim='d')
        rfiles =list(np.append(rf1,rf2))
    else:
#        try:
#            print('trying to read normally')
#            rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim,combine='by_coords')
#        except ValueError as ve:
        #rfiles = glob.glob(config['rfiles']+"*")
        rvar = xr.open_mfdataset(rfiles,autoclose=True,combine='nested',concat_dim='d',preprocess=reduce_dim)
        #rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim)

    try:
        rvar = rvar.rename({'x0':'x'})
        rvar = rvar.rename({'y0':'y'})
        rvar = rvar.rename({'z0':'z'})
    except:
        print('Dims do not need renaming')
    #print('Current dimensions:',rvar.dims)

    if drop_vars == True:
        print("dropping extra variables for memory!")
        rvar= rvar.drop(['vrad03','vdop02','elev03','elev02','vdop03','vang02','vang03','vrad02','zhh02','zhh03','zdr02','zdr03','kdp02','kdp03','rhohv02','rhohv03'])
 
    print('Radar files ready.\n')
    time.sleep(3)

    # =====
    # (3) Get datetime objects from radar file names.
    # =====

    tm = []
    for d in rfiles:
        #print(d)
        dformat = config['rdate_format']
        base = os.path.basename(d)
        radcdate=np.str(base[config['rdstart']:config['rdend']])
        date=datetime.datetime.strptime(radcdate,dformat)
        tm.append(date)

    # =====
    # (4) 
    # =====

    if config['dd_on']==True:
        print('In your config file, dd_on is set to True.')
        time.sleep(3)
        with open(config['dfiles'], 'r') as g:
            dfiles1 = g.read().splitlines()
        #dfiles1 = glob.glob(config['dd_files']+"*")
        tmd = []
        for d in dfiles1:
            dformat = config['ddate_format']
            base = os.path.basename(d)
#            print('dd base',base,config['ddoff'],config['ddadd'])
            radcdate = base[config['ddstart']:config['ddend']]
#            print('dformat is',dformat,radcdate)
            if dformat == '%H%M':
            
                hr=int(base[config['ddstart']:config['ddstart']+2])
                mn=int(base[config['ddstart']+2:config['ddstart']+4])
                #print('hr','mn',hr,mn)
            #    print radcdate
                #date=datetime.datetime.strptime(radcdate,dformat)
                #print(config['date'])
                dstart=datetime.datetime.strptime(config['date'],'%Y%m%d')
                dat2 = datetime.datetime(dstart.year,dstart.month,dstart.day,hr,mn)
            else:
                dat2=datetime.datetime.strptime(radcdate,dformat)
                #dstart=datetime.datetime.strptime(config['date'],'%Y%m%d')
            tmd.append(dat2)

        print('Matching Dual-Doppler')
        dmatch = find_dd_match(rfiles,dfiles1,tm,tmd)
        #print('dmatch is ',dmatch)
        dvar = xr.open_mfdataset(dfiles1,concat_dim='d')

        # NEW! MultiDop names velocity fields in long-form. Shorten fieldnames in dopp files here for plotting labels.
        Uname = 'U'
        Vname = 'V'
        Wname = 'W'
        dvar = dvar.rename({config['uname']:Uname})
        dvar = dvar.rename({config['vname']:Vname})
        dvar = dvar.rename({config['wname']:Wname})

        wvar = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        wvar.fill(np.nan)

        unew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        unew.fill(np.nan)

        vnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        vnew.fill(np.nan)

        conv = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        conv.fill(np.nan)
        
        # NEW! MultiDop only works if distance values are in metres, not km. Need a condition to convert back to km so that doppler and radar distances are comparable.
        if np.array_equal(dvar.variables['x'].values, 1000.0*rvar.variables['x'].values):
            dvar['x'] = rvar['x']
            dvar['y'] = rvar['y']
            dvar['z'] = rvar['z']

        xsubmin = np.where(rvar.variables['x']==np.min(dvar.variables['x']))[0][0]
        xsubmax = np.where(rvar.variables['x']==np.max(dvar.variables['x']))[0][0]

        ysubmin = np.where(rvar.variables['y']==np.min(dvar.variables['y']))[0][0]
        ysubmax = np.where(rvar.variables['y']==np.max(dvar.variables['y']))[0][0]

        zsubmin = np.where(rvar.variables['z']==np.min(dvar.variables['z']))[0][0]
        zsubmax = np.where(rvar.variables['z']==np.max(dvar.variables['z']))[0][0]        

        for q,d in enumerate(dmatch.keys()):
            #print(q,'i outer',dmatch[d])
            if dmatch[d] is not None:
                #print('good, dmatch is not none')
                dfile = dmatch[d]
                if dfile in dfiles1:
                    #print(dfile)
                    i = dfiles1.index(dfile)
                    #print(i,'i inner')
                    #wvar[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['wname']].sel(d=i)
                    wvar[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Wname][i,:,:,:]
                    unew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Uname][i,:,:,:]
                    vnew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Vname][i,:,:,:]
                    #conv[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['convname']][i,:,:,:]

        rvar[Wname] = (['d','z','y','x'],wvar)
        rvar[Uname] = (['d','z','y','x'],unew)
        rvar[Vname] = (['d','z','y','x'],vnew)
        #rvar[config['convname']] = (['d','z','y','x'],conv)

    print('\nSending data to RadarData...')
   
    rdata = RadarData.RadarData(rvar,tm,ddata = None,dz=config['dz_name'],zdr=config['dr_name'],kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],u=Uname,v=Vname,w=Wname,conv=config['convname'],x=config['xname'],rr=config['rr_name'],band = config['band'],vr = config['vr_name'],lat_r=config['lat'],lon_r=config['lon'],y=config['yname'],z=config['zname'],lat=config['latname'], lon=config['lonname'],lat_0=config['lat'],lon_0=config['lon'],exper=config['exper'],mphys=config['mphys'],radar_name =config['radarname'],z_thresh=0,conv_types=config['conv_types'],strat_types=config['strat_types'],color_blind=config['cb_friendly'])

    if config['snd_on'] == True:
        print('In your config file, snd_on is set to True.')
        time.sleep(3)
        smatch = find_snd_match(config)
        #print("rfiles",rfiles[0])
        sfile = smatch[rfiles[0]]
        print('Matching Sounding')
    else:
        smatch = None
                                          
    if smatch is not None:
        print ('Found sounding match!',sfile,'\n')
        snd = SkewT.Sounding(sfile)
        rdata.add_sounding_object(snd) # this will add the sounding object to the radar object
                    # and then will take the heights and temps
        rdata.interp_sounding()

    if config['convert_Tk_Tc'] == True:
        print('converting T')
        rdata.convert_t()
    #print 'Calculating polarimetric fields like HID and rain...'
    #if config['pol_on'] == True:
    if config['mask_model'] == True:
        print('masking model data')
        rdata.mask_model()
        
    rdata.calc_pol_analysis()
#    print(config['cs_z'],'in 312 cs_z')
    rdata.calc_cs_shy(cs_z=config['cs_z'])
    rdata.raintype=rdata.data['CSS'].values#    rdata.set_hid()
    
    if config['comb_vicr'] == True:
        whvi = np.where(rdata.hid == 6)
        rdata.hid[whvi] = 3
   

    #Do some quick masking of the data####
#     mask = np.zeros([rdata.data.dims['d'],rdata.data.dims['z'],rdata.data.dims['y'],rdata.data.dims['x']])
#     whbad = np.logical_or(np.logical_or(np.logical_or(np.logical_or(rdata.data[rdata.dz_name].values>-20.,rdata.data[rdata.zdr_name].values>-2.),rdata.data[rdata.kdp_name].values<10.),rdata.data[rdata.zdr_name].values<10.),rdata.data[rdata.dz_name].values<70.)
#     whbad2= np.where(~whbad)
#     mask[whbad] = 1
#     if np.nanmin(rdata.data['CSS'].values)<1.:
#         mask[rdata.data['CSS'].values<=0] = 0
#     else:
#         mask[np.isnan(rdata.data['CSS'].values)]=0
# 
#     rdata.data['CSS'] = rdata.data['CSS'].where(mask ==1)
#     rdata.data[rdata.dz_name].values[whbad2] = np.nan
#     rdata.data[rdata.zdr_name].values[whbad2] = np.nan
#     rdata.data[rdata.kdp_name].values[whbad2] = np.nan
#     rdata.data[rdata.rho_name].values[whbad2] = np.nan
#     rdata.data[rdata.w_name].values[whbad2] = np.nan

    return rdata, config, Uname, Vname, Wname
