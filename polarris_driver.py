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
from copy import deepcopy
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
        radcdate=str(base[config['sdstart']:config['sdend']])
        #print('radcdate',radcdate)
        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config['sdate_format'])
        sdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
#            print cname
        base = os.path.basename(cname)
#        print('base',base)
        radcdate=str(base[config['rdstart']:config['rdend']])
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

def find_wrfpol_match(config):
    rdum =[]
    with open(config['rfiles']) as f: 
        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    
    with open(config['wfiles']) as f:
        slist = f.read().splitlines()
    
    wdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
        radcdate=str(base[config['wdstart']:config['wdend']])
        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config['wdate_format'])
        wdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
        base = os.path.basename(cname)
        radcdate=str(base[config['rdstart']:config['rdend']])
        dates = datetime.datetime.strptime('{r}'.format(r=radcdate),config['rdate_format'])
        sdt = datetime.datetime.strptime(config['sdatetime'],config['sdatetime_format'])
        edt = datetime.datetime.strptime(config['edatetime'],config['edatetime_format'])
        if (dates <= edt) and (dates >= sdt):
            mv = match_snd(dates,wdates)
            if mv != 'no':
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
        except KeyError as k:
        #except ValueError as e:
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

    print('\nReading '+str(configfile[0])+'...')
    with open(configfile[0]) as f:
        lines1 = [mm for mm in (line.replace('\t',' ') for line in f) if mm]
        lines2 = [nn for nn in (line.strip() for line in lines1) if nn] # NEW! Allow new lines in config file - can be skipped over!
        for line in lines2: #f:
            if not line.startswith("#"):
                key, val, comment = line.split('==')
                vval = val.replace(" ","")
                numck = hasNumbers(vval)
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'drop_vars' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'vr_name' or key.replace(" ", "") == 'mphys' or key.replace(" ", "") == 'xname' or key.replace(" ", "") == 'yname' or key.replace(" ", "") == 'zname' or key.replace(" ", "") == 'latname' or key.replace(" ", "") == 'lonname':
                    numck = False
                if key.replace(" ", "") == 'exper': # or key.replace(" ", "") == 'ptype':
                    vval = vval.strip("''")
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

    print('Read-in complete.\n')

    # =====
    # (2) Find input radar files and concatenate the data. Rename x, y, z variables. 
    # =====

    print('Station/experiment: '+config['exper'])
    print('Input: '+config['mphys'].upper())
    print('Start: '+config['sdatetime'])
    print('End: '+config['edatetime'])
    time.sleep(3)

    drop_vars = config['drop_vars']
    sdatetime = int(config['sdatetime'][0:8]+config['sdatetime'][9:13])
    edatetime = int(config['edatetime'][0:8]+config['edatetime'][9:13])
    rfiles = []
    with open(config['rfiles'], 'r') as f:
        allrfiles = f.read().splitlines()
        for rfile in allrfiles:
            fullname = os.path.basename(rfile)
            filedatestr = fullname[config['rdstart']:config['rdend']].replace('_','').replace(':','').replace('-','')
            filedate = int(filedatestr[0:-2])
            if filedate >= sdatetime and filedate <= edatetime:
                rfiles.append(rfile)
    
    if rfiles == []:
        print("\nOops! There is no radar data for the dates given in your config file. Exiting...\n")
        sys.exit(1)

    if config['exper'] == 'MC3E' and config['mphys'] == 'obs':
        print("special handling for ",config['exper'])

        file = open(config['rfiles'], "r")
        rf1=[]
        rf2=[]
        for line in file:
            if re.search('vtzms', line):
                rf1.append(line.rstrip('\n'))
            else:
                rf2.append(line.rstrip('\n'))
        if not rf2:
            rvar = xr.open_mfdataset(rf1,autoclose=True,combine='nested',concat_dim='d',preprocess=fix_my_data)
        else:
            rvar1 = xr.open_mfdataset(rf1,autoclose=True,combine='nested',compat='override',preprocess=fix_my_data)
            rvar2 = xr.open_mfdataset(rf2,autoclose=True,concat_dim='d')
            rvar = xr.concat((rvar1,rvar2),dim='d')
            rfiles = list(np.append(rf1,rf2))         
    else:
        rvar = xr.open_mfdataset(rfiles,autoclose=True,combine='nested',concat_dim='d',preprocess=reduce_dim)
        #rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim,combine='by_coords')
        #rvar = xr.open_mfdataset(rfiles,autoclose=True,concat_dim='d',preprocess=reduce_dim)

    if config['type'].startswith('wrf'):
        refvals = deepcopy(rvar[config['dz_name']].values)
        refvals[refvals < float(config['refthresh'])] = np.nan
        newref = xr.DataArray(refvals, dims=['d','z','y','x'], name=config['dz_name'])
        rvar[config['dz_name']] = newref
        
        parsers = ['dr_name','kd_name','rh_name','vr_name']
        for v in parsers:
            newvals = deepcopy(rvar[config[v]].values)
            newvals[newvals == -999.0] = np.nan
            newvals[np.isnan(refvals)] = np.nan
            newvar = xr.DataArray(newvals, dims=['d','z','y','x'], name=config[v])
            rvar[config[v]] = newvar
    
    # =====
    # (3) Get datetime objects from radar file names.
    # =====

    tm = []
    for d in rfiles:
        dformat = config['rdate_format']
        base = os.path.basename(d)
        radcdate=str(base[config['rdstart']:config['rdend']])
        date=datetime.datetime.strptime(radcdate,dformat)
        tm.append(date)

    rvar = rvar.rename({config['xname']:'x'})
    rvar = rvar.rename({config['yname']:'y'})
    
    if config['type'].startswith('obs'):
        rvar = rvar.rename({config['zname']:'z'})
    elif config['type'].startswith('wrf'):
        if 'd' in rvar['hgt'].dims:
            hgt = rvar['hgt'].values[0,:]
        else:
            hgt = rvar['hgt'].values
        newz = xr.DataArray(hgt, coords={'z': hgt})
        rvar['z'] = newz

    if drop_vars:
        print("dropping extra variables for memory!")
        rvar= rvar.drop(['vrad03','vdop02','elev03','elev02','vdop03','vang02','vang03','vrad02','zhh02','zhh03','zdr02','zdr03','kdp02','kdp03','rhohv02','rhohv03'])
 
    print('Radar files ready.')
    time.sleep(3)

    # =====
    # (4) 
    # =====

    if config['dd_on']:
        print('In your config file, dd_on is set to True.')
        time.sleep(3)
        with open(config['dfiles'], 'r') as g:
            dfiles1 = g.read().splitlines()
        tmd = []
        for d in dfiles1:
            dformat = config['ddate_format']
            base = os.path.basename(d)
            radcdate = base[config['ddstart']:config['ddend']]
            if dformat == '%H%M':
                hr=int(base[config['ddstart']:config['ddstart']+2])
                mn=int(base[config['ddstart']+2:config['ddstart']+4])
                dstart=datetime.datetime.strptime(config['date'],'%Y%m%d')
                dat2 = datetime.datetime(dstart.year,dstart.month,dstart.day,hr,mn)
            else:
                dat2=datetime.datetime.strptime(radcdate,dformat)
            tmd.append(dat2)

        print('Matching Dual-Doppler')
        dmatch = find_dd_match(rfiles,dfiles1,tm,tmd)
        try:
            dvar = xr.open_mfdataset(dfiles1,concat_dim='d')
        except ValueError as ve:
            print('Trying nested instead of concat_dim to read DD files')
            dvar = xr.open_mfdataset(dfiles1,combine='nested',concat_dim='d')
            nf= len(dfiles1)
        
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
            if dmatch[d] is not None:
                dfile = dmatch[d]
                if dfile in dfiles1:
                    i = dfiles1.index(dfile)
                    #wvar[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['wname']].sel(d=i)
                    wvar[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Wname][i,:,:,:]
                    unew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Uname][i,:,:,:]
                    vnew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[Vname][i,:,:,:]
                    #conv[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['convname']][i,:,:,:]

        rvar[Wname] = (['d','z','y','x'],wvar)
        rvar[Uname] = (['d','z','y','x'],unew)
        rvar[Vname] = (['d','z','y','x'],vnew)

    else:
        Uname = None
        Vname = None
        Wname = None

    print('\nSending data to RadarData...')
 
    if config['wrft_on']:
        print('In your config file, wrft_on is set to True.')
        time.sleep(3)
        if not 't_air' in list(rvar.keys()):
            wmatch = find_wrfpol_match(config)
            if len(wmatch) > 0:
                print('Found POLARRIS-f files!')
                try:
                    tvar = xr.open_mfdataset(list(wmatch.values()),concat_dim='d')
                except ValueError as ve:
                    tvar = xr.open_mfdataset(list(wmatch.values()),combine='nested',concat_dim='d')
                rvar[config['t_name']] = tvar['t_air']-273.15
        else:
            rvar[config['t_name']].values = deepcopy(rvar[config['t_name']])-273.15

    rdata = RadarData.RadarData(rvar,tm,ddata = None,dz=config['dz_name'],zdr=config['dr_name'],kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],u=Uname,v=Vname,w=Wname,conv=config['convname'],rr=config['rr_name'],band = config['band'],vr = config['vr_name'],lat_r=config['lat'],lon_r=config['lon'],lat=config['latname'], lon=config['lonname'],lat_0=config['lat'],lon_0=config['lon'],exper=config['exper'],mphys=config['mphys'],z_thresh=0,conv_types=config['conv_types'],strat_types=config['strat_types'],color_blind=config['cb_friendly'])

    if config['snd_on']:
        print('In your config file, snd_on is set to True.')
        time.sleep(3)
        smatch = find_snd_match(config)
        if len(smatch) > 0:
            sfile = smatch[rfiles[0]]
            print ('Found sounding match!',sfile,'\n')
            snd = SkewT.Sounding(sfile)
            rdata.add_sounding_object(snd) # this will add the sounding object to the radar object and then will take the heights and temps
            rdata.interp_sounding()
    
    if config['mask_model']:
        print('masking model data')
        rdata.mask_model()
    
    rdata.calc_pol_analysis(tm,config)
    rdata.calc_cs_shy(cs_z=config['cs_z'])
    rdata.raintype=rdata.data['CSS'].values
    
    if config['comb_vicr']:
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
