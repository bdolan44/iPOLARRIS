"""
This is the configruation file for setting up iPOLARRIS. Herein, the speicifcs of the dataset need to be defined such as the experiment, location and type of reading, etc.
Written by Brenda Dolan
May 2017
bdolan@atmos.colostate.edu
"""
from __future__ import print_function
import glob
import os
import sys
import re
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import RadarData
import datetime
import matplotlib.pyplot as plt
import numpy as np
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
            print ('Found DD match!', dfile)
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
    with open(config['radar_files']) as f: 
    #    dum.append(foo(f.readline()))
    #    dum.append(foo(f.readline()))

        for line in f:
            dat = (line)
            rdum.append(foo(dat))
    #print('sfiles:',config['sfiles'])
    slist = sorted(glob.glob('{p}*{s}_*.txt'.format(p=config['sfiles'],s=(config['sstat']))))
    sdates=[]
    for v,sname in enumerate(slist):

        base = os.path.basename(sname)
#            print base
        radcdate=np.str(base[13:13+10])
        #print('radcdate',radcdate)
        dates=datetime.datetime.strptime('{r}'.format(r=radcdate),config['sdate_format'])
        sdates.append(dates)

    msfiles = {}

    for v,cname in enumerate(rdum):
#            print cname
        base = os.path.basename(cname)
        radcdate=np.str(base[config['doff']:config['doff']+15])
        dates=datetime.datetime.strptime(radcdate,config['rdate_format'])
        if (dates >= config['etime']) and (dates <= config['stime']):
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

    config = {}
    print('ready to roll')
    with open(configfile[0]) as f:
        for line in f:
            #print line
            if not line.startswith("#"):
                #print('line',line)
                key, val, comment = line.split('==')
                vval = val.replace(" ","")
                numck = hasNumbers(vval)
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'dz_name' or key.replace(" ", "") == 'drop_vars' or key.replace(" ", "") == 'extrax' or key.replace(" ", "") == 'radarname' or key.replace(" ", "") == 'dr_name' or key.replace(" ", "") == 'kd_name' or key.replace(" ", "") == 'rh_name' or key.replace(" ", "") == 'vr_name' or key.replace(" ", "") == 'mphys':
                    numck = False
                if key.replace(" ", "") == 'exper' or key.replace(" ", "") == 'extra' or  key.replace(" ", "") == 'ptype' or key.replace(" ", "") == 'extrax':
                    vval = vval.strip("''")
                #print numck
                #print vval,key
                if key.replace(" ", "") == 'image_dir':
                    numck = True
                if key.replace(" ", "") == 'radar_files':
                    numck = True

                if numck is True or vval == 'None' or vval == 'True' or vval == 'False':
                    try:
                        config[(key.replace(" ", ""))] = eval(vval)
                    except:
                        if "datetime" in vval:
                            config[(key.replace(" ", ""))] = vval
                else:
                    config[(key.replace(" ", ""))] = vval
            
    print(config['radar_files'])
    drop_vars=config['drop_vars']
    with open(config['radar_files'], 'r') as f:
        rfiles = f.read().splitlines()
    #rfiles= glob.glob('*.nc')
    print((config['exper']),(config['mphys']))
    if config['exper'] == 'MC3E'  and config['mphys'] == 'obs':
        print("special handling for ",config['exper'])

        file = open(config['radar_files'], "r")
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
            print("trying nesting")
            rvar = xr.open_mfdataset(rfiles,autoclose=True,combine='nested',concat_dim='d',preprocess=reduce_dim)
    try:
        rvar = rvar.rename({'x0':'x'})
        rvar = rvar.rename({'y0':'y'})
        rvar = rvar.rename({'z0':'z'})
    except:
        print('Dims do not need renaming')
    print('Current dimensions:',rvar.dims)

    if drop_vars == True:
        print("dropping extra variables for memory!")
        rvar= rvar.drop(['vrad03','vdop02','elev03','elev02','vdop03','vang02','vang03','vrad02','zhh02','zhh03','zdr02','zdr03','kdp02','kdp03','rhohv02','rhohv03'])
    lon_0 = config['lon']
    lat_0 = config['lat']

    lat_r = config['lat']
    lon_r = config['lon']

    if config['snd_on'] == True:
        smatch = find_snd_match(config)
        #print("rfiles",rfiles[0])
        sfile = smatch[rfiles[0]]
        print('matching sounding')
    else:
        smatch = None



    tm = []
    for d in rfiles:
        print(d)
        dformat = config['wdate_format']
        base = os.path.basename(d)
        radcdate=np.str(base[config['time_parse'][0]:config['time_parse'][1]])
        date=datetime.datetime.strptime(radcdate,dformat)
        tm.append(date)

    if config['dd_on']==True:
        with open(config['dd_files'], 'r') as f:
            dfiles1 = f.read().splitlines()
        tmd = []
        for d in dfiles1:
            dformat = config['ddate_format']
            base = os.path.basename(d)
#            print('dd base',base,config['ddoff'],config['ddadd'])
            radcdate = base[config['ddoff']:config['ddadd']]
#            print (radcdate)
#            print('dformat is',dformat,radcdate)
            if dformat == '%H%M':
            
                hr=int(base[config['ddoff']:config['ddoff']+2])
                mn=int(base[config['ddoff']+2:config['ddoff']+4])
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

        wvar = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        wvar.fill(np.nan)

        unew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        unew.fill(np.nan)

        vnew = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        vnew.fill(np.nan)

        conv = np.zeros([rvar.dims['d'],rvar.dims['z'],rvar.dims['y'],rvar.dims['x']])
        conv.fill(np.nan)
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
                    i = dfiles1.index(dfile)
                    #print(i,'i inner')
                    wvar[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['wname']].sel(d=i)
                    unew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['uname']].sel(d=i)
                    vnew[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['vname']].sel(d=i)
                    conv[q,zsubmin:zsubmax+1,ysubmin:ysubmax+1,xsubmin:xsubmax+1] = dvar[config['convname']].sel(d=i)


        rvar[config['wname']] = (['d','z','y','x'],wvar)
        rvar[config['uname']] = (['d','z','y','x'],unew)
        rvar[config['vname']] = (['d','z','y','x'],vnew)
        rvar[config['convname']] = (['d','z','y','x'],conv)

    print('sending data to RadarData!')
    rdata = RadarData.RadarData(rvar,tm,ddata = None,dz =config['dz_name'],zdr=config['dr_name'],
                                                  kdp=config['kd_name'],rho=config['rh_name'],temp=config['t_name'],
                                                  u=config['uname'],v=config['vname'],w=config['wname'],conv=config['convname'],x=config['xname'],
                                                  rr=config['rr_name'],band = config['band'],vr = config['vr_name'],lat_r=lat_r,lon_r=lon_r,
                                                  y=config['yname'],z=config['zname'],lat=config['latname'], lon=config['lonname'],lat_0=lat_0,lon_0=lon_0,
                                                  exper=config['exper'],mphys=config['mphys'],radar_name =config['radarname'],
                                                  z_thresh=0,conv_types =  config['conv_types'],
                                                   strat_types = config['strat_types'])
                                               
    if smatch is not None:
        print ('Smatch',sfile)
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


    return rdata, config
