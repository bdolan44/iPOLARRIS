import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import numpy as np
#import RadarData
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from polarris_config import get_data
import RadarData
import GeneralFunctions as GF
from matplotlib import colors
plt.style.use('presentation')

def label_subplots(fig, xoff = 0.0, yoff = 0.02, nlabels = None,**kwargs):
    letters = ['a', 'b', 'c', 'd','e', 'f', 'g', 'h','l', 'm','n','o','p','q','r']
    figaxes = fig.get_axes()
    if nlabels is None: nlabels = len(figaxes)

    for fa in range(nlabels):
        xbox = figaxes[fa].get_position()
        xmin, ymax = xbox.xmin, xbox.ymax
    # this is the position I want
        if letters[fa] != '-':
            fig.text(xmin+xoff, ymax+yoff, '({})'.format(letters[fa]),**kwargs)




def plot_cfad_int(dat1,config,typ='dz',n1=None):
    fig, ax = plt.subplots(1,1,figsize=(8,6))
#    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}_{t}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,t=config['extra'])

    dat1cnt = np.shape(dat1['{t}cfad'.format(t=typ)])[0]

    cfad1_all = np.nansum(np.array(dat1['{t}cfad'.format(t=typ)]),axis=0)/dat1cnt
#    print dat1['hts'][0]

    if typ == 'w' or typ == 'wc' or typ == 'ws':
        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=ax,cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

    else:
        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ[0:2])],ax=ax,cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)
    ax.set_title(n1)

    ax.set_ylim(0,18)

    plt.tight_layout()

    plt.savefig('{id}CFAD_{tp}_{s}_int.{t}'.format(id=config['image_dir'],s=n1,t=config['ptype'],tp=typ.upper()),dpi=200)
    plt.clf()

def plot_hid_int(dat1,config,typ='hid',n1 = None):
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ht1sum = np.nansum(dat1['{t}cfad'.format(t=typ)], axis=0)
    dat1cnt = np.nanmax(ht1sum, axis=0) / 100.
    if n1 is None:
        n1 = '{e}_{x}_{t}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,t=config['extra'])

    fig, ax = GF.plot_hid_cdf(np.nansum(dat1['{t}cfad'.format(t=typ)], axis=0) / dat1cnt, dat1['hidhts'][0], ax=ax,
                              rconf=dat1['rconf'])

    #    print dat1cnt, dat1['hts']
#    fig, ax = plt.subplots(1,1,figsize=(12,8))
    fig, ax = GF.plot_hid_cdf(np.nansum(np.array(dat1['{t}cfad'.format(t=typ)]),axis=0)/dat1cnt,dat1['hidhts'][0],rconf=dat1['rconf'])
    ax.set_title(n1)

    ax.set_ylim(0,18)

    plt.tight_layout()
    
    plt.savefig('{id}CFAD_{h}_{s}_int.{t}'.format(id=config['image_dir'],h=typ.upper(),s=n1,t=config['ptype']),bbox_inches='tight', pad_inches=0.01,dpi=200)
    plt.clf()
    
def plot_hid_prof_int(dat1,config,typ='hid',n1 = None,n2 = None):
    fig, ax = plt.subplots(1,1,figsize=(12,8))
    if n1 is None:
        n1 = '{e}_{x}_{t}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,t=config['extra'])

    tw_water_vert1 = np.nansum(dat1['water_vert'],axis=0)
    tw_graup_vert1 = np.nansum(dat1['graup_vert'],axis=0)
    tw_hail_vert1 = np.nansum(dat1['hail_vert'],axis=0)
    tw_snow_vert1 = np.nansum(dat1['snow_vert'],axis=0)

    hts = dat1['hidhts'][0]

    lw=5
    ax.plot(tw_water_vert1,hts,color='blue',label='Water',lw=lw)
    ax.plot(tw_graup_vert1,hts,color='green',label='Graupel',lw=lw)
    ax.plot(tw_hail_vert1,hts,color='red',label='Hail',lw=lw)
    ax.plot(tw_snow_vert1,hts,color='goldenrod',label='Snow',lw=lw)
    ax.set_title('{e1} Hydromeor Freq.'.format(e1=n1),fontsize=20)
    ax.set_xlabel('Frequency',fontsize=18)
    ax.legend(loc='best',fontsize=18)
    ax.set_ylabel('Height (km)',fontsize=18)
    ax.set_ylim(0,20)
    plt.tight_layout()
    plt.savefig('{d}{e1}_hid_vert_int.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
    plt.clf()

def plot_joint_int(dat1,config,typ='zzdr',n1= None,n2=None):
    fig, ax = plt.subplots(1,1,figsize=(12,9))
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if typ == 'zzdr':              
#        print np.shape(dat1['edgzzdr'][0][0]),np.shape(np.nansum(dat1['histzzdr'],axis=0))
        cb6 = ax.contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],np.nansum(dat1['histzzdr'],axis=0))
        ax.set_xlabel('Zdr')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zzdr_int_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = ax.contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        ax.set_xlabel('Kdp')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zkdp_int_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = ax.contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        ax.set_xlabel('W (m/s)')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zw_int_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'wr':
#        print np.shape(dat1['edgwr'])
        cb6 = ax.contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],np.nansum(dat1['histwr'],axis=0).T)
        ax.set_ylabel('RR (mm/hr)')
        ax.set_xlabel('W (M/s)')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_wr_int_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

def plot_upwidth_int(dat1,config,n1= None):
        #print np.max(m_warea_wrf)
        plt.plot(np.nanmean(dat1['warea'],axis=0),dat1['wareat'][0],color='k',lw=5)
        plt.ylim(20,-60)
        plt.xlabel('Updraft Width (km$^2$)')
        plt.ylabel('Temperature (deg C)')
        plt.savefig('{d}{e1}_upwidth_int_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

def plot_uppercent_compare(dat1,dat2,config,n1= None,n2=None):
    fig, ax = plt.subplots(1,1,figsize=(8,8))
#    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    ax.plot(np.nanmean(dat1['w99'],axis=0),dat1['wperht'][0],lw=5,color='red',label='{e} 99th'.format(e=dat1['rconf'].exper))
    ax.plot(np.nanmean(dat1['w90'],axis=0),dat1['wperht'][0],lw=5,color='black',label='{e} 90th'.format(e=dat1['rconf'].exper))
    ax.plot(np.nanmean(dat1['w50'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='{e} 50th'.format(e=dat1['rconf'].exper))
    ax.set_xlabel('Vertical Wind Percentiles (m s$^{-1}$)')
    ax.set_ylabel('Height (km)')
    ax.set_title('{n} {n1}'.format(n=n1,n1=n2))

    ax.plot(np.nanmean(dat2['w99'],axis=0),dat2['wperht'][0],lw=5,color='red',label='{e} 99th'.format(e=dat2['rconf']),ls='--'.format(e=dat2['rconf'].exper))
    ax.plot(np.nanmean(dat2['w90'],axis=0),dat2['wperht'][0],lw=5,color='black',label='{e} 90th'.format(e=dat2['rconf']),ls='--'.format(e=dat2['rconf'].exper))
    ax.plot(np.nanmean(dat2['w50'],axis=0),dat2['wperht'][0],lw=5,color='gold',label='{e} 50th'.format(e=dat2['rconf']),ls='--'.format(e=dat2['rconf'].exper))
    # ax.set_xlabel('Vertical Wind Percentiles (m s$^{-1}$)')
    # ax.set_ylabel('Height (km)')
    # ax.set_title('{n}'.format(n=n2))
    ax.legend(loc='best')
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(16)

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    plt.savefig('{id}{e1}_{e2}_wpercent_compare_{s}_{x}.{t}'.format(id=config['image_dir'],e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,s=st_diff,x=config['extra'],t=config['ptype']),dpi=200)
    plt.clf()

def plot_uppercent_compare_updn(dat1, dat2, config, n1=None, n2=None):

    fig, ax = plt.subplots(1,3,figsize=(18,6))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    axf[0].plot(np.nanmean(dat1['w99u'],axis=0),dat1['wperht'][0],lw=5,color='red',label='99th',ls='-')
    axf[0].plot(np.nanmean(dat1['w90u'],axis=0),dat1['wperht'][0],lw=5,color='black',label='90th',ls='-')
    axf[0].plot(np.nanmean(dat1['w50u'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th',ls='-')
    axf[0].plot(np.nanmean(dat1['w99d'],axis=0),dat1['wperht'][0],lw=5,color='red',label='99th',ls='-')
    axf[0].plot(np.nanmean(dat1['w90d'],axis=0),dat1['wperht'][0],lw=5,color='black',label='90th',ls='-')
    axf[0].plot(np.nanmean(dat1['w50d'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th',ls='-')
    axf[0].legend(loc='best')

    axf[0].set_xlabel('Up/Down Percentiles (m s$^{-1}$)',fontsize=16)

    axf[0].set_ylabel('Height (km)',fontsize=16)
    axf[0].set_title('{n}'.format(n=n1),fontsize=16)

    axf[1].plot(np.nanmean(dat2['w99u'],axis=0),dat2['wperht'][0],lw=5,color='red',label='99th',ls='--')
    axf[1].plot(np.nanmean(dat2['w90u'],axis=0),dat2['wperht'][0],lw=5,color='black',label='90th',ls='--')
    axf[1].plot(np.nanmean(dat2['w50u'],axis=0),dat2['wperht'][0],lw=5,color='gold',label='50th',ls='--')
    axf[1].plot(np.nanmean(dat2['w99d'],axis=0),dat2['wperht'][0],lw=5,color='red',ls='--')
    axf[1].plot(np.nanmean(dat2['w90d'],axis=0),dat2['wperht'][0],lw=5,color='black',ls='--')
    axf[1].plot(np.nanmean(dat2['w50d'],axis=0),dat2['wperht'][0],lw=5,color='gold',ls='--')

    axf[1].set_xlabel('Up/Down Percentiles (m s$^{-1}$)',fontsize=16)
    axf[1].set_ylabel('Height (km)',fontsize=16)
    axf[1].set_title('{n}'.format(n=n2),fontsize=16)
    axf[1].legend(loc='best')

    axf[2].plot(np.nanmean(dat2['w99u'],axis=0),dat2['wperht'][0],lw=5,color='red',label='99th',ls='--')
    axf[2].plot(np.nanmean(dat2['w90u'],axis=0),dat2['wperht'][0],lw=5,color='black',label='90th',ls='--')
    axf[2].plot(np.nanmean(dat2['w50u'],axis=0),dat2['wperht'][0],lw=5,color='gold',label='50th',ls='--')
    axf[2].plot(np.nanmean(dat2['w99d'],axis=0),dat2['wperht'][0],lw=5,color='red',ls='--')
    axf[2].plot(np.nanmean(dat2['w90d'],axis=0),dat2['wperht'][0],lw=5,color='black',ls='--')
    axf[2].plot(np.nanmean(dat2['w50d'],axis=0),dat2['wperht'][0],lw=5,color='gold',ls='--')

    axf[2].plot(np.nanmean(dat1['w99u'],axis=0),dat1['wperht'][0],lw=5,color='red',label='99th',ls='-')
    axf[2].plot(np.nanmean(dat1['w90u'],axis=0),dat1['wperht'][0],lw=5,color='black',label='90th',ls='-')
    axf[2].plot(np.nanmean(dat1['w50u'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th',ls='-')
    axf[2].plot(np.nanmean(dat1['w99d'],axis=0),dat1['wperht'][0],lw=5,color='red',ls='-')
    axf[2].plot(np.nanmean(dat1['w90d'],axis=0),dat1['wperht'][0],lw=5,color='black',ls='-')
    axf[2].plot(np.nanmean(dat1['w50d'],axis=0),dat1['wperht'][0],lw=5,color='gold',ls='-')


    axf[2].set_xlabel('Up/Down Percentiles (m s$^{-1}$)',fontsize=16)
    axf[2].set_ylabel('Height (km)',fontsize=16)
    axf[2].set_title('{n} {n2}'.format(n=n1,n2=n2),fontsize=16)
#    axf[2].legend(loc='best')


    for a in axf:
        for label in (a.get_xticklabels() + a.get_yticklabels()):
            label.set_fontsize(16)

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    print st_diff
#    plt.savefig('test.png')
    plt.savefig('{id}{s}_updownstats_{x}.{t}'.format(id=config['image_dir'],s=st_diff,t=config['ptype'],x=config['extra']),dpi=200)
#
def plot_uppercent(dat1,config,n1= None):
    fig, ax = plt.subplots(1,2,figsize=(18,8))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

    axf[0].plot(np.nanmean(dat1['w99'],axis=0),dat1['wperht'][0],lw=5,color='black',label='99th')
    axf[0].plot(np.nanmean(dat1['w90'],axis=0),dat1['wperht'][0],lw=5,color='red',label='90th')
    axf[0].plot(np.nanmean(dat1['w50'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th')
    axf[0].set_xlabel('Vertical Wind Percentiles (m s$^{-1}$)',fontsize=16)
    axf[0].set_ylabel('Height (km)',fontsize=16)
    axf[0].set_title('{n}'.format(n=n1),fontsize=20)
    axf[0].legend(loc='best')
    axf[1].plot(np.nanmean(dat1['w99u'],axis=0),dat1['wperht'][0],lw=5,color='black',label='99th up')
    axf[1].plot(np.nanmean(dat1['w90u'],axis=0),dat1['wperht'][0],lw=5,color='red',label='90th up')
    axf[1].plot(np.nanmean(dat1['w50u'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th up')
    axf[1].plot(np.nanmean(dat1['w99d'],axis=0),dat1['wperht'][0],lw=5,color='black',label='99th down',ls='--')
    axf[1].plot(np.nanmean(dat1['w90d'],axis=0),dat1['wperht'][0],lw=5,color='red',label='90th down',ls='--')
    axf[1].plot(np.nanmean(dat1['w50d'],axis=0),dat1['wperht'][0],lw=5,color='gold',label='50th down',ls='--')
    axf[1].legend(loc='best')
    axf[1].set_xlabel('Up/Down Percentiles (m s$^{-1}$)',fontsize=16)

    axf[1].set_ylabel('Height (km)',fontsize=16)
    axf[1].set_title('{n}'.format(n=n1),fontsize=20)

    for a in axf:
        for label in (a.get_xticklabels() + a.get_yticklabels()):
            label.set_fontsize(16)


    plt.tight_layout()
    plt.savefig('{id}{e}_vvelstats_{x}.{t}'.format(id=config['image_dir'],e=dat1['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=200)
#    plt.clf()


def plot_upwidth(dat1,dat2,config,n1= None,n2=None):
    fig, ax = plt.subplots(1,2,figsize=(18,8))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    axf[0].plot(np.nanmean(dat1['warea'],axis=0),dat1['wareat'][0],color='k',lw=5)
    axf[0].set_ylim(20,-60)
    axf[0].set_xlabel('Updraft Width (km$^2$)')
    axf[0].set_ylabel('Temperature (deg C)')
    axf[0].set_title('{n}'.format(n=n1))

    axf[1].plot(np.nanmean(dat2['warea'],axis=0),dat2['wareat'][0],color='k',lw=5)
    axf[1].set_ylim(20,-60)
    axf[1].set_xlabel('Updraft Width (km$^2$)')
    axf[1].set_ylabel('Temperature (deg C)')
    axf[1].set_title('{n}'.format(n=n2))

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    plt.savefig('{id}{e1}_{e2}_wpwidth_compare_{s}_{x}.{t}'.format(id=config['image_dir'],e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,s=st_diff,x=config['extra'],t=config['ptype']),dpi=200)
    plt.clf()

def plot_joint_comp(dat1,dat2,config,typ='zzdr',n1= None,n2=None):
    fig, ax = plt.subplots(1,3,figsize=(18,8))
    axf = ax.flatten()

    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    if typ == 'zzdr':              
#        print np.shape(dat1['edgzzdr'][0][0]),np.shape(np.nansum(dat1['histzzdr'],axis=0))
        cb6 = axf[0].contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],np.nansum(dat1['histzzdr'],axis=0))
        axf[0].set_xlabel('Zdr')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])

        cb6 = axf[1].contourf(dat2['edgzzdr'][0][1][:-1],dat2['edgzzdr'][0][0][:-1],np.nansum(dat2['histzzdr'],axis=0))
        axf[1].set_xlabel('Zdr')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=config['extra'],m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histzzdr'],axis=0)-np.nansum(dat2['histzzdr'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')

        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('dBZ',fontsize=18)
        axf[2].set_xlabel('Zdr',fontsize = 18)
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.savefig('{d}{e1}_{e2}_zzdr_comp_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = axf[0].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        axf[0].set_xlabel('Kdp')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax[0])

        cb6 = axf[1].contourf(dat2['edgkdz'][0][1][:-1],dat2['edgkdz'][0][0][:-1],np.nansum(dat2['histkdz'],axis=0))
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=config['extra'],m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histkdz'],axis=0)-np.nansum(dat2['histkdz'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('Kdp')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zkdp_comp_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = axf[0].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        axf[0].set_xlabel('W (m/s)')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=config['extra'],m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])
        
        cb6 = axf[1].contourf(dat2['edgzw'][0][1][:-1],dat2['edgzw'][0][0][:-1],np.nansum(dat2['histzw'],axis=0))
        axf[1].set_xlabel('W (m/s)')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=config['extra'],m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])
        
        diffdat = np.nansum(dat1['histzw'],axis=0)-np.nansum(dat2['histzw'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('W (m/s)')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zw_comp_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'wr':
        cb6 = axf[0].contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],np.nansum(dat1['histwr'],axis=0).T)
        axf[0].set_ylabel('RR (mm/hr)')
        axf[0].set_xlabel('W (M/s)')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])
        
        cb6 = axf[1].contourf(dat2['edgwr'][0][0][:-1],dat2['edgwr'][0][1][:-1],np.nansum(dat2['histwr'],axis=0).T)
        axf[1].set_ylabel('RR (mm/hr)')
        axf[1].set_xlabel('W (M/s)')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=extra,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histwr'],axis=0)-np.nansum(dat2['histwr'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],diffdat.T,levels,cmap='bwr',extend='both')
        axf[2].set_ylabel('RR (mm/hr)')
        axf[2].set_xlabel('W (M/s)')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])
        
        plt.savefig('{d}{e1}_{e2}_wr_comp_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

def plot_cfad_compare(dat1,dat2,config,typ='dz',n1 = None,n2 = None,n3= None):
    fig, ax = plt.subplots(1,3,figsize=(18,8))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}_{t}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,t=config['extra'])
    if n2 is None:
        n2 = '{e}_{x}_{t}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys,t=config['extra'])
    if n3 is None:
        n3 = '{e}{m1}-{x}{m2}_{t}'.format(e=dat1['rconf'].exper,x=dat2['rconf'].exper,m1=dat1['rconf'].mphys,m2=dat2['rconf'].mphys,t=config['extra'])

    dat1cnt = np.shape(dat1['{t}cfad'.format(t=typ)])[0]
    dat2cnt = np.shape(dat2['{t}cfad'.format(t=typ)])[0]
#

    cfad1_all = np.sum(dat1['{t}cfad'.format(t=typ)],axis=0)/dat1cnt
    cfad2_all = np.sum(dat2['{t}cfad'.format(t=typ)],axis=0)/dat2cnt

#     print np.nanmax(cfad1_all)
#     print np.nanmax(cfad2_all)
    
    if typ == 'w':
        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=axf[0],cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad2_all, hts = dat2['hts'][0],  bins = dat2['{t}bins'.format(t=typ)],ax=axf[1],cfad_on = 0,rconf = dat2['rconf'],tspan = dat2['time'],maxval=20,cont=True,levels = True)

    else:
        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=axf[0],cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad2_all, hts = dat2['hts'][0],  bins = dat2['{t}bins'.format(t=typ)],ax=axf[1],cfad_on = 0,rconf = dat2['rconf'],tspan = dat2['time'],maxval=20,cont=True,levels = True)
    axf[0].set_title('{n}'.format(n=n1))
    axf[1].set_title('{n}'.format(n=n2))

    axf[0].set_ylim(0,18)
    axf[1].set_ylim(0,18)

    if len(dat1['hts'][0]) != len(dat2['hts'][0]):
        hvals = [dat1['hts'][0],dat2['hts'][0]]
        vals = np.array([cfad1_all, cfad2_all])

        lens = [len(dat1['hts'][0]),len(dat2['hts'][0])]
        sz = np.max(lens)
        arg = np.argmax(lens)
        cfad_new1=np.zeros_like(vals[arg])
        cfad_new2=np.zeros_like(vals[arg])
    #    print np.shape(vals[arg]),np.shape(cfad_new2)

        for i,h in enumerate(hvals[arg]):
            close_h = np.argmin(np.abs(h-hvals[0][:]))
            cfad_new1[i,:] = vals[0][close_h,:]
        
            close_h = np.argmin(np.abs(h-hvals[1][:]))
            cfad_new2[i,:] = vals[1][close_h,:]
        hts = hvals[arg]
#        print 'calc new CFADs', np.nanmax(cfad_new1), np.nanmax(cfad_new2),  dat1cnt, dat2cnt
        diff_cfad = cfad_new1-cfad_new2
    else:
        diff_cfad = cfad1_all - cfad2_all
        hts = dat1['hts'][0]

    
    cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
    maxa = np.nanpercentile(np.abs(cfad_ma),98)
    levels=np.linspace(-1*maxa,maxa,50)
#    print typ, maxa
    if typ=='w':
        maxa = np.around(np.nanpercentile(np.abs(cfad_ma), 96),decimals=1)
        nor = np.around(np.nanpercentile(np.abs(cfad_ma),93),decimals=1)
        delt = np.around((maxa+maxa)/50,decimals=2)
        print maxa,nor,delt
        levels = np.arange(-1 * maxa, maxa+delt, delt)
        cb=axf[2].contourf(dat1['{t}bins'.format(t=typ)][:-1],hts,cfad_ma,levels=levels,norm=colors.Normalize(vmin=-1*nor,vmax=nor),cmap='bwr',extend='both')

    else:

        cb=axf[2].contourf(dat1['{t}bins'.format(t=typ)][:-1],hts,cfad_ma,levels,cmap='bwr',extend='both')

    plt.colorbar(cb,ax=axf[2])
    axf[2].set_ylabel('Height (km MSL)',fontsize=18)

    if typ == 'drc' or typ == 'drs' or typ == 'dr':
        varn = 'DR'
    elif typ == 'dzc' or typ == 'dzs'  or typ == 'dz':
        varn = 'DZ'
    elif typ == 'kdc' or typ == 'kds'  or typ == 'kd':
        varn = 'KD'
    elif typ == 'wsvar' or typ == 'wsvar'  or typ == 'w':
        varn = 'Wvar'
    else:
        varn = typ
    try:

        axf[2].set_xlabel('{n} {u}'.format(n=dat1['rconf'].names[varn],u=dat1['rconf'].units[varn]),fontsize = 18)
    except:
    #     print 'Exception!'
        axf[2].set_xlabel('{tp}'.format(tp=typ.upper()),fontsize = 18)
    axf[2].set_title('{v}'.format(v=n3))

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)

    plt.savefig('{id}CFAD_{tp}_{s}_{x}.{t}'.format(id=config['image_dir'],s=st_diff,t=config['ptype'],x=config['extra'],tp=typ.upper()),dpi=200)
#    plt.clf()

def plot_hid_2panel(dat1,dat2,config,typ='hid',n1 = None,n2 = None,):
    dat1cnt = np.shape(dat1['hts'])[0]
    dat2cnt = np.shape(dat2['hts'])[0]
    if n1 is None:
        n1 = '{e}_{k}_{x}_{t}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,k=typ,t=config['extra'])
    if n2 is None:
        n2 = '{e}_{k}_{x}_{t}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys,k=typ,t=config['extra'])

    fig, ax = plt.subplots(1,2,figsize=(18,8))
    axf = ax.flatten()
    
    ht1sum = np.nansum(dat1['{t}cfad'.format(t=typ)],axis=0)
    dat1cnt = np.nanmax(ht1sum,axis=0)/100.

    ht2sum = np.nansum(dat2['{t}cfad'.format(t=typ)],axis=0)
    dat2cnt = np.nanmax(ht2sum,axis=0)/100.
    
    
    fig, ax = GF.plot_hid_cdf(np.nansum(dat1['{t}cfad'.format(t=typ)],axis=0)/dat1cnt,dat1['hidhts'][0],ax=axf[0],rconf=dat1['rconf'])
    axf[0].set_title(n1)
    fig, ax = GF.plot_hid_cdf(np.nansum(dat2['{t}cfad'.format(t=typ)],axis=0)/dat2cnt,dat2['hidhts'][0],ax=axf[1],rconf=dat2['rconf'])
    axf[1].set_title(n2)

    #        self.HID_barplot_colorbar(fig)  # call separate HID colorbar function for bar plots


    axf[0].set_ylim(0,18)
    axf[1].set_ylim(0,18)

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    plt.savefig('{id}CFAD_{h}_{s}_{x}.{t}'.format(id=config['image_dir'],h=typ.upper(),s=st_diff,x=config['extra'],t=config['ptype']),bbox_inches='tight',dpi=200)
    plt.clf()




def plot_hid_profile(dat1,dat2,config,typ='hid',n1 = None,n2 = None):
    fig, ax = plt.subplots(1,3,figsize=(18,9))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    tw_water_vert1 = np.nanmean(dat1['water_vert'],axis=0)
    tw_graup_vert1 = np.nanmean(dat1['graup_vert'],axis=0)
    tw_hail_vert1 = np.nanmean(dat1['hail_vert'],axis=0)
    tw_snow_vert1 = np.nanmean(dat1['snow_vert'],axis=0)

    tw_water_vert2 = np.nanmean(dat2['water_vert'],axis=0)
    tw_graup_vert2 = np.nanmean(dat2['graup_vert'],axis=0)
    tw_hail_vert2 = np.nanmean(dat2['hail_vert'],axis=0)
    tw_snow_vert2 = np.nanmean(dat2['snow_vert'],axis=0)
    hts1 = dat1['hidhts'][0]
    hts2 = dat2['hidhts'][0]

    lw=3
    axf[0].plot(tw_water_vert1,hts1,color='blue',label='Water',lw=lw)
    axf[0].plot(tw_graup_vert1,hts1,color='green',label='Graupel',lw=lw)
    axf[0].plot(tw_hail_vert1,hts1,color='red',label='Hail',lw=lw)
    axf[0].plot(tw_snow_vert1,hts1,color='goldenrod',label='Snow',lw=lw)
    axf[0].set_title('{e1} Hydromeor Freq.'.format(e1=n1),fontsize=20)
    axf[0].set_xlabel('Frequency',fontsize=18)
    axf[0].legend(loc='best',fontsize=18)
    axf[0].set_ylabel('Height (km)',fontsize=18)
    axf[0].set_ylim(0,20)

    axf[1].plot(tw_water_vert2,hts2,color='blue',label='Water',lw=lw)
    axf[1].plot(tw_graup_vert2,hts2,color='green',label='Graupel',lw=lw)
    axf[1].plot(tw_hail_vert2,hts2,color='red',label='Hail',lw=lw)
    axf[1].plot(tw_snow_vert2,hts2,color='goldenrod',label='Snow',lw=lw)
    axf[1].set_title('{e1} Hydromeor Freq.'.format(e1=n2),fontsize=20)
    axf[1].set_xlabel('Frequency',fontsize=18)
    axf[1].legend(loc='best',fontsize=18)
    axf[1].set_ylabel('Height (km)',fontsize=18)
    axf[1].set_ylim(0,20)

    if len(dat1['hts'][0]) != len(dat2['hts'][0]):
        hvals = [dat1['hts'][0],dat2['hts'][0]]
        wvals = np.array([tw_water_vert1,tw_water_vert2])
        gvals = np.array([tw_graup_vert1,tw_graup_vert2])
        hailvals = np.array([tw_hail_vert1,tw_hail_vert2])
        svals = np.array([tw_snow_vert1,tw_snow_vert2])

        lens = [len(dat1['hts'][0]),len(dat2['hts'][0])]
        sz = np.max(lens)
        arg = np.argmax(lens)
        wvals_new1=np.zeros_like(wvals[arg])
        wvals_new2=np.zeros_like(wvals[arg])

        gvals_new1=np.zeros_like(gvals[arg])
        gvals_new2=np.zeros_like(gvals[arg])

        hailvals_new1=np.zeros_like(hailvals[arg])
        hailvals_new2=np.zeros_like(hailvals[arg])

        svals_new1=np.zeros_like(svals[arg])
        svals_new2=np.zeros_like(svals[arg])

    #    print np.shape(vals[arg]),np.shape(cfad_new2)

        for i,h in enumerate(hvals[arg]):
            close_h = np.argmin(np.abs(h-hvals[0][:]))
            wvals_new1[i] = wvals[0][close_h]
        
            close_h = np.argmin(np.abs(h-hvals[1][:]))
            wvals_new2[i] = wvals[1][close_h]
            
            close_h = np.argmin(np.abs(h-hvals[0][:]))
            gvals_new1[i] = gvals[0][close_h]
        
            close_h = np.argmin(np.abs(h-hvals[1][:]))
            gvals_new2[i] = gvals[1][close_h]
            
            close_h = np.argmin(np.abs(h-hvals[0][:]))
            hailvals_new1[i] = hailvals[0][close_h]
        
            close_h = np.argmin(np.abs(h-hvals[1][:]))
            hailvals_new2[i] = hailvals[1][close_h]

            close_h = np.argmin(np.abs(h-hvals[0][:]))
            svals_new1[i] = svals[0][close_h]
        
            close_h = np.argmin(np.abs(h-hvals[1][:]))
            svals_new2[i] = svals[1][close_h]

        hts = hvals[arg]
        diff_water = wvals_new1-wvals_new2
        diff_graup = gvals_new1-gvals_new2
        diff_hail = hailvals_new1-hailvals_new2
        diff_snow = svals_new1-svals_new2
    else:
        diff_water = tw_water_vert1-tw_water_vert2
        diff_graup = tw_graup_vert1-tw_graup_vert2
        diff_hail = tw_hail_vert1-tw_hail_vert2
        diff_snow = tw_snow_vert1-tw_snow_vert2
        hts = dat1['hts'][0]



    axf[2].plot(diff_water,hts,color='blue',label='Water',lw=lw)
    axf[2].plot(diff_graup,hts,color='green',label='Graupel',lw=lw)
    axf[2].plot(diff_hail,hts,color='red',label='Hail',lw=lw)
    axf[2].plot(diff_snow,hts,color='goldenrod',label='Snow',lw=lw)
    axf[2].set_title('{e1}-{e2} Hydromeor Freq.'.format(e1=n1,e2=n2),fontsize=20)
    axf[2].set_xlabel('Frequency',fontsize=18)
    axf[2].legend(loc='best',fontsize=18)
    axf[2].set_ylabel('Height (km)',fontsize=18)
    axf[2].set_ylim(0,20)

    plt.savefig('{d}{e1}_{e2}_hid_vert_compare_{x}.{t}'.format(d=config['image_dir'],e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=config['extra'],t=config['ptype']),dpi=300)
    plt.clf() 



def plot_upstat(dat1,dat2,config,typ='hid',n1 = None,n2 = None):
    fig, ax = plt.subplots(1,3,figsize=(18,9))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad2_all, hts = dat2['hts'][0],  bins = dat2['{t}bins'.format(t=typ)],ax=axf[1],cfad_on = 0,rconf = dat2['rconf'],tspan = dat2['time'],maxval=20,cont=True,levels = True)
    axf[0].set_title(n1)
    axf[2].set_title(n2)

    axf[0].set_ylim(0,18)
    axf[1].set_ylim(0,18)

    diff_cfad = cfad1_all - cfad2_all
    cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)

    levels=np.arange(-2,2.1,0.1)
    cb=axf[2].contourf(dat1['{t}bins'.format(typ)][:-1],dat1['hts'][0],cfad_ma,levels,cmap='bwr',extend='both')

    plt.colorbar(cb,ax=axf[2])
    axf[2].set_ylabel('Height (km MSL)',fontsize=18)

        
def make_single_pplots(rdat,flags,config,y=None):
#    print 'in make_singl_pplots'
    tspan= [rdat.date,rdat.date]
    ts = tspan[0]
#    print ts
#    print rdat.exper
    te = tspan[1]
#    print ts, te
    title_string = '{e} {t} {d1:%Y%m%d-%H%M%S} {x}'.format(e=rdat.exper,t=rdat.mphys,d1=ts,x=config['extra'])
    xlim = config['xlim']
    ylim = config['ylim']
    y = config['y']
    z = config['z']

    if flags['cfad_mpanel_flag'] == True:
#        print 'Working on Cfad mpanel'
        fig, ax = plt.subplots(2,2,figsize=(18,12))
        axf = ax.flatten()

        if config['wname'] in rdat.data.variables.keys():
            rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
        rdat.cfad_plot(rdat.dz_name,ax = axf[1],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan)
        rdat.cfad_plot(rdat.zdr_name,ax= axf[2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan)
        rdat.cfad_plot(rdat.kdp_name,ax = axf[3],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
        plt.tight_layout()
#        print "{s:%Y%m%d%H%M%S}".format(s=ts[0])
        plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M%S}_{r}_{m}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,m=rdat.mphys,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

        if config['plot_cs'] == True:
            fig, ax = plt.subplots(2,2,figsize=(18,12))
            axf = ax.flatten()

            if config['wname'] in rdat.data.variables.keys():
                rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='convective')
            rdat.cfad_plot(rdat.dz_name,ax = axf[1],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='convective')
            rdat.cfad_plot(rdat.zdr_name,ax= axf[2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='convective')
            rdat.cfad_plot(rdat.kdp_name,ax = axf[3],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='convective')
            extrahold = config['extra']
            config['extra']='{e}_convective'.format(e=extrahold)
            plt.tight_layout()
    #        print "{s:%Y%m%d%H%M%S}".format(s=ts[0])
            plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M%S}_{r}_{m}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,m=rdat.mphys,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
            config['extra']=extrahold
            plt.clf()

            fig, ax = plt.subplots(2,2,figsize=(18,12))
            axf = ax.flatten()
            if config['wname'] in rdat.data.variables.keys():
                rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='stratiform')
            rdat.cfad_plot(rdat.dz_name,ax = axf[1],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='stratiform')
            rdat.cfad_plot(rdat.zdr_name,ax= axf[2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='stratiform')
            rdat.cfad_plot(rdat.kdp_name,ax = axf[3],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='stratiform')
            extrahold = config['extra']
            config['extra']='{e}_stratiform'.format(e=extrahold)
            plt.tight_layout()
    #        print "{s:%Y%m%d%H%M%S}".format(s=ts[0])
            plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M%S}_{r}_{m}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,m=rdat.mphys,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
            config['extra']=extrahold
            plt.clf()


        
    if flags['cfad_individ_flag'] == True:
        fig, ax = plt.subplots(1,1,figsize=(18,12))
#        axf = ax.flatten()
        if config['wname'] in rdat.data.variables.keys():

            rdat.cfad_plot(rdat.w_name,ax = ax,bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
            plt.tight_layout()
            plt.savefig('{d}{p}_CFAD_W_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
            plt.clf()
        

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.dz_name,ax = ax,bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_dBZ_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.zdr_name,ax= ax,bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_Zdr_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.kdp_name,ax = ax,bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
        plt.savefig('{d}{p}_CFAD_Kdp_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.tight_layout()
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.rho_name,ax = ax,z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_RHO_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()
        
    if flags['hid_cfad_flag'] == True:
        fig, ax = rdat.plot_hid_cdf()
        plt.savefig('{d}{p}_CFAD_HID_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        
        plt.clf()
        
    if flags['joint_flag'] == True:

        fig, ax = plt.subplots(2,2,figsize=(12,12))
        axf = ax.flatten()

        zzdr_wrf,ed = rdat.hist2d(varx=rdat.dz_name,vary=rdat.zdr_name,binsx=config['dzbins'],binsy=config['drbins'])
        rdat.plot_2dhist(zzdr_wrf,ed,ax=axf[0])
        axf[0].set_xlabel('Zdr')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title(title_string)
        zkdp_wrf,edk = rdat.hist2d(varx=rdat.dz_name,vary=rdat.kdp_name,binsx=config['dzbins'],binsy=config['kdbins'])
        rdat.plot_2dhist(zkdp_wrf,edk,ax=axf[1])
        axf[1].set_title(title_string)
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')


        zw_wrf,edw = rdat.hist2d(varx=rdat.dz_name,vary=rdat.w_name,binsx=config['dzbins'],binsy=config['wbins'])
        rdat.plot_2dhist(zw_wrf,edw,ax=axf[2])
        axf[2].set_title(title_string)
        axf[2].set_xlabel('W')
        axf[2].set_ylabel('dBZ')

        zr_wrf,edr = rdat.hist2d(varx='RRB',vary=rdat.w_name,binsx=config['rrbins'],binsy=config['wbins'],xthr=0.00000)
        cb6 = rdat.plot_2dhist(zr_wrf,edr,ax=axf[3],cbon=True)
        axf[3].set_title(title_string)
        axf[3].set_xlabel('W')
        axf[3].set_ylabel('RRB')
        axf[3].set_ylim(0,50)
        
        plt.savefig('{d}{p}_2dPDF_4panel_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()
        

    if flags['hid_prof'] == True:

        hts, mwrf_water_vert = rdat.hid_vertical_fraction(config['hidwater'],z_resolution =config['z_resolution'])
        hts, mwrf_graup_vert = rdat.hid_vertical_fraction(config['hidgraup'],z_resolution =config['z_resolution'])
        hts, mwrf_hail_vert = rdat.hid_vertical_fraction(config['hidhail'],z_resolution =config['z_resolution'])
        hts, mwrf_snow_vert = rdat.hid_vertical_fraction(config['hidsnow'],z_resolution =config['z_resolution'])
        lw = 4
        plt.plot(mwrf_water_vert,hts,color='b',label='water',lw=lw)
        plt.plot(mwrf_graup_vert,hts,color='g',label='graupel',lw=lw)
        plt.plot(mwrf_hail_vert,hts,color='r',label='hail',lw=lw)
        plt.plot(mwrf_snow_vert,hts,color = 'yellow',label='snow',lw=lw)
        plt.xlabel('Frequency (%)')
        plt.ylabel('Height (km)')
        plt.title(title_string)
        plt.legend(loc = 'best')
        plt.savefig('{d}{p}_HID_prof_{s:%Y%m%d%H%M%S}_{r}_{x}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype']),dpi=300)
        plt.clf()
        
    if flags['all_cappi']== True:
        #z=2.0
        #print xlim
        if config['cappi_multi'] == True:
            fig,ax = rdat.cappi_multiplot(ts=rdat.date,xlim=config['xlim'],ylim=config['ylim'],z=config['z'],res = config['cappi_vectres'],varlist = eval(config['cappi_vars']),vectors = eval(config['cvectors']),contours = eval(config['cappi_contours']))
    #        plt.tight_layout()
            label_subplots(fig,yoff=0.01,xoff=0.01,size=16,nlabels=6)
            plt.savefig('{d}{p}_polcappi_6panel_{s:%Y%m%d%H%M%S}_{r}_{x}_{z}km.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype'],z=config['z']),dpi=300)
            plt.clf()

        else:
            for v in eval(config['cappi_vars']):
                fig,ax = rdat.cappi_multiplot(ts=rdat.date,xlim=config['xlim'],ylim=config['ylim'],z=config['z'],res = config['cappi_vectres'],varlist = v,vectors = eval(config['cvectors']),contours = eval(config['cappi_contours']))
    #        plt.tight_layout()
                #label_subplots(fig,yoff=0.01,xoff=0.01,size=16,nlabels=1)
                plt.savefig('{d}{p}_polcappi_{v}_{s:%Y%m%d%H%M%S}_{r}_{x}_{z}km.{t}'.format(d=config['image_dir'],v=v,p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype'],z=config['z']),dpi=300)
                plt.clf()

        
    if flags['all_xsec']== True:
        #y=-12.5
        fig, ax = rdat.xsec_multiplot(ts=rdat.date,y=config['y'],vectors=eval(config['rvectors']),res = config['rhi_vectres'],xlim=config['xlim'],varlist=eval(config['rhi_vars']))
        
        #plt.tight_layout()
        nvars = len(eval(config['rhi_vars']))
        if nvars <=6:
            yof = 0.01
        else:
            yof=-0.02
        
        label_subplots(fig,yoff=yof,xoff=0.01,size=16,nlabels=nvars)
        plt.savefig('{d}{p}_polrhi_{v}panel_{s:%Y%m%d%H%M%S}_{r}_{x}_{y}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],v=nvars,t=config['ptype'],y=config['y']),dpi=300)
        plt.clf()


    if config['wname'] in rdat.data.variables.keys():
        if flags['up_width'] == True:
            tmp, m_warea_wrf = rdat.updraft_width_profile(thresh_dz=True)
            #print np.max(m_warea_wrf)
            plt.plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
            plt.ylim(20,-60)
            plt.xlabel('Updraft Width (km$^2$)')
            plt.ylabel('Temperature (deg C)')
            plt.title(title_string)
            plt.savefig('{d}{p}_upwidth_{s:%Y%m%d%H%M%S}_{r}_{x}_{y}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype'],y=config['y']),dpi=300)
            plt.clf()

    if flags['qr_cappi'] == True:
        rdat.cappi_multiplot(z=config['z'],ts=rdat.date,xlim=config['xlim'],ylim=config['ylim'],varlist=eval(config['mix_vars']))
        plt.savefig('{d}{p}_qcappi_6panel_{s:%Y%m%d%H%M%S}_{r}_{x}_{z}km.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype'],z=config['z']),dpi=300)
        plt.clf()
        
    if flags['qr_rhi'] == True:
        rdat.xsec_multiplot(ts=rdat.date,y=config['y'],xlim=config['xlim'],varlist=eval(config['mix_vars']))
        plt.savefig('{d}{p}_qrhi_6panel_{s:%Y%m%d%H%M%S}_{r}_{x}_{y}.{t}'.format(d=config['image_dir'],p=rdat.exper,s=rdat.date,r=rdat.radar_name,x=config['extra'],t=config['ptype'],y=config['y']),dpi=300)
        plt.clf()
        


