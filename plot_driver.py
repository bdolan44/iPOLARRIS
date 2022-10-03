from __future__ import print_function
import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import numpy as np
import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from copy import deepcopy
import RadarData
import GeneralFunctions as GF
from matplotlib import colors
#plt.style.use('./presentation.mplstyle')
#plt.style.use('default')

from matplotlib.dates import DateFormatter,HourLocator
dayFormatter = DateFormatter('%H%M')      # e.g., 12
hourFormatter = DateFormatter('%H')      # e.g., 12

from matplotlib.colors import from_levels_and_colors
import cartopy.crs as ccrs
import matplotlib.ticker as ticker


def label_subplots(fig, xoff = 0.0, yoff = 0.02, nlabels = None,**kwargs):
    letters = ['a', 'b', 'c', 'd','e', 'f', 'g', 'h','l', 'm','n','o','p','q','r']
    figaxes = fig.get_axes()
    if nlabels is None: nlabels = len(figaxes)

    for fa in range(nlabels):
        xbox = figaxes[fa].get_position()
        xmin, ymax = xbox.xmin, xbox.ymax
    # this is the position I want
        if letters[fa] != '-':
            fig.text(xmin+xoff, ymax+yoff, '({})'.format(letters[fa]),**kwargs) #,transform=figaxes[fa].transAxes)



def plot_cfad_int(dat1,config,typ='dz',n1=None):
    fig, ax = plt.subplots(1,1,figsize=(8,6))
#    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

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

    plt.savefig('{id}CFAD_{tp}_{s}_int.{t}'.format(id=outdir,s=n1,t=config['ptype'],tp=typ.upper()),dpi=200)
    plt.clf()

def plot_hid_int(dat1,config,typ='hid',n1 = None):
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ht1sum = np.nansum(dat1['{t}cfad'.format(t=typ)], axis=0)
    dat1cnt = np.nanmax(ht1sum, axis=0) / 100.
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

    fig, ax = GF.plot_hid_cdf(np.nansum(dat1['{t}cfad'.format(t=typ)], axis=0) / dat1cnt, dat1['hidhts'][0], ax=ax,
                              rconf=dat1['rconf'])

    #    print dat1cnt, dat1['hts']
#    fig, ax = plt.subplots(1,1,figsize=(12,8))
    fig, ax = GF.plot_hid_cdf(np.nansum(np.array(dat1['{t}cfad'.format(t=typ)]),axis=0)/dat1cnt,dat1['hidhts'][0],rconf=dat1['rconf'])
    ax.set_title(n1)

    ax.set_ylim(0,18)

    plt.tight_layout()
    
    plt.savefig('{id}CFAD_{h}_{s}_int.{t}'.format(id=outdir,h=typ.upper(),s=n1,t=config['ptype']),bbox_inches='tight', pad_inches=0.01,dpi=200)
    plt.clf()
    
def plot_hid_prof_int(dat1,config,typ='hid',n1 = None,n2 = None):
    fig, ax = plt.subplots(1,1,figsize=(12,8))
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

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
    plt.savefig('{d}{e1}_hid_vert_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
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
        ax.set_title('{e} {m}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zzdr_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = ax.contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        ax.set_xlabel('Kdp')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zkdp_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = ax.contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        ax.set_xlabel('W (m/s)')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zw_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'wr':
#        print np.shape(dat1['edgwr'])
        cb6 = ax.contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],np.nansum(dat1['histwr'],axis=0).T)
        ax.set_ylabel('RR (mm/hr)')
        ax.set_xlabel('W (M/s)')
        ax.set_title('{e} {m}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_wr_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

def plot_upwidth_int(dat1,config,n1= None):
        #print max(m_warea_wrf)
        plt.plot(np.nanmean(dat1['warea'],axis=0),dat1['wareat'][0],color='k',lw=5)
        plt.ylim(20,-60)
        plt.xlabel('Updraft Width (km$^2$)')
        plt.ylabel('Temperature (deg C)')
        plt.savefig('{d}{e1}_upwidth_int.{t}'.format(d=outdir,e1=dat1['rconf'].exper,t=config['ptype']),dpi=300)
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
    plt.savefig('{id}{e1}_{e2}_wpercent_compare_{s}.{t}'.format(id=outdir,e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,s=st_diff,t=config['ptype']),dpi=200)
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
    print (st_diff)
#    plt.savefig('test.png')
    plt.savefig('{id}{s}_updownstats.{t}'.format(id=outdir,s=st_diff,t=config['ptype']),dpi=200)
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
    plt.savefig('{id}{e}_vvelstats.{t}'.format(id=outdir,e=dat1['rconf'].exper,t=config['ptype']),dpi=200)
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
    plt.savefig('{id}{e1}_{e2}_wpwidth_compare_{s}.{t}'.format(id=outdir,e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,s=st_diff,t=config['ptype']),dpi=200)
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
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])

        cb6 = axf[1].contourf(dat2['edgzzdr'][0][1][:-1],dat2['edgzzdr'][0][0][:-1],np.nansum(dat2['histzzdr'],axis=0))
        axf[1].set_xlabel('Zdr')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histzzdr'],axis=0)-np.nansum(dat2['histzzdr'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')

        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('dBZ',fontsize=18)
        axf[2].set_xlabel('Zdr',fontsize = 18)
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.savefig('{d}{e1}_{e2}_zzdr_comp.{t}'.format(d=outdir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = axf[0].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        axf[0].set_xlabel('Kdp')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax[0])

        cb6 = axf[1].contourf(dat2['edgkdz'][0][1][:-1],dat2['edgkdz'][0][0][:-1],np.nansum(dat2['histkdz'],axis=0))
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histkdz'],axis=0)-np.nansum(dat2['histkdz'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('Kdp')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zkdp_comp.{t}'.format(d=outdir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = axf[0].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        axf[0].set_xlabel('W (m/s)')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])
        
        cb6 = axf[1].contourf(dat2['edgzw'][0][1][:-1],dat2['edgzw'][0][0][:-1],np.nansum(dat2['histzw'],axis=0))
        axf[1].set_xlabel('W (m/s)')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])
        
        diffdat = np.nansum(dat1['histzw'],axis=0)-np.nansum(dat2['histzw'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('W (m/s)')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zw_comp.{t}'.format(d=outdir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

    if typ == 'wr':
        cb6 = axf[0].contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],np.nansum(dat1['histwr'],axis=0).T)
        axf[0].set_ylabel('RR (mm/hr)')
        axf[0].set_xlabel('W (M/s)')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])
        
        cb6 = axf[1].contourf(dat2['edgwr'][0][0][:-1],dat2['edgwr'][0][1][:-1],np.nansum(dat2['histwr'],axis=0).T)
        axf[1].set_ylabel('RR (mm/hr)')
        axf[1].set_xlabel('W (M/s)')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histwr'],axis=0)-np.nansum(dat2['histwr'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],diffdat.T,levels,cmap='bwr',extend='both')
        axf[2].set_ylabel('RR (mm/hr)')
        axf[2].set_xlabel('W (M/s)')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])
        
        plt.savefig('{d}{e1}_{e2}_wr_comp.{t}'.format(d=outdir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,t=config['ptype']),dpi=300)
        plt.clf()

def plot_difference_cfad(rdat1,rdat2,var1,var2,config1,bins=np.arange(0,82,2),xlab=None,savefig=True,n1=None,n2=None,n3=None,cscfad=None, nor=False,ylim=None,xlim=None):
    ###Pass a value for nor in order to normalize the colorbar over a standard range rather than normalization across values within the data. To normalize within the data, pass a value of False
    r1cdf,r1bins,r1ht = rdat1.cfad(var1,ret_z=1,z_resolution=1.0,value_bins=bins,cscfad=cscfad)
    r2cdf,r2bins,r2ht = rdat2.cfad(var2,ret_z=1,z_resolution=1.0,value_bins=bins,cscfad=cscfad)
    print('In plot_driver, csfad is ',cscfad)
    if n1 is None:
        n1 = rdat1.exper
    if n2 is None:
        n2 = rdat2.exper
    if n3 is None:
        n3 = '{a}-{b}'.format(a=rdat1.exper,b=rdat2.exper)

    fig, axf = plot_cfad_compare(r1cdf,r2cdf,r1ht,r2ht,r1bins,r2bins,config1,xlab=xlab,n1=n1,n2=n2,n3=n3,typ='dz',nor=nor,ylim=ylim,xlim=xlim)
    #if cscfad is not False:
    #    plt.suptitle('{c} {l}'.format(c=cscfad,l=lonvar),y=1.05,fontsize=30)
    #else:
    #    plt.suptitle('{l}'.format(l=lonvar),y=1.05,fontsize=30)
    
    #axf[0].set_xlabel('{l} bin'.format(l=lonvar))
    #axf[1].set_xlabel('{l} bin'.format(l=lonvar))
    #axf[2].set_xlabel('{l} bin'.format(l=lonvar))
    #if savefig == True:
    #    if cscfad is not False:
    #        plt.savefig('{d}CFAD_diff_{e1}_{e2}_{c}{l}.{p}'.format(p=config['ptype'],d=config1['image_dir'],c=cscfad,e1=rdat1.exper,e2=rdat2.exper,l=var1),dpi=400,bbox_inches='tight')
    #    else:
    #        plt.savefig('{d}CFAD_diff_{e1}_{e2}_{l}.{p}'.format(p=config['ptype'],d=config1['image_dir'],e1=rdat1.exper,e2=rdat2.exper,l=var1),dpi=400,bbox_inches='tight')
    #    return fig, axf
    #else:
    #    
    return fig,axf

def plot_cfad_compare(dat1,dat2,ht1,ht2,bin1,bin2,config,typ='dz',xlab = None, n1 = None,n2 = None,n3= None,savefig=False,nor=False,ylim=False,xlim=False):
    ###Pass a value for nor in order to normalize the colorbar over a standard range rather than normalization across values within the data. To normalize within the data, pass a value of False
    fig, ax = plt.subplots(1,4,figsize=(14,8),gridspec_kw={'wspace': 0.1,'hspace': 0.05,'width_ratios': [4,4,1.2,4],\
        'top':1., 'bottom':0., 'left':0., 'right':1.})
    axf = ax.flatten()

    #dat1cnt = np.shape(dat1)[0]
    #dat2cnt = np.shape(dat2)[0]

    #cfad1_all = np.sum(dat1,axis=0)/dat1cnt
    #cfad2_all = np.sum(dat2,axis=0)/dat2cnt
    
    cfad1_all = dat1
    cfad2_all = dat2


#     print np.nanmax(cfad1_all)
#     print np.nanmax(cfad2_all)
    
    if typ == 'w':
        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),data = cfad1_all, hts = ht1,  bins = bin1,ax=axf[0],cfad_on = 0,rconf =config,tspan = dat1['time'],maxval=20,cont=True,levels = True,cbyes=0, xlim=xlim, ylim=ylim, xlab=xlab)

        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad2_all, hts =ht2,  bins = bin2,ax=axf[1],cfad_on = 0,rconf = confing,tspan = dat2['time'],maxval=20,cont=True,levels = True,cbyes=1, xlim=xlim, ylim=ylim, xlab=xlab)

    else:
        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad1_all, hts = ht1,  bins = bin1,ax=axf[0],cfad_on = 0,rconf = config,tspan = config['sdatetime']+'_'+config['edatetime'],maxval=20,cont=True,levels = True,cbyes=0, xlim=xlim, ylim=ylim, xlab=xlab)

        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad2_all, hts = ht2,  bins = bin2,ax=axf[1],cfad_on = 0,rconf = config,tspan = config['sdatetime']+'_'+config['edatetime'],maxval=20,cont=True,levels = True,cbyes=1, xlim=xlim, ylim=ylim, xlab=xlab)
    
    #axf[0].set_title('{n}'.format(n=n1))
    #axf[1].set_title('{n}'.format(n=n2))

    #axf[0].set_ylim(0,18)
    #axf[1].set_ylim(0,18)
    #axf[2].set_ylim(0,18)

    if len(ht1) != len(ht2):
        print('fixing heights')
        hvals = [ht1,ht2]
        vals = np.array([cfad1_all, cfad2_all])

        lens = [len(ht1),len(ht2)]
        sz = max(lens)
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
        hts = ht1
    
    cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
    maxa = np.nanpercentile(np.abs(cfad_ma),98)
    levels=np.linspace(-1*maxa,maxa,50)
    
    axf[2].remove()

    if nor is False:
    #    print typ, maxa
        if typ=='w':
            maxa = np.around(np.nanpercentile(np.abs(cfad_ma), 96),decimals=1)
            nor = np.around(np.nanpercentile(np.abs(cfad_ma),93),decimals=1)
            delt = np.around((maxa+maxa)/50,decimals=2)
            print( maxa,nor,delt)
            levels = np.arange(-1 * maxa, maxa+delt, delt)
            cb=axf[3].contourf(bin1[:-1],hts,cfad_ma,levels=levels,norm=colors.Normalize(vmin=-1*nor,vmax=nor),cmap='bwr',extend='both')

        else:
            nor = np.around(np.nanpercentile(np.abs(cfad_ma),98),decimals=1)
            cb=axf[3].contourf(bin1[:-1],hts,cfad_ma,levels,cmap='bwr',norm=colors.Normalize(vmin=-1.*nor,vmax=nor),extend='both')
    else:
        nor=nor
        cb=axf[3].contourf(bin1[:-1],hts,cfad_ma,levels,cmap='bwr',norm=colors.Normalize(vmin=-1.*nor,vmax=nor),extend='both')
    
    axf[3].set_xlim(xlim)
    axf[3].set_ylim(ylim)
    axf[3].set_xlabel(xlab,fontsize=16)
    axf[3].set_yticks([])
    axf[3].set_yticklabels([])
    axf[3].tick_params(axis='x', which='major', labelsize=16)
    axf[3].tick_params(axis='y', which='major', labelsize=0)

    lur,bur,wur,hur = axf[3].get_position().bounds
    cbar_ax_dims = [lur+wur+0.02,bur,0.02,hur]
    cbar_ax = fig.add_axes(cbar_ax_dims)
    cbt = plt.colorbar(cb,cax=cbar_ax)
    cbt.set_ticks(np.arange(-1*nor,nor+1,0.5*nor))
    cbt.ax.tick_params(labelsize=16)
    cbt.set_label('Relative Difference (%)', fontsize=16, rotation=270, labelpad=15)

    #cb3 = plt.colorbar(cb,ax=axf[2])
    #cb3.set_label('Relative difference (%)',fontsize=16,rotation=270,labelpad=20)
    #cb3.ax.tick_params(labelsize=16)
    #cb3.set_ticks(np.linspace(-1.*nor,nor,9))
#    print('nor',nor)
    #axf[2].set_ylabel('Height (km MSL)',fontsize=18)

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
    #try:

    #    axf[2].set_xlabel('{n} {u}'.format(n=dat1['rconf'].names[varn],u=dat1['rconf'].units[varn]),fontsize = 18)
    #except:
    #     print 'Exception!'
    #axf[2].set_xlabel('{tp}'.format(tp=typ.upper()),fontsize = 18)
    #axf[2].set_title('{v}'.format(v=n3))


    #plt.tight_layout()
    #if savefig==True:
    #
    #    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    #
    #    plt.savefig('{id}CFAD_{tp}_{s}.{t}'.format(id=outdir,s=st_diff,t=config['ptype'],tp=typ.upper()),dpi=200)
    #else:
       
    return fig, axf
#    plt.clf()

def plot_hid_2panel(dat1,dat2,config,typ='hid',n1 = None,n2 = None,):
    dat1cnt = np.shape(dat1['hts'])[0]
    dat2cnt = np.shape(dat2['hts'])[0]
    if n1 is None:
        n1 = '{e}_{k}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys,k=typ)
    if n2 is None:
        n2 = '{e}_{k}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys,k=typ)

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
    plt.savefig('{id}CFAD_{h}_{s}.{t}'.format(id=outdir,h=typ.upper(),s=st_diff,t=config['ptype']),bbox_inches='tight',dpi=200)
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
        sz = max(lens)
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

    plt.savefig('{d}{e1}_{e2}_hid_vert_compare.{t}'.format(d=outdir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,t=config['ptype']),dpi=300)
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

        
def make_single_pplots(rdat,config,y=None):
    
    tspan= [rdat.date[0],rdat.date[-1]]
    tms = np.array(rdat.date)
    outpath = config['image_dir']

    if (config['pol_compare'] | config['all3']):
     
        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel CFADs for various polarimetric vars.\n')
    
        outdir = outpath+'pol_compare/'
        os.makedirs(outdir,exist_ok=True)

        if config['wname'] in rdat.data.variables.keys():
            numr,numc = 2,2
            figsize=(16,12)
        else:
            numr,numc = 1,3
            figsize=(16,8)
        
        fig, ax = plt.subplots(numr,numc,figsize=figsize,gridspec_kw={'wspace': 0.05, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        axf = ax.flatten()

        if config['wname'] in rdat.data.variables.keys():
            dum =rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,ylab=True)
            print('Panel 1: '+rdat.w_name)

        dum =rdat.cfad_plot(rdat.dz_name,ax = axf[numr*numc-3],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,ylab=True if numr==1 else False)
        print('Panel '+str(numr*numc-2)+': '+rdat.dz_name)

        dum =rdat.cfad_plot(rdat.zdr_name,ax= axf[numr*numc-2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,ylab=True if numr==2 else False)
        print('Panel '+str(numr*numc-1)+': '+rdat.zdr_name)

        dum =rdat.cfad_plot(rdat.kdp_name,ax = axf[numr*numc-1],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan)
        print('Panel '+str(numr*numc)+': '+rdat.kdp_name)

        lur1,bur1,wur1,hur1 = axf[1].get_position().bounds
        lur2,bur2,wur2,hur2 = axf[-1].get_position().bounds
        cbar_ax_dims = [lur2+wur2+0.02,bur2-0.001,0.03,bur1+hur1]
        cbar_ax = fig.add_axes(cbar_ax_dims)
        cbt = plt.colorbar(dum[-2],cax=cbar_ax)
        cbt.ax.tick_params(labelsize=20)
        cbt.set_label('Frequency (%)', fontsize=22, rotation=270, labelpad=20)
        cbt.set_ticks(dum[-1])
        cbt.set_ticklabels(dum[-1])

        axf[0].text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=24, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left
        
        plt.savefig('{d}{p}_CFAD_4panel.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype']),dpi=400,bbox_inches='tight')
        plt.clf()

        print('\nDone! Saved to '+outdir)
        print('\nIN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel CONVECTIVE CFADs for various polarimetric vars.\n')
        
        if config['wname'] in rdat.data.variables.keys():
            numr,numc = 2,2
            figsize=(16,12)
        else:
            numr,numc = 1,3
            figsize=(16,8)
        
        fig, ax = plt.subplots(2,2,figsize=(18,12),constrained_layout=True)
        axf = ax.flatten()

        if config['wname'] in rdat.data.variables.keys():
            dum =rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='convective',cont=True)
            print('Panel 1: '+rdat.w_name)
        
        dum =rdat.cfad_plot(rdat.dz_name,ax = axf[numr*numc-3],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='convective',cont=True)
        print('Panel '+str(numr*numc-2)+': '+rdat.dz_name)
        
        dum =rdat.cfad_plot(rdat.zdr_name,ax = axf[numr*numc-2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='convective',cont=True)
        print('Panel '+str(numr*numc-1)+': '+rdat.zdr_name)
        
        dum =rdat.cfad_plot(rdat.kdp_name,ax = axf[numr*numc-1],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='convective',cont=True)
        print('Panel '+str(numr*numc)+': '+rdat.kdp_name)
        
        extrax ='conv'
        #plt.tight_layout()
#        print "{s:%Y%m%d%H%M%S}".format(s=ts[0])
        plt.savefig('{d}{p}_CFAD_4panel_{x}.{t}'.format(d=outdir,p=rdat.exper,x=extrax,t=config['ptype']),dpi=400,bbox_inches='tight')
        plt.clf()
        
        print('\nDone! Saved to '+outdir)
        print('\nIN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel STRATIFORM CFADs for various polarimetric vars.\n')

        fig, ax = plt.subplots(2,2,figsize=(18,12))
        axf = ax.flatten()

        if config['wname'] in rdat.data.variables.keys():
            dum =rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=config['wbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='stratiform')
            print('Panel 1: '+rdat.w_name)
        
        dum =rdat.cfad_plot(rdat.dz_name,ax = axf[numr*numc-3],bins=config['dzbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='stratiform')
        print('Panel '+str(numr*numc-2)+': '+rdat.dz_name)
        
        dum =rdat.cfad_plot(rdat.zdr_name,ax= axf[numr*numc-2],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan= tspan,cscfad='stratiform')
        print('Panel '+str(numr*numc-1)+': '+rdat.zdr_name)
        
        dum =rdat.cfad_plot(rdat.kdp_name,ax = axf[numr*numc-1],bins=config['drbins'],z_resolution=config['z_resolution'],levels='levs',tspan = tspan,cscfad='stratiform')
        print('Panel '+str(numr*numc)+': '+rdat.kdp_name)
        
        extrax ='strat'
        #plt.tight_layout()
#        print "{s:%Y%m%d%H%M%S}".format(s=ts[0])
        #plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M%S}_{r}_{m}_{x}.{t}'.format(d=outdir,p=rdat.exper,s=tstart,m=rdat.mphys,r=rdat.band+'-band',t=config['ptype'],x=extrax),dpi=300) 
        plt.savefig('{d}{p}_CFAD_4panel_{x}.{t}'.format(d=outdir,p=rdat.exper,x=extrax,t=config['ptype']),dpi=400,bbox_inches='tight')
        plt.clf()

        print('Done! Saved to '+outdir)
        print('Moving on.\n')
 

    if (config['cfad_multi'] | config['all3']):
 
        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel CFADs for various polarimetric vars.\n')
        
        outdir = outpath+'cfad_multi/'
        os.makedirs(outdir,exist_ok=True)
 
        zmax = config['zmax']
        st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
        en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')
        if st.startswith(en): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' UTC'
        else:
            if st[0:8].startswith(en[0:8]): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+'-'+en[9:11]+':'+en[11:13]+' UTC'
            else: dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' - '+en[0:4]+'-'+en[4:6]+'-'+en[6:8]+' '+en[9:11]+':'+en[11:13]+' UTC'

        if not zmax == '':
            fig, ax = rdat.cfad_multiplot(varlist = eval(config['cfad_vars']),z_resolution=config['z_resolution'],zmax=zmax)
        else:
            fig, ax = rdat.cfad_multiplot(varlist = eval(config['cfad_vars']),z_resolution=config['z_resolution'])
        
        nvars=0
        for var in eval(config['cfad_vars']):
            if var in rdat.data.variables.keys():
                nvars+=1
        
        if nvars <=6:
            yof = 0.01
        else:
            yof = -0.02
        yof = -0.01
        xof = 0.01
        
        label_subplots(fig,yoff=yof,xoff=xof,size=16,nlabels=nvars,horizontalalignment='left',verticalalignment='top',color='k',bbox=dict(facecolor='w', edgecolor='w', pad=2.0),weight='bold')

        axf = ax.flatten()
        axf[0].text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left
        axf[2].text(1, 1, dtlab, horizontalalignment='right', verticalalignment='bottom', size=20, color='k', zorder=10, weight='bold', transform=axf[2].transAxes) # (a) Top-left
 
        if config['ptype'].startswith('mp4'):
            plt.savefig('{d}{p}_CFAD_{t1}-{t2}.png'.format(d=outdir,p=rdat.exper,t1=st,t2=en),dpi=400,bbox_inches='tight')
        else: 
            plt.savefig('{d}{p}_CFAD_{t1}-{t2}.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype'],t1=st,t2=en),dpi=400,bbox_inches='tight')

        plt.close()

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')


    if (config['cfad_individ'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating individual CFADs for various polarimetric vars.\n')
        outdir = outpath+'cfad_individ/'
        os.makedirs(outdir,exist_ok=True)

        allvars = eval(config['cfad_vars'])

        zmax = config['zmax']
        st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
        en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')

        if st.startswith(en): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' UTC'
        else:
            if st[0:8].startswith(en[0:8]): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+'-'+en[9:11]+':'+en[11:13]+' UTC'
            else: dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' - '+en[0:4]+'-'+en[4:6]+'-'+en[6:8]+' '+en[9:11]+':'+en[11:13]+' UTC'

        for i,v in enumerate(eval(config['cfad_vars'])):
            
            if v is None:
                continue 
            else:
                
                if v.startswith('HID'):
 
                    print(v)

                    if not zmax == '':
                        fig, ax = rdat.plot_hid_cdf(cbar=2,z_resolution=config['z_resolution'],zmax=zmax)
                    else:
                        fig, ax = rdat.plot_hid_cdf(cbar=2,z_resolution=config['z_resolution'])

                    ax.text(0,1,'{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=ax.transAxes)
                    ax.text(1,1, dtlab, horizontalalignment='right', verticalalignment='bottom', size=18, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                   
                    if config['ptype'].startswith('mp4'):
                        plt.savefig('{d}{p}_HID_CFAD_{t1}-{t2}.png'.format(d=outdir,p=rdat.exper,t1=st,t2=en),dpi=400,bbox_inches='tight')
                    else: 
                        plt.savefig('{d}{p}_HID_CFAD_{t1}-{t2}.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype'],t1=st,t2=en),dpi=400,bbox_inches='tight')

                    plt.close()
                
                else:

                    if not rdat.cfbins[v] is '': # and v in rdat.data.variables.keys():
                        
                        print(v)

                        if not zmax == '':
                            cfad, hts, pl, fig, ax = rdat.cfad_plot(v,bins=rdat.cfbins[v],z_resolution=config['z_resolution'],levels=1,zmax=zmax)
                        else:
                            cfad, hts, pl, fig, ax = rdat.cfad_plot(v,bins=rdat.cfbins[v],z_resolution=config['z_resolution'],levels=1)
                            
                        ax.text(0,1,'{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=ax.transAxes)
                        ax.text(1,1, dtlab, horizontalalignment='right', verticalalignment='bottom', size=18, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                      
                        if config['ptype'].startswith('mp4'):
                            plt.savefig('{d}{p}_{v}_CFAD_{t1}-{t2}.png'.format(d=outdir,p=rdat.exper,v=rdat.names_uc[v],t1=st,t2=en),dpi=400,bbox_inches='tight')
                        else:
                            plt.savefig('{d}{p}_{v}_CFAD_{t1}-{t2}.{t}'.format(d=outdir,p=rdat.exper,v=rdat.names_uc[v],t=config['ptype'],t1=st,t2=en),dpi=400,bbox_inches='tight')

                        plt.close()

                    else:
                        
                        continue

                print('\nDone! Saved to '+outdir)
                print('Moving on.\n')
         

    if (config['hist_multi'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel histograms comparing various polarimetric vars.\n')
        outdir = outpath+'hist_multi/'
        os.makedirs(outdir,exist_ok=True)
 
        if config['wname'] in rdat.data.variables.keys():
            numr,ncol = 2,2
            figsize = (16,14)
            wspace = 0.25
        else:
            numr,ncol = 1,2
            figsize = (12,8)
            wspace = 0.5
        
        fig, ax = plt.subplots(numr,ncol,figsize=figsize,gridspec_kw={'wspace': wspace, 'hspace': 0.2, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.})
        axf = ax.flatten()

        zzdr_wrf,ed = rdat.hist2d(varx=rdat.dz_name,vary=rdat.zdr_name,binsx=config['dzbins'],binsy=config['drbins'])
        rdat.plot_2dhist(zzdr_wrf,ed,ax=axf[0])
        axf[0].set_xlabel(rdat.zdr_name+' '+rdat.units[rdat.zdr_name],fontsize=26)
        axf[0].set_ylabel(rdat.dz_name+' '+rdat.units[rdat.dz_name],fontsize=26,labelpad=0)
        #axf[0].set_title(title_string)
        axf[0].text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=12, color='k', zorder=10, weight='bold', transform=axf[0].transAxes) # (a) Top-left
        print('Panel 1: '+rdat.zdr_name+' vs. '+rdat.dz_name)

        zkdp_wrf,edk = rdat.hist2d(varx=rdat.dz_name,vary=rdat.kdp_name,binsx=config['dzbins'],binsy=config['kdbins'])
        rdat.plot_2dhist(zkdp_wrf,edk,ax=axf[1])
        #axf[1].set_title(title_string)
        axf[1].set_xlabel(rdat.kdp_name+' '+rdat.units[rdat.kdp_name],fontsize=26)
        axf[1].set_ylabel(rdat.dz_name+' '+rdat.units[rdat.dz_name],fontsize=26,labelpad=0)
        print('Panel 2: '+rdat.kdp_name+' vs. '+rdat.dz_name)

        if config['wname'] in rdat.data.variables.keys():
            zw_wrf,edw = rdat.hist2d(varx=rdat.dz_name,vary=rdat.w_name,binsx=config['dzbins'],binsy=config['wbins'])
            rdat.plot_2dhist(zw_wrf,edw,ax=axf[2])
            #axf[2].set_title(title_string)
            axf[2].set_xlabel(rdat.w_name+' '+rdat.units[rdat.w_name],fontsize=26)
            axf[2].set_ylabel(rdat.dz_name+' '+rdat.units[rdat.dz_name],fontsize=26,labelpad=0)
            print('Panel 3: '+rdat.w_name+' vs. '+rdat.dz_name)

            zr_wrf,edr = rdat.hist2d(varx=rdat.rr_name,vary=rdat.w_name,binsx=config['rrbins'],binsy=config['wbins'],xthr=0.00000)
            cb6 = rdat.plot_2dhist(zr_wrf,edr,ax=axf[3])
            #axf[3].set_title(title_string)
            axf[3].set_xlabel(rdat.w_name+' '+rdat.units[rdat.w_name],fontsize=26)
            axf[3].set_ylabel(rdat.rr_name+' '+rdat.units[rdat.rr_name],fontsize=26,labelpad=10)
            axf[3].set_ylim(0,50)
            print('Panel 4: '+rdat.w_name+' vs. '+rdat.rr_name)

        plt.savefig('{d}{p}_2dPDF_4panel.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype']),dpi=400,bbox_inches='tight')
        plt.clf()

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
 

    if (config['hid_prof'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating vertical profiles of water, graupel, hail and snow.')

        outdir = outpath+'vertical_profile/'
        os.makedirs(outdir,exist_ok=True)
 
        fig, ax = plt.subplots(1,1,figsize=(18,12))
        
        hts, mwrf_water_vert = rdat.hid_vertical_fraction(config['hidwater'],z_resolution =config['z_resolution'])
        hts, mwrf_graup_vert = rdat.hid_vertical_fraction(config['hidgraup'],z_resolution =config['z_resolution'])
        hts, mwrf_hail_vert = rdat.hid_vertical_fraction(config['hidhail'],z_resolution =config['z_resolution'])
        hts, mwrf_snow_vert = rdat.hid_vertical_fraction(config['hidsnow'],z_resolution =config['z_resolution'])
        lw = 4
        plt.plot(mwrf_water_vert,hts,color='b',label='water',lw=lw)
        plt.plot(mwrf_graup_vert,hts,color='g',label='graupel',lw=lw)
        plt.plot(mwrf_hail_vert,hts,color='r',label='hail',lw=lw)
        plt.plot(mwrf_snow_vert,hts,color = 'yellow',label='snow',lw=lw)
        ax.tick_params(axis='both',labelsize=22)
        plt.xlabel('Frequency (%)',fontsize=24)
        plt.ylabel('Height (km)',fontsize=24)
        #plt.title(title_string)
        plt.legend(loc='best',fontsize=22)

        ax.text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=24, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        plt.savefig('{d}{p}_HID_vertprof.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype']),dpi=400,bbox_inches='tight')
        plt.clf()

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
 

    if (config['up_width'] | config['all3']):
            
        if config['wname'] in rdat.data.variables.keys():
        
            print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating vertical profile of updraft width as a function of temperature.\n')
            
            tmp, m_warea_wrf = rdat.updraft_width_profile(thresh_dz=True)
            #print max(m_warea_wrf)

            fig, ax = plt.subplots(1,1,figsize=(18,12))

            plt.plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
            plt.ylim(20,-60)
            plt.xlabel('Updraft Width (km$^2$)',fontsize=24)
            plt.ylabel('Temperature (deg C)',fontsize=24)
            ax.tick_params(axis='both',labelsize=22)
            #plt.title(title_string)

            ax.text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=26, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left

            plt.savefig('{d}{p}_updraft_width_{y}_vertprof.{t}'.format(d=outdir,p=rdat.exper,t=config['ptype'],y=config['y']),dpi=400,bbox_inches='tight')
            plt.clf()

            print('\nDone! Saved to '+outdir)
            print('Moving on.\n')
     

    if (config['cappi_multi'] | config['all3']):
 
        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel CAPPIs for various polarimetric vars.')
        print('Plotting CAPPIs by time for variables '+str(eval(config['cappi_vars']))+'...')
        
        outdir = outpath+'cappi_multi/'
        os.makedirs(outdir,exist_ok=True)
 
        if not config['z'] == '': zspan = list(eval(str([config['z']])))
        else: zspan = rdat.data[rdat.z_name].values
        
        for z in zspan: 
            
            print('\nz = '+str(z))
            xlim = config['xlim']
            ylim = config['ylim']

            for ii in range(len(tms)):
                
                ts = tms[ii]
                print(ts)
                
                fig = rdat.cappi_multiplot(ts=ts,xlim=xlim,ylim=config['ylim'],z=z,res = config['cappi_vectres'],varlist = eval(config['cappi_vars']),latlon=config['latlon'],statpt=True,dattype=config['type']) #eval(config['cappi_contours']))
                #fig = rdat.cappi_multiplot(ts=ts,xlim=config['xlim'],ylim=config['ylim'],z=config['z'],res = config['cappi_vectres'],varlist = eval(config['cappi_vars']),vectors = eval(config['cvectors']),contours = None,statpt=True) #eval(config['cappi_contours']))

                nvars = len(eval(config['cappi_vars']))
                if nvars <=6:
                    yof = 0.01
                else:
                    yof = -0.02
                yof = -0.01
                xof = 0.01
                
                label_subplots(fig,yoff=yof,xoff=xof,size=16,nlabels=nvars,horizontalalignment='left',verticalalignment='top',color='k',bbox=dict(facecolor='w', edgecolor='w', pad=2.0),weight='bold')
 
                if not config['ptype'].startswith('mp4'):
                    plt.savefig('{i}{e}_multi_cappi_{t:%Y%m%d_%H%M%S}_{h}.{p}'.format(p=config['ptype'],i=outdir,e=rdat.exper,h=z,t=ts),dpi=400,bbox_inches='tight')
                else: 
                    if len(rdat.date) < 6:
                        plt.savefig('{i}{e}_multi_cappi_{t:%Y%m%d_%H%M%S}_{h}.png'.format(i=outdir,e=rdat.exper,h=z,t=ts),dpi=400,bbox_inches='tight')
                    else:
                        plt.savefig(outdir+'/fig'+str(ii).zfill(3)+'.png',dpi=400,bbox_inches='tight')
                
                plt.close()

            if config['ptype'].startswith('mp4') and len(rdat.date) >= 6:

                st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
                en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')

                os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdat.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_multi_cappi_{t1}-{t2}_{h}.mp4'.format(p=config['ptype'],e=rdat.exper,i=outdir,t1=st,t2=en,h=z))
                
                plt.close()
   
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
 
    
    if (config['cappi_individ'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating individual CAPPIs for various polarimetric vars.')
        print('Plotting CAPPIs by time...\n')
 
        if not config['z'] == '': zspan = list(eval(str([config['z']])))
        else: zspan = rdat.data[rdat.z_name].values

        for i,v in enumerate(eval(config['cappi_vars'])):
           
            print(v)

            if v is None:
                continue
            else:
        
                for z in zspan: 
                    
                    print('\nz = '+str(z))
            
                    for ii in range(len(tms)):
                        ts = tms[ii]
                        print(ts)

                        outdir = outpath+'cappi_individ/'+rdat.names_uc[v]+'/'
                        os.makedirs(outdir,exist_ok=True)
                        
                        if v.startswith('HID'): cbar = 2
                        else: cbar = 1
                        
                        #fig, ax = rdat.cappi(str(v),ts=ts,xlim=config['xlim'],ylim=config['ylim'],cbar=cbar,z=config['z'],res =config['cappi_vectres'],vectors = eval(config['cvectors'])[i],statpt=True)
                        fig, ax = rdat.cappi(v,ts=ts,xlim=config['xlim'],ylim=config['ylim'],cbar=cbar,z=z,latlon=config['latlon'],statpt=True,dattype=config['type'])
                        
                        ax.text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                        ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=ts), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                        ax.text(0.99, 0.99, 'z = {a} km'.format(a=config['z']), horizontalalignment='right',verticalalignment='top', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='none', pad=0.0))
                        
                        if not config['ptype'].startswith('mp4'):
                            plt.savefig('{i}{e}_{v}_individ_cappi_{t:%Y%m%d_%H%M%S}_{h}.{p}'.format(p=config['ptype'],i=outdir,e=rdat.exper,h=config['z'],t=ts,v=rdat.names_uc[v]),dpi=400,bbox_inches='tight')
                        else: 
                            if len(rdat.date) < 6:
                                plt.savefig('{i}{e}_{v}_individ_cappi_{t:%Y%m%d_%H%M%S}_{h}.png'.format(p=config['ptype'],i=outdir,e=rdat.exper,h=config['z'],t=ts,v=rdat.names_uc[v]),dpi=400,bbox_inches='tight')
                            else:
                                plt.savefig(outdir+'/fig'+str(ii).zfill(3)+'.png',dpi=400,bbox_inches='tight')
                        
                        plt.close()

                    if config['ptype'].startswith('mp4') and len(rdat.date) >= 6:

                        st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
                        en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')

                        os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdat.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_{v}_individ_cappi_{t1}-{t2}_{h}.mp4'.format(p=config['ptype'],e=rdat.exper,i=outdir,v=rdat.names_uc[v],t1=st,t2=en,h=z))
                    
            print('')

        print('Done! Saved to '+outdir)
        print('Moving on.\n')


    if (config['rhi_multi'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating multi-panel RHIs for various polarimetric vars.')
        print('Plotting RHIs at y = '+str(config['y'])+'km north of the radar by time for variables '+str(eval(config['cappi_vars']))+'...')
 
        outdir = outpath+'rhi_multi/'
        os.makedirs(outdir,exist_ok=True)

        if not config['y'] == '': yspan = list([config['y']])
        else: yspan = rdat.data[rdat.y_name].values[25::50]

        for y in yspan:

            print('\ny = '+str(y))

            for ii in range(len(tms)):
                ts = tms[ii]
                print(ts)

                if rdat.w_name is not None:
                    fig = rdat.xsec_multiplot(ts=ts,y=y,vectors=eval(config['rvectors']),res = config['rhi_vectres'],xlim=config['xlim'],zmax=config['zmax'],varlist=eval(config['rhi_vars']),latlon=config['latlon'])
                else:
                    fig = rdat.xsec_multiplot(ts=ts,y=y,xlim=config['xlim'],zmax=config['zmax'],varlist=eval(config['rhi_vars']),latlon=config['latlon'])
                
                nvars = len(eval(config['rhi_vars']))-1*(rdat.w_name is None)
                if nvars <= 6:
                    yof = 0.01
                else:
                    yof=-0.02
                yof = -0.01
                xof = 0.01
                
                label_subplots(fig,yoff=yof,xoff=xof,size=16,nlabels=nvars,horizontalalignment='left',verticalalignment='top',color='k',bbox=dict(facecolor='w', edgecolor='w', pad=2.0),weight='bold')
                
                if not config['ptype'].startswith('mp4'):
                    plt.savefig('{i}{e}_multi_rhi_{t:%Y%m%d_%H%M%S}_{h}.{p}'.format(p=config['ptype'],i=outdir,e=rdat.exper,h=y,t=ts),dpi=400,bbox_inches='tight')
                else: 
                    if len(rdat.date) < 6:
                        plt.savefig('{i}{e}_multi_rhi_{t:%Y%m%d_%H%M%S}_{h}.png'.format(i=outdir,e=rdat.exper,h=y,t=ts),dpi=400,bbox_inches='tight')
                    else:
                        plt.savefig(outdir+'/fig'+str(ii).zfill(3)+'.png',dpi=400,bbox_inches='tight')
                
                plt.close()

            if config['ptype'].startswith('mp4') and len(rdat.date) >= 6:

                st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
                en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')

                os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdat.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_multi_rhi_{t1}-{t2}_{h}.mp4'.format(p=config['ptype'],e=rdat.exper,i=outdir,t1=st,t2=en,h=y))
            
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
 

    if (config['rhi_individ'] | config['all3']):

        print('IN PLOT_DRIVER.MAKE_SINGLE_PLOTS... creating individual RHIs for various polarimetric vars.')
        print('Plotting RHIs at y = '+str(config['y'])+'km north of the radar by time.\n')

        if not config['y'] == '': yspan = list(eval(str([config['y']])))
        else: yspan = rdat.data[rdat.y_name].values[0::50]
     
        for i,v in enumerate(eval(config['rhi_vars'])):

            if v is None:
                continue
            else:
                
                print(v)
            
                for y in yspan:

                    print('\ny = '+str(y))

                    for ii in range(len(tms)):
                        ts = tms[ii]
                        print(ts)
           
                        outdir = outpath+'rhi_individ/'+rdat.names_uc[v]+'/'
                        os.makedirs(outdir,exist_ok=True)
                        
                        if v.startswith('HID'): cbar = 2
                        else: cbar = 1
                        
                        fig, ax = rdat.xsec(v,ts=ts,y=y,res = config['rhi_vectres'],xlim=config['xlim'],cbar=cbar,zmax=config['zmax'],latlon=config['latlon'])
                        #fig, ax = rdat.xsec(v,ts=ts,y=config['y'],vectors=eval(config['rvectors'])[i],res = config['rhi_vectres'],xlim=config['xlim'],cbar=cbar,zmax=config['zmax'])
                        
                        ax.text(0, 1, '{e} {r}'.format(e=rdat.exper,r=rdat.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                        ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=ts), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                        ax.text(0.99, 0.99, 'y = {a} km'.format(a=y), horizontalalignment='right',verticalalignment='top', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='none', pad=0.0))
                        
                        if not config['ptype'].startswith('mp4'):
                            plt.savefig('{i}{e}_{v}_individ_rhi_{t:%Y%m%d_%H%M%S}_{h}.{p}'.format(p=config['ptype'],i=outdir,e=rdat.exper,h=y,t=ts,v=rdat.names_uc[v]),dpi=400,bbox_inches='tight')
                        else:
                            if len(rdat.date) < 6:
                                plt.savefig('{i}{e}_{v}_individ_rhi_{t:%Y%m%d_%H%M%S}_{h}.png'.format(i=outdir,e=rdat.exper,h=y,t=ts,v=rdat.names_uc[v]),dpi=400,bbox_inches='tight')
                            else:
                                plt.savefig(outdir+'/fig'+str(ii).zfill(3)+'.png',dpi=400,bbox_inches='tight')
                        
                        plt.close()

                    if config['ptype'].startswith('mp4') and len(rdat.date) >= 6:
                        st = rdat.date[0].strftime('%Y%m%d_%H%M%S')
                        en = rdat.date[-1].strftime('%Y%m%d_%H%M%S')

                        os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdat.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_{v}_individ_rhi_{t1}-{t2}_{h}.mp4'.format(p=config['ptype'],e=rdat.exper,i=outdir,v=rdat.names_uc[v],t1=st,t2=en,h=y))

            print('')

        print('Done! Saved to '+outdir)
        print('Moving on.\n')


    if config['qr_cappi']:

        for ts in tms:
        
            print ("qr_cappi")
            fig = rdat.cappi_multiplot(z=config['z'],ts=ts,xlim=config['xlim'],ylim=config['ylim'],varlist=eval(config['mix_vars']))
            plt.savefig('{d}{p}_qcappi_6panel_{s:%Y%m%d%H%M%S}_{r}_{z}km.{t}'.format(d=outdir,p=rdat.exper,s=ts,r=rdat.band+'-band',t=config['ptype'],z=config['z']),dpi=300)
            plt.clf()
    

    if config['qr_rhi']:
        
        for ts in tms:
        
            print ("qr_rhi")
            fig = rdat.xsec_multiplot(ts=ts,y=config['y'],xlim=config['xlim'],varlist=eval(config['mix_vars']))
            plt.savefig('{d}{p}_qrhi_6panel_{s:%Y%m%d%H%M%S}_{r}_{y}.{t}'.format(d=outdir,p=rdat.exper,s=ts,r=rdat.band+'-band',t=config['ptype'],y=config['y']),dpi=300)
            plt.clf()
    

def subset_convstrat(data,rdat,zlev=1):
    cssum =(rdat.data[rdat.cs_name].max(dim='z'))
    stratsub=data.sel(z=slice(zlev,zlev+1)).where(cssum==1)
    convsub=data.sel(z=slice(zlev,zlev+1)).where(cssum==2)
    allsub=data.sel(z=slice(zlev,zlev+1)).where(cssum>0)
    return stratsub,convsub,allsub
    
def plot_timeseries(data,tm,ax,ls = '-',cs=False,rdat=None,thresh=-50,typ='',zlev=1,make_zeros=False,areas=False,domain_rel=False):
    data.values[data.values<thresh] = np.nan
    if make_zeros==True:
        data = data.fillna(0.0)
    if cs == True:
        
        sdat,cdat,adat = subset_convstrat(data,rdat,zlev=zlev)
#       print ('plotting')
        if areas == True:
            if domain_rel == True:
                gridsize = np.float(data.coords['x'].size*data.coords['y'].size)
                print(gridsize,'gridsize')
            else:
                gridsize = adat.count(dim=['z','y','x'])
            ax.plot(np.array(tm),adat.count(dim=['z','y','x'])/gridsize,color='k',label='Total {e}'.format(e=typ),ls=ls)
            ax.plot(np.array(tm),cdat.count(dim=['z','y','x'])/gridsize,color='r',label='Conv {e}'.format(e=typ),ls=ls)
            ax.plot(np.array(tm),sdat.count(dim=['z','y','x'])/gridsize,color='b',label='strat {e}'.format(e=typ),ls=ls)

        else:
#             print('lottin data')
#             print(tm)
#             print(np.nanmax(adat))
            ax.plot(np.array(tm),adat.mean(dim=['z','y','x'],skipna=True),color='k',label='Total {e}'.format(e=typ),ls=ls)
            ax.plot(np.array(tm),cdat.mean(dim=['z','y','x'],skipna=True),color='r',label='Conv {e}'.format(e=typ),ls=ls)
            ax.plot(np.array(tm),sdat.mean(dim=['z','y','x'],skipna=True),color='b',label='strat {e}'.format(e=typ),ls=ls)
        
        ax.legend(loc='best',fontsize=14)

    else:
        ax.plot(np.array(tm),data.where(data>thresh).mean(dim=['z','y','x'],skipna=True),color='k',label='Total')

    #ax.xaxis.set_major_formatter(hourFormatter)
    #ax.xaxis.set_major_locator(HourLocator(interval=1))
    #d=plt.setp(plt.gca().get_xticklabels(), rotation=45, horizontalalignment='right')
    ax.tick_params(axis='both',labelsize=14)
    ax.set_xlabel('Time (UTC)',fontsize=16)

    return ax#,adat,cdat,sdat

def plot_quartiles(data,q1,q2,q3,z,ax,c1='goldenrod',c2='r',c3='k',split_updn=False,ls = '-',typ=''):
    if split_updn == True:
        
        pdat=data.load()
        pdat.values[pdat.values<-100] = np.nan
        
        wup = pdat.where(pdat.values>0)
        wdn = pdat.where(pdat.values<0)
    
        wup50 = wup.quantile(q2,dim=['x','y','d'])
        wup90 = wup.quantile(q1,dim=['x','y','d'])
        wup99 = wup.quantile(q3,dim=['x','y','d'])


        wdn50 = wdn.quantile(1-q2,dim=['x','y','d'])
        wdn90 = wdn.quantile(1-q1,dim=['x','y','d'])
        wdn99 = wdn.quantile(1-q3,dim=['x','y','d'])                                               
        if 'd' in z.dims:
            zdat = z.sel(d=0).values        
        else:
            zdat = z.values
        ax.plot(wup50,zdat,color=c2,label='50th {e}'.format(e=typ),ls=ls)
        ax.plot(wup90,zdat,color=c1,label='90th {e}'.format(e=typ),ls=ls)
        ax.plot(wup99,zdat,color=c3,label='99th {e}'.format(e=typ),ls=ls)
        ax.legend(loc='best',fontsize=14)
        ax.plot(wdn50,zdat,color=c2)
        ax.plot(wdn90,zdat,color=c1)
        ax.plot(wdn99,zdat,color=c3)

    else:
        pdat =data.load()
        pdat.values[pdat.values<-100] = np.nan
        
    
        wup50 = pdat.quantile(q2,dim=['x','y','d'])
        wup90 = pdat.quantile(q1,dim=['x','y','d'])
        wup99 = pdat.quantile(q3,dim=['x','y','d'])

        if 'd' in z.dims:
            zdat = z.sel(d=0).values        
        else:
            zdat = z.values
        ax.plot(wup50,zdat,color=c2,label='50th {e}'.format(e=typ),ls=ls)
        ax.plot(wup90,zdat,color=c1,label='90th {e}'.format(e=typ),ls=ls)
        ax.plot(wup99,zdat,color=c3,label='99th {e}'.format(e=typ),ls=ls)
        ax.legend(loc='best',fontsize=14)

    ax.tick_params(axis='both',labelsize=14)
    ax.set_ylabel('Height (km)',fontsize=16)

    return ax
    
    
def plot_verprof(data,z,ax,c='r',lab='',split_updn=False,ls = '-',typ='',thresh=-50):
    if split_updn == True:
        
        pdat=data.load().copy()
        pdat.values[pdat.values<-100] = np.nan
        
        wup = pdat.where(data>0)
        wdn = pdat.where(data<0)
    
        wup50 = wup.mean(dim=['x','y','d'])

        wdn50 = wdn.mean(dim=['x','y','d'])

        if 'd' in z.dims:
            zdat = z.sel(d=0).values        
        else:
            zdat = z.values
        ax.plot(wup50,zdat,color=c,label='{l} {e}'.format(l=lab,e=typ),ls=ls)
        ax.plot(wdn50,zdat,color=c,label='{l} {e}'.format(l=lab,e=typ),ls=ls)
        ax.legend(loc='best',fontsize=14)

    else:
        pdat =data.load().copy()
        pdat.values[pdat.values<thresh] = np.nan
        
    
        wup50 = pdat.mean(dim=['x','y','d'])
        if 'd' in z.dims:
            zdat = z.sel(d=0).values        
        else:
            zdat = z.values
        ax.plot(wup50,zdat,color=c,label='{l} {e}'.format(l=lab,e=typ),ls=ls)
        ax.legend(loc='best',fontsize=14)

    ax.tick_params(axis='both',labelsize=14)
    ax.set_ylabel('Height (km)',fontsize=16)
    return ax
    
def hid_cdf(data,rdat,hts,species,z_resolution=1.0, ret_z=0,pick=None,z_ind =0, mask = None):
    # vertical HID_cdf with bar plots I think
    delz = hts[1]-hts[0]
    if np.mod(z_resolution, delz) != 0:
            print ('Need even multiple of vertical resolution: {d.1f}'.format(d = delz))
            return
    hold = deepcopy(data)

    if mask is not None:
#        print 'maskind HID data'
        hold[mask] = -1

    multiple = np.int(z_resolution/delz)

    # loop thru the species and just call the vertical hid volume
    all_vols = []
    for sp in range(len(species)):
        #print sp
        htsn, tdat = GF.vertical_hid_volume(hold,hts,delz,[sp+1], z_resolution=z_resolution, pick=pick,z_ind=0) # need the +1
        all_vols.append(tdat)
        
    print('htsn',np.shape(htsn))
    all_vols = np.array(all_vols)
    all_cdf = np.zeros_like(all_vols)
#9        print np.shape(all_vols)
#3    print min(all_vols)
    # shape is 10,16, which is nspecies x nheights
    # need to do cdf on each level
    all_vols[all_vols == np.nan] = 0.0
#    print max(all_vols)
    for iz in range(all_vols.shape[1]):
        # loop thru the vertical
#        print all_vols[:,iz]
#        print iz
        level_cum_vol = np.cumsum((all_vols[:, iz]))
#        if level_cum_vol[-1] != 0:
        all_cdf[:, iz] = 100.0*level_cum_vol/level_cum_vol[-1]
#        else:
#            all_cdf[:, iz] = 100.0*level_cum_vol/1.

#    all_cdf[np.isnan(all_cdf)] = 0.0
#    print max(all_cdf)
    if ret_z == 1:
    
        return htsn,all_cdf,all_vols
    else:
        return htsn,all_cdf#, all_vols

def plot_hid_cdf(data, hts,rconf=None, ax=None, pick=None):
    # this will just plot it
    if rconf is None:
        print ("sorry, need rconf to run properly")
        return
    #print np.shape(data)
    if ax is None:
        fig, ax = plt.subplots(1,1)
    else:
        fig = ax.get_figure()   

    fig.subplots_adjust(left = 0.07, top = 0.93, right = 0.8, bottom = 0.1)
    print("starting loop",np.shape(data),np.shape(hts))

    for i, vl in enumerate(hts):
        #print vl,i
#            print self.data[self.z_name].data[vl]
        #print data[0,:]
#        print vl, rconf.hid_colors[1],data[0,i]
        ax.barh(vl, data[0, i], left = 0., edgecolor = 'none', color = rconf.hid_colors[1]) 
#        print vl

        for spec in range(1, len(rconf.species)): # now looping thru the species to make bar plot
#             print rconf.hid_colors[spec+1]
#            print data[spec-1,i]
#            print spec, data[spec,i], data[spec-1,i]
            if data[spec-1,i] == np.nan:
                print( 'shoot')
            ax.barh(vl, data[spec, i], left = data[spec-1, i], \
            color = rconf.hid_colors[spec+1], edgecolor = 'none')
    ax.set_xlim(0,100)
    ax.set_xlabel('Cumulative frequency (%)')
    ax.set_ylabel('Height (km MSL)')
    # now have to do a custom colorbar?
    GF.HID_barplot_colorbar(rconf,fig)  # call separate HID colorbar function for bar plots

    return fig, ax 

def plot_hid_comparison_cfad(rdat1,rdat2,z_res=1.0,config=None,n1=None,n2=None,n3=None,cscfad=None,savefig=True):
    hidtest1,hts1 = rdat1.hid_cdf(z_resolution=z_res,cscfad=cscfad)
    hidtest2,hts2 = rdat2.hid_cdf(z_resolution=z_res,cscfad=cscfad)
    
    fig, bigax = plt.subplots(1,2,figsize=(14,8),gridspec_kw={'hspace': 0.05, 'wspace': 0.1})
    axf = bigax.flatten()
  
    fig,ax = rdat1.plot_hid_cdf(ylab=1,cbar=2,ax=axf[0])
    fig,ax = rdat2.plot_hid_cdf(ylab=0,cbar=0,ax=axf[1])
    '''    
    if n1 is None:
        n1 = rdat1.exper
    if n2 is None:
        n2 = rdat2.exper
    
    axf[0].set_title(n1)
    axf[1].set_title(n2)
    
    if n3 is None and cscfad is not None:
        n3=cscfad
    else:
        n3='ALL'
    plt.tight_layout()
    if savefig == True:
        if cscfad is not None:
            plt.savefig('{d}CFAD_diff_{e1}_{e2}_{c}_HID.{p}'.format(p=config['ptype'],d=outdir,c=cscfad,e1=rdat1.exper,e2=rdat2.exper),dpi=400,bbox_inches='tight')
        else:
            plt.savefig('{d}CFAD_diff_{e1}_{e2}_HID.{p}'.format(p=config['ptype'],d=outdir,e1=rdat1.exper,e2=rdat2.exper),dpi=400,bbox_inches='tight')
        return fig, axf
    else:
        return fig,axf
    '''
    return fig,axf

def cfad(data,rdat,zvals, var='zhh01',nbins=30,value_bins=None, multiple=1,ret_z=0,z_resolution=1.0,cscfad = False):
# pick a variable and do a CFAD for the cell
    if value_bins is None: # set a default if nothing is there
        value_bins = np.linspace(rdat.lims[var][0], rdat.lims[var][1], 20)
        nbins = len(value_bins)
    else:
        value_bins = np.arange(0,nbins+1,1)
    
    if 'd' in zvals.dims:
        zvals=zvals.sel(d=0)
    
    sz=len(zvals.values)
    print (sz,multiple)
    print('Shape data',np.shape(data))
    looped = np.arange(0, sz, multiple)
    cfad_out = np.zeros((sz//multiple, nbins-1))

    for ivl, vl in enumerate(looped[:-1]):
        v = zvals.values[vl]
        v2 = zvals.values[vl+multiple]
        try:
            dum = np.squeeze(data.sel(z=slice(v,v2)).values)
        except:
            dum = np.squeeze(data.sel(z=slice(vl,vl+multiple)).values)
        #print('shape dum', np.shape(dum))
        dum2 = np.where(np.isfinite(dum))
        lev_hist, edges = np.histogram(np.ravel(dum[dum2]), bins=value_bins, density=True) 
        lev_hist = 100.0*lev_hist/np.sum(lev_hist)
        if max(lev_hist) > 0:
            cfad_out[ivl, :] = lev_hist

    if ret_z == 1:
    
        return cfad_out,zvals.values[0][looped]
    else:
        return cfad_out,value_bins

def plot_cfad(fig,cfad,hts,vbins, ax, maxval=10.0, above=2.0, below=15.0, bins=None, 
        log=False, pick=None, z_resolution=1.0,levels=None,tspan =None,cont = False, rconf = None,mask = None,**kwargs):


    if hts is None:
        print ('please provide nominal heights to cfad_plot')
        return
    try:
        if 'd' in hts.dims:
            hts=hts.sel(d=0)
    except AttributeError as e:
        hts=hts
    if log:
        norm = colors.LogNorm(vmin=1e-5, vmax=1e2)
    else:
        norm = None

    # plot the CFAD
    cfad_ma = np.ma.masked_where(cfad==0, cfad)
    #print('CFAD shape',np.shape(cfad_ma))

    levs = [0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.]
    cols = ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet']
    if cont is True:
        try:
            #print(np.shape(cfad_ma),np.shape(hts),'ln 1283')
            pc = ax.contourf(vbins[0:-1],hts,cfad_ma,levs,colors=cols,extend = 'both')
        except Exception as e:
            print( 'Can not plot with exception {e}'.format(e=e))
            return ax
    else:

        if levels is not None:
            cmap, norm = from_levels_and_colors(levs,cols) # mention levels and colors here
            #print cmap
            pc = ax.pcolormesh(vbins, hts, cfad_ma, norm=norm, cmap=cmap)
        else:
            
#            pc = ax.pcolormesh(vbins, hts, cfad_ma, vmin=0, vmax=maxval, norm=norm, **kwargs)
            pc = ax.pcolormesh(vbins, hts, cfad_ma, vmin=0, vmax=maxval, norm=norm, **kwargs)

    #cb = plt.colorbar(pc, ax=ax)
    #cb.set_label('Frequency (%)',fontsize=16)
    lur,bur,wur,hur = ax.get_position().bounds
    cbar_ax_dims = [lur+wur+0.02,bur-0.001,0.03,hur]
    cbar_ax = fig.add_axes(cbar_ax_dims)
    cbt = plt.colorbar(pc,cax=cbar_ax)
    cbt.ax.tick_params(labelsize=16)
    cbt.set_label('Frequency (%)', fontsize=16, rotation=270, labelpad=20)

    ax.tick_params(axis='both',labelsize=14)
    ax.set_ylabel('Height (km MSL)',fontsize=16)
#        try:
    if rconf is not None:
        if var == 'DRC' or var == 'DRS':
            varn = rconf.zdr_name
        elif var == 'DZC' or var == 'DZS':
            varn = rconf.dz_name
        elif var == 'KDC' or var == 'KDS':
            varn = rconf.kdp_name
        elif var == 'WSvar' or var == 'WCvar':
            varn = rconf.w_name
        else:
            varn = var
#        print 'ln192',varn
        if varn in rconf.names.keys():
            
            ax.set_xlabel('{n} {u}'.format(n=rconf.names[varn], u=rconf.units[varn]))
            #print rconf.print_title(tm=tspan)
#            ax.set_title("{d}".format(d=rconf.print_title(tm=tspan)))
    #            ax.set_title('%s %s %s CFAD' % (self.print_date(), self.band+'-band', self.longnames[var]))
        else:
            ax.set_xlabel('{n}'.format(n=var))
            #print rconf.print_title(tm=tspan)
#            ax.set_title("{d}".format(d=rconf.print_title(tm=tspan)))
#        except:
#            pass

    return ax

######################################
##### plot_driver.PLOT_COMPOSITE #####
######################################

# Description: plot_composite overlays a Cartopy basemap with a colormesh plot of a given variable.

def plot_composite(rdat,var,time,resolution='10m',cs_over=False,statpt=False):
    
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    
    # (1) Create an array that is a copy of the variable of interest. Remove bad data.
    dat = deepcopy(rdat.data[var].sel(d=time))
    # dat = np.ma.masked_below(rdat.data[rdat.zdr_name].sel(d=time).values,-2)
    whbad = np.where(rdat.data['CSS'].sel(d=time).values<0)
    dat.values[whbad] = np.nan
    dat = np.squeeze(dat.values)
    dzcomp = np.nanmax(dat,axis=0)
    cs_arr = np.squeeze(np.nanmax(np.squeeze(rdat.data['CSS'].sel(d=time).values),axis=0))

    # (2) Find lat/lon array for the basemap. If not found, calculate it using get_latlon_fromxy().
    if not rdat.lat_name in rdat.data.keys():
        print('No latitude. Calculating....')
        rdat.get_latlon_fromxy()
        lats = rdat.data['lat']
        lons = rdat.data['lon']
    else:
        if 'd' in rdat.data[rdat.lat_name].dims:
            lats = rdat.data[rdat.lat_name].sel(d=time).values
            lons= rdat.data[rdat.lon_name].sel(d=time).values
        else:
            lats = rdat.data[rdat.lat_name].values
            lons= rdat.data[rdat.lon_name].values
        
    # Specifies the detail level of the map.
    # Options are '110m' (default), '50m', and '10m'
    # print(min(lons),max(lons))
    
    # (3) Extract the range of values in the variable of interest and derive a colourmap.
    if var in rdat.lims.keys(): # If the variable exists in the dataset:
        #print( 'var:',var)
        range_lim = rdat.lims[var][1] - rdat.lims[var][0]
    #          print np.shape(data), np.shape(xdat),np.shape(ydat)
    #            print 'in var',var
        #print **kwargs
        vmin=rdat.lims[var][0]
        vmax=rdat.lims[var][1]
        cmap=rdat.cmaps[var]
    else: # If not
        print ('unrecognized var',var)
        dat = rdat.data[var].data
        dat[dat<-900.0]=np.nan
        range_lim = np.nanmax(dat) - np.nanmin(dat)
        vmin=np.nanmin(dat)
        vmax=np.nanmax(dat)
        cmap = plt.cm.gist_ncar

    # (4) Initiate a new figure with Mercator projection.
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    
    lnspc = 2
    ltspc = 1
    
    # (5) Overlay variable data as a colormesh plot (with colorbar) using PlateCarree projection.
    cb = ax.pcolormesh(lons,lats,dzcomp,vmin=vmin,vmax=vmax,cmap=cmap,transform=ccrs.PlateCarree())
    if cs_over == True: ax.contour(lons,lats,cs_arr,levels=[0,1,2,3],linewidths=3,colors=['black','black'],transform=ccrs.PlateCarree())
    if statpt: ax.plot(rdat.lon_0,rdat.lat_0,markersize=12,marker='^',color='k',transform=ccrs.PlateCarree())
   
    lur,bur,wur,hur = ax.get_position().bounds
    cbar_ax_dims = [lur+wur+0.02,bur-0.001,0.03,hur]
    cbar_ax = fig.add_axes(cbar_ax_dims)
    cbt = plt.colorbar(cb,cax=cbar_ax)
    cbt.ax.tick_params(labelsize=16)
    if var.startswith('REF'): labtxt = 'Composite Reflectivity (dBZ)'
    else: labtxt = var
    cbt.set_label(labtxt, fontsize=16, rotation=270, labelpad=20)

    # (6) Make the figure look pretty! 
    ax.coastlines(resolution=resolution)

    minlon = np.floor(lons[(0,0)])
    maxlon = np.ceil(lons[(0,-1)])
    minlat = np.floor(lats[(0,0)])
    maxlat = np.ceil(lats[(-1,0)])
    lonticks = np.linspace(minlon,maxlon,int((maxlon-minlon)/lnspc+1))
    latticks = np.linspace(minlat,maxlat,int((maxlat-minlat)/ltspc+1))
    gl.xlocator = ticker.FixedLocator(lonticks)
    gl.ylocator = ticker.FixedLocator(latticks)
    ax.set_extent([minlon-0.001,maxlon+0.001,minlat-0.001,maxlat+0.001])
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1.5, alpha=0.75, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.x_inline = False
    gl.y_inline = False
    gl.xlabel_style = {'size': 16, 'color': 'black'}#,'rotation':-15}
    gl.ylabel_style = {'size': 16, 'color': 'black'}#,'rotation':-15}
    
    newax = fig.add_axes(ax.get_position(), frameon=False)
    newax.tick_params(axis='x', labelsize=0, length=0, pad=15)
    newax.tick_params(axis='y', labelsize=0, length=0, pad=45)
    newax.set_xlabel('Longitude',fontsize=16)
    newax.set_ylabel('Latitude',fontsize=16)

    return fig, ax
