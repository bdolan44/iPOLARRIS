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


def plot_cfad_int(dat1,typ='dz',image_dir='./',ptype = 'png',n1 = None,extra='wrf'):
    fig, ax = plt.subplots(1,1,figsize=(12,8))
#    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

    dat1cnt = np.shape(dat1['hts'])[0]

    cfad1_all = np.sum(dat1['{t}cfad'.format(t=typ)],axis=0)/dat1cnt

    if typ == 'w':
        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=ax,cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

    else:
        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=ax,cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)
    ax.set_title(n1)

    ax.set_ylim(0,18)

    plt.tight_layout()

    plt.savefig('{id}CFAD_{tp}_{s}_int.{t}'.format(id=image_dir,s=n1,t=ptype,tp=typ.upper()),dpi=200)
    plt.clf()

def plot_hid_int(dat1,typ='hid',image_dir = './',ptype = 'png',extra= 'wrf',n1 = None):
    dat1cnt = np.shape(dat1['hts'])[0]
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)

    fig, ax = plt.subplots(1,1,figsize=(12,8))
    fig, ax = GF.plot_hid_cdf(np.nansum(dat1['hidcfad'],axis=0)/dat1cnt,dat1['hidhts'][0],ax=ax,rconf=dat1['rconf'])
    ax.set_title(n1)

    ax.set_ylim(0,18)

    plt.tight_layout()
    plt.savefig('{id}CFAD_{h}_{s}_int_{x}.{t}'.format(id=image_dir,h=typ.upper(),s=n1,x=extra,t=ptype),dpi=200)
    plt.clf()
def plot_hid_prof_int(dat1,typ='hid',image_dir = './',ptype = 'png',extra='wrf',n1 = None,n2 = None):
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
    plt.savefig('{d}{e1}_hid_vert_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
    plt.clf()

def plot_joint_int(dat1,typ='zzdr',image_dir ='./',ptype='png',extra='wrf',n1= None,n2=None):
    fig, ax = plt.subplots(1,1,figsize=(12,9))
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if typ == 'zzdr':              
#        print np.shape(dat1['edgzzdr'][0][0]),np.shape(np.nansum(dat1['histzzdr'],axis=0))
        cb6 = ax.contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],np.nansum(dat1['histzzdr'],axis=0))
        ax.set_xlabel('Zdr')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zzdr_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = ax.contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        ax.set_xlabel('Kdp')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zkdp_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = ax.contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        ax.set_xlabel('W (m/s)')
        ax.set_ylabel('dBZ')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_zw_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

    if typ == 'wr':
        cb6 = ax.contourf(dat1['edgwr'][0][0][:-1],dat1['edgwr'][0][1][:-1],np.nansum(dat1['histwr'],axis=0).T)
        ax.set_ylabel('RR (mm/hr)')
        ax.set_xlabel('W (M/s)')
        ax.set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax)
        plt.savefig('{d}{e1}_wr_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

def plot_upwidth_int(dat1,image_dir ='./',ptype='png',extra='wrf',n1= None):
        #print np.max(m_warea_wrf)
        plt.plot(np.nanmean(dat1['warea'],axis=0),dat1['wareat'][0],color='k',lw=5)
        plt.ylim(20,-60)
        plt.xlabel('Updraft Width (km$^2$)')
        plt.ylabel('Temperature (deg C)')
        plt.savefig('{d}{e1}_upwidth_int_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

def plot_upwidth(dat1,dat2,image_dir ='./',ptype='png',extra='wrf',n1= None,n2=None):
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

    axf[1].plot(np.nanmean(dat2['warea'],axis=0),dat2['wareat'][0],color='k',lw=5)
    axf[1].set_ylim(20,-60)
    axf[1].set_xlabel('Updraft Width (km$^2$)')
    axf[1].set_ylabel('Temperature (deg C)')

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    plt.savefig('{id}{e1}_{e2}_wpwidth_compare_{s}_{x}.{t}'.format(id=image_dir,e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,s=st_diff,x=extra,t=ptype),dpi=200)
    plt.clf()

def plot_joint_comp(dat1,dat2,typ='zzdr',image_dir ='./',ptype='png',extra='wrf',n1= None,n2=None):
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
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])

        cb6 = axf[1].contourf(dat2['edgzzdr'][0][1][:-1],dat2['edgzzdr'][0][0][:-1],np.nansum(dat2['histzzdr'],axis=0))
        axf[1].set_xlabel('Zdr')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=extra,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histzzdr'],axis=0)-np.nansum(dat2['histzzdr'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzzdr'][0][1][:-1],dat1['edgzzdr'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')

        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('dBZ',fontsize=18)
        axf[2].set_xlabel('Zdr',fontsize = 18)
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.savefig('{d}{e1}_{e2}_zzdr_comp_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

    if typ == 'zkdp':              
        cb6 = axf[0].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],np.nansum(dat1['histkdz'],axis=0))
        axf[0].set_xlabel('Kdp')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=ax[0])

        cb6 = axf[1].contourf(dat2['edgkdz'][0][1][:-1],dat2['edgkdz'][0][0][:-1],np.nansum(dat2['histkdz'],axis=0))
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=extra,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])

        diffdat = np.nansum(dat1['histkdz'],axis=0)-np.nansum(dat2['histkdz'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgkdz'][0][1][:-1],dat1['edgkdz'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('Kdp')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zkdp_comp_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

    if typ == 'zw':
        cb6 = axf[0].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],np.nansum(dat1['histzw'],axis=0))
        axf[0].set_xlabel('W (m/s)')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title('{e} {m} {x}'.format(e=dat1['rconf'].exper,x=extra,m=dat1['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[0])
        
        cb6 = axf[1].contourf(dat2['edgzw'][0][1][:-1],dat2['edgzw'][0][0][:-1],np.nansum(dat2['histzw'],axis=0))
        axf[1].set_xlabel('W (m/s)')
        axf[1].set_ylabel('dBZ')
        axf[1].set_title('{e} {m} {x}'.format(e=dat2['rconf'].exper,x=extra,m=dat2['rconf'].mphys))
        plt.colorbar(cb6,ax=axf[1])
        
        diffdat = np.nansum(dat1['histzw'],axis=0)-np.nansum(dat2['histzw'],axis=0)
        levels=np.arange(-2,2.1,0.1)
        cb=axf[2].contourf(dat1['edgzw'][0][1][:-1],dat1['edgzw'][0][0][:-1],diffdat,levels,cmap='bwr',extend='both')
        axf[2].set_xlabel('W (m/s)')
        axf[2].set_ylabel('dBZ')
        axf[2].set_title('{d}-{v}'.format(d=n1,v=n2))
        plt.colorbar(cb,ax=axf[2])

        plt.savefig('{d}{e1}_{e2}_zw_comp_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=extra,t=ptype),dpi=300)
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
        
        plt.savefig('{d}{e1}_{e2}_wr_comp_{x}.{t}'.format(d=image_dir,e1=dat1['rconf'].exper,e2=dat2['rconf'].exper,x=extra,t=ptype),dpi=300)
        plt.clf()

def plot_cfad_compare(dat1,dat2,typ='dz',image_dir='./',ptype = 'png',n1 = None,n2 = None,extra='wrf'):
    fig, ax = plt.subplots(1,3,figsize=(18,8))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    dat1cnt = np.shape(dat1['hts'])[0]
    dat2cnt = np.shape(dat2['hts'])[0]

    cfad1_all = np.sum(dat1['{t}cfad'.format(t=typ)],axis=0)/dat1cnt
    cfad2_all = np.sum(dat2['{t}cfad'.format(t=typ)],axis=0)/dat2cnt

    if typ == 'w':
        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=axf[0],cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

        fig, ax = GF.cfad_plot('{t}var'.format(t=typ.upper()),cfad = cfad2_all, hts = dat2['hts'][0],  bins = dat2['{t}bins'.format(t=typ)],ax=axf[1],cfad_on = 0,rconf = dat2['rconf'],tspan = dat2['time'],maxval=20,cont=True,levels = True)

    else:
        fig, ax = GF.cfad_plot(typ.upper(),cfad = cfad1_all, hts = dat1['hts'][0],  bins = dat1['{t}bins'.format(t=typ)],ax=axf[0],cfad_on = 0,rconf = dat1['rconf'],tspan = dat1['time'],maxval=20,cont=True,levels = True)

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
    axf[2].set_xlabel(dat1['rconf'].names['{tp}'.format(tp=typ.upper())],fontsize = 18)

    axf[2].set_title('{d-{v}'.format(d=n1,v=n2))

    plt.tight_layout()

    plt.savefig('{id}CFAD_{tp}_{s}.{t}'.format(id=image_dir,s=st_diff,t=ptype,tp=typ.upper()),dpi=200)
    plt.clf()

def plot_hid_2panel(dat1,dat2,typ='hid',image_dir = './',ptype = 'png',extra= 'wrf',n1 = None,n2 = None,):
    dat1cnt = np.shape(dat1['hts'])[0]
    dat2cnt = np.shape(dat2['hts'])[0]
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    fig, ax = plt.subplots(1,2,figsize=(18,8))
    axf = ax.flatten()
    fig, ax = GF.plot_hid_cdf(np.nansum(dat1['hidcfad'],axis=0)/dat1cnt,dat1['hidhts'][0],ax=axf[0],rconf=dat1['rconf'])
    axf[0].set_title(n1)
    fig, ax = GF.plot_hid_cdf(np.nansum(dat2['hidcfad'],axis=0)/dat2cnt,dat2['hidhts'][0],ax=axf[1],rconf=dat2['rconf'])
    axf[1].set_title(n2)

    axf[0].set_ylim(0,18)
    axf[1].set_ylim(0,18)

    plt.tight_layout()
    st_diff = '{e1}-{e2}'.format(e1=dat1['rconf'].exper,e2=dat2['rconf'].exper)
    plt.savefig('{id}CFAD_{h}_{s}_{x}.{t}'.format(id=image_dir,h=typ.upper(),s=st_diff,x=extra,t=ptype),dpi=200)
    plt.clf()




def plot_hid_profile(dat1,dat2,typ='hid',image_dir = './',ptype = 'png',extra='wrf',n1 = None,n2 = None):
    fig, ax = plt.subplots(1,3,figsize=(18,9))
    axf = ax.flatten()
    if n1 is None:
        n1 = '{e}_{x}'.format(e=dat1['rconf'].exper,x=dat1['rconf'].mphys)
    if n2 is None:
        n2 = '{e}_{x}'.format(e=dat2['rconf'].exper,x=dat2['rconf'].mphys)

    tw_water_vert1 = dat1['water_vert']
    tw_graup_vert1 = dat1['graup_vert']
    tw_hail_vert1 = dat1['hail_vert']
    tw_snow_vert1 = dat1['snow_vert']

    tw_water_vert2 = dat2['water_vert']
    tw_graup_vert2 = dat2['graup_vert']
    tw_hail_vert2 = dat2['hail_vert']
    tw_snow_vert2 = dat2['snow_vert']
    hts = dat1['hidhts']

    lw=3
    axf[0].plot(tw_water_vert1,hts,color='blue',label='Water',lw=lw)
    axf[0].plot(tw_graup_vert1,hts,color='green',label='Graupel',lw=lw)
    axf[0].plot(tw_hail_vert1,hts,color='red',label='Hail',lw=lw)
    axf[0].plot(tw_snow_vert1,hts,color='goldenrod',label='Snow',lw=lw)
    axf[0].set_title('{e1} Hydromeor Freq.'.format(e1=n1),fontsize=20)
    axf[0].set_xlabel('Frequency',fontsize=18)
    axf[0].legend(loc='best',fontsize=18)
    axf[0].set_ylabel('Height (km)',fontsize=18)
    axf[0].set_ylim(0,20)

    axf[1].plot(tw_water_vert2,hts,color='blue',label='Water',lw=lw)
    axf[1].plot(tw_graup_vert2,hts,color='green',label='Graupel',lw=lw)
    axf[1].plot(tw_hail_vert2,hts,color='red',label='Hail',lw=lw)
    axf[1].plot(tw_snow_vert2,hts,color='goldenrod',label='Snow',lw=lw)
    axf[1].set_title('{e1} Hydromeor Freq.'.format(e1=n2),fontsize=20)
    axf[1].set_xlabel('Frequency',fontsize=18)
    axf[1].legend(loc='best',fontsize=18)
    axf[1].set_ylabel('Height (km)',fontsize=18)
    axf[1].set_ylim(0,20)

    diff_water = tw_water_vert1-tw_water_vert2
    diff_graup = tw_graup_vert1-tw_graup_vert2
    diff_hail = tw_hail_vert1-tw_hail_vert2
    diff_snow = tw_snow_vert1-tw_snow_vert2
    

    axf[2].plot(diff_water,hts,color='blue',label='Water',lw=lw)
    axf[2].plot(diff_graup,hts,color='green',label='Graupel',lw=lw)
    axf[2].plot(diff_hail,hts,color='red',label='Hail',lw=lw)
    axf[2].plot(diff_snow,hts,color='goldenrod',label='Snow',lw=lw)
    axf[2].set_title('{e1}-{e2} Hydromeor Freq.'.format(e1=n1,e2=n2),fontsize=20)
    axf[2].set_xlabel('Frequency',fontsize=18)
    axf[2].legend(loc='best',fontsize=18)
    axf[2].set_ylabel('Height (km)',fontsize=18)
    axf[2].set_ylim(0,20)

    plt.savefig('{d}{e1}_{e2}_hid_vert_compare_{x}.{t}'.format(d=dir,e2=dat1['rconf'].exper,e1=dat2['rconf'].exper,x=extra,t=ptype),dpi=300)
    plt.clf() 

def plot_upstat(dat1,dat2,typ='hid',image_dir = './',ptype = 'png',extra='wrf',n1 = None,n2 = None):
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

        
def make_single_pplots(rdat,flags,dir='./',exp='TWPICE',ty='png',extra='test',z=2.0,y=None):
    tspan= [rdat.date,rdat.date]
    ts = tspan[0]
#    print ts
#    print rdat.exper
    te = tspan[1]
    title_string = '{e} {t} {d1:%Y%m%d-%H%M} {x}'.format(e=rdat.exper,t=rdat.mphys,d1=ts[0],x=extra)

    if exp == 'TWPICE':
        xlim =[129.5,132.5]
        ylim=[-13.5,-10.5]
        if y is None:
            y=-12.5
        if z != 2.0:
            z=z
    if exp == 'MC3E':
        xlim = [-99,-95.5]
        ylim = [35.0,37.5]

        if y is None:
            y=35.7
        if z != 2.0:
            z=z


    if flags['cfad_4panel_flag'] == True:
        fig, ax = plt.subplots(2,2,figsize=(18,12))
        axf = ax.flatten()

        cbins = np.arange(-25,26,0.5)
        dzbins = np.arange(-10,70,1)
        drbins = np.arange(-2,6,0.1)
        kdbins = np.arange(-2,4,0.05)
        rrbins = np.logspace(0.01,100.01,30)

        rdat.cfad_plot(rdat.w_name,ax = axf[0],bins=cbins,z_resolution=1.0,levels='levs',tspan = tspan)
        rdat.cfad_plot(rdat.dz_name,ax = axf[1],bins=dzbins,z_resolution=1.0,levels='levs',tspan= tspan)
        rdat.cfad_plot(rdat.zdr_name,ax= axf[2],bins=drbins,z_resolution=1.0,levels='levs',tspan= tspan)
        rdat.cfad_plot(rdat.kdp_name,ax = axf[3],bins=drbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.tight_layout()
#        print "{s:%Y%m%d%H%M}".format(s=ts[0])
        plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()
        
    if flags['cfad_individ_flag'] == True:
        fig, ax = plt.subplots(1,1,figsize=(18,12))
#        axf = ax.flatten()

        cbins = np.arange(-25,26,0.5)
        dzbins = np.arange(-10,70,1)
        rhbins = np.arange(0.5,1.01,0.01)
        drbins = np.arange(-2,6,0.1)
        kdbins = np.arange(-2,4,0.05)
        rrbins = np.logspace(0.01,100.01,30)

        rdat.cfad_plot(rdat.w_name,ax = ax,bins=cbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_W_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()
        

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.dz_name,ax = ax,bins=dzbins,z_resolution=1.0,levels='levs',tspan= tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_dBZ_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.zdr_name,ax= ax,bins=drbins,z_resolution=1.0,levels='levs',tspan= tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_Zdr_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.kdp_name,ax = ax,bins=drbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.savefig('{d}{p}_CFAD_Kdp_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.tight_layout()
        plt.clf()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.rho_name,ax = ax,bins=rhbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_RHO_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()
        
    if flags['hid_cfad_flag'] == True:
        fig, ax = rdat.plot_hid_cdf()
        plt.savefig('{d}{p}_CFAD_HID_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        
        plt.clf()
        
    if flags['joint_flag'] == True:
        rrbins = np.linspace(1,142,71)

        fig, ax = plt.subplots(2,2,figsize=(12,12))
        axf = ax.flatten()

        zzdr_wrf,ed = rdat.hist2d(varx=rdat.dz_name,vary=rdat.zdr_name,binsx=dzbins,binsy=drbins)
        rdat.plot_2dhist(zzdr_wrf,ed,ax=axf[0])
        axf[0].set_xlabel('Zdr')
        axf[0].set_ylabel('dBZ')
        axf[0].set_title(title_string)
        zkdp_wrf,edk = rdat.hist2d(varx=rdat.dz_name,vary=rdat.kdp_name,binsx=dzbins,binsy=drbins)
        rdat.plot_2dhist(zkdp_wrf,edk,ax=axf[1])
        axf[1].set_title(title_string)
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')


        zw_wrf,edw = rdat.hist2d(varx=rdat.dz_name,vary=rdat.w_name,binsx=dzbins,binsy=cbins)
        rdat.plot_2dhist(zw_wrf,edw,ax=axf[2])
        axf[2].set_title(title_string)
        axf[2].set_xlabel('W')
        axf[2].set_ylabel('dBZ')

        zr_wrf,edr = rdat.hist2d(varx='RRB',vary=rdat.w_name,binsx=rrbins,binsy=cbins,xthr=0.00000)
        cb6 = rdat.plot_2dhist(zr_wrf,edr,ax=axf[3],cbon=True)
        axf[3].set_title(title_string)
        axf[3].set_xlabel('W')
        axf[3].set_ylabel('RRB')
        axf[3].set_ylim(0,50)
        
        plt.savefig('{d}{p}_2dPDF_4panel_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()
        

    if flags['hid_prof'] == True:
        hidwater = [1,2,10]
        hidgraup = [7,8]
        hidhail = [9]
        hidsnow =[3,4,5,6]

        hts, mwrf_water_vert = rdat.hid_vertical_fraction(hidwater,z_resolution =0.5)
        hts, mwrf_graup_vert = rdat.hid_vertical_fraction(hidgraup,z_resolution =0.5)
        hts, mwrf_hail_vert = rdat.hid_vertical_fraction(hidhail,z_resolution =0.5)
        hts, mwrf_snow_vert = rdat.hid_vertical_fraction(hidsnow,z_resolution =0.5)

        plt.plot(mwrf_water_vert,hts,color='b',label='water')
        plt.plot(mwrf_graup_vert,hts,color='g',label='graupel')
        plt.plot(mwrf_hail_vert,hts,color='r',label='hail')
        plt.plot(mwrf_snow_vert,hts,color = 'yellow',label='snow')
        plt.xlabel('Frequency (%)')
        plt.ylabel('Height (km)')
        plt.title(title_string)
        plt.legend(loc = 'best')
        plt.savefig('{d}{p}_HID_prof_{s:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty),dpi=300)
        plt.clf()
        
    if flags['all_cappi']== True:
        #z=2.0
        #print xlim
        rdat.cappi_multiplot(ts=ts[0],xlim=xlim,ylim=ylim,z=2.0)
        plt.savefig('{d}{p}_polcappi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{z}km.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty,z=z),dpi=300)
        plt.clf()
        
    if flags['all_xsec']== True:
        #y=-12.5
        rdat.xsec_multiplot(ts=ts[0],y=y,vectors=True,res = [15,2],xlim=xlim)    
        plt.savefig('{d}{p}_polrhi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)
        plt.clf()
        
    if flags['up_width'] == True:
        tmp, m_warea_wrf = rdat.updraft_width_profile(thresh_dz=True)
        #print np.max(m_warea_wrf)
        plt.plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
        plt.ylim(20,-60)
        plt.xlabel('Updraft Width (km$^2$)')
        plt.ylabel('Temperature (deg C)')
        plt.title(title_string)
        plt.savefig('{d}{p}_upwidth_{s:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)
        plt.clf()
        
    if flags['qr_cappi'] == True:
        mx_vars = ['qc','qr','qg','qi','qh']
        #print type(rdat)
        rdat.cappi_multiplot(z=2.0,ts=ts[0],xlim=xlim,ylim=ylim,varlist=mx_vars)
        plt.savefig('{d}{p}_qcappi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{z}km.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty,z=z),dpi=300)
        plt.clf()
        
    if flags['qr_rhi'] == True:
        mx_vars = ['qc','qr','qg','qi','qh']
        rdat.xsec_multiplot(ts=ts[0],y=y,xlim=xlim,varlist=mx_vars)
        plt.savefig('{d}{p}_qrhi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.exper,s=ts[0],r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)
        plt.clf()
        


