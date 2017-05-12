import glob
import os
import sys
from netCDF4 import Dataset
import pandas as pd
import xarray as xr
import numpy as np
#import RadarData
import datetime
import matplotlib.pyplot as plt

from polarris_config import get_data
import RadarData


def make_single_pplots(rdat,tspan,flags,dir,exp='TWPICE',ty='png',extra='',z=2.0,y=-12.5):
    ts = tspan[0]
    te = tspan[1]

    if exp == 'TWPICE':
        xlim =[129.5,132.5]
        ylim=[-13.5,-10.5]
        if y != -12.5:
            y=y
        if z != 2.0:
            z=z
    if exp == 'MC3E':
        xlim = [-99,-95.5]
        ylim = [35.0,37.5]

        if y != -12.5:
            y=y
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
        plt.savefig('{d}{p}_CFAD_4panel_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        
        
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
        plt.savefig('{d}{p}_CFAD_W_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        rdat.cfad_plot(rdat.dz_name,ax = ax,bins=dzbins,z_resolution=1.0,levels='levs',tspan= tspan)
        plt.tight_layout()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        plt.savefig('{d}{p}_CFAD_dBZ_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        rdat.cfad_plot(rdat.zdr_name,ax= ax,bins=drbins,z_resolution=1.0,levels='levs',tspan= tspan)
        plt.tight_layout()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        plt.savefig('{d}{p}_CFAD_Zdr_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        rdat.cfad_plot(rdat.kdp_name,ax = ax,bins=drbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.tight_layout()

        fig, ax = plt.subplots(1,1,figsize=(18,12))
        plt.savefig('{d}{p}_CFAD_Kdp_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        rdat.cfad_plot(rdat.rho_name,ax = ax,bins=rhbins,z_resolution=1.0,levels='levs',tspan = tspan)
        plt.tight_layout()
        plt.savefig('{d}{p}_CFAD_RHO_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        
    if flags['hid_cfad_flag'] == True:
        fig, ax = rdat.plot_hid_cdf()
        plt.savefig('{d}{p}_CFAD_HID_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
        
    if flags['joint_flag'] == True:
        rrbins = np.linspace(1,142,71)

        fig, ax = plt.subplots(2,2,figsize=(12,12))
        axf = ax.flatten()

        zzdr_wrf,ed = rdat.hist2d(varx=rdat.dz_name,vary=rdat.zdr_name,binsx=dzbins,binsy=drbins)
        rdat.plot_2dhist(zzdr_wrf,ed,ax=axf[0])
        axf[0].set_xlabel('Zdr')
        axf[0].set_ylabel('dBZ')
        zkdp_wrf,edk = rdat.hist2d(varx=rdat.dz_name,vary=rdat.kdp_name,binsx=dzbins,binsy=drbins)
        rdat.plot_2dhist(zkdp_wrf,edk,ax=axf[1])
        axf[1].set_xlabel('Kdp')
        axf[1].set_ylabel('dBZ')


        zw_wrf,edw = rdat.hist2d(varx=rdat.dz_name,vary=rdat.w_name,binsx=dzbins,binsy=cbins)
        rdat.plot_2dhist(zw_wrf,edw,ax=axf[2])
        axf[2].set_xlabel('W')
        axf[2].set_ylabel('dBZ')

        zr_wrf,edr = rdat.hist2d(varx='RRB',vary=rdat.w_name,binsx=rrbins,binsy=cbins,xthr=0.00000)
        cb6 = rdat.plot_2dhist(zr_wrf,edr,ax=axf[3],cbon=True)
        axf[3].set_xlabel('W')
        axf[3].set_ylabel('RRB')
        axf[3].set_ylim(0,50)
        
        plt.savefig('{d}{p}_2dPDF_4panel_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)

    if flags['hid_prof'] == True:
        hidwater = [1,2,10]
        hidgraup = [7,8]
        hidhail = [9]
        hidsnow =[3,4,5,6]

        hts, mwrf_water_vert = rdat.hid_vertical_fraction(hidwater,z_resolution =0.5)
        hts, mwrf_graup_vert = rdat.hid_vertical_fraction(hidgraup,z_resolution =0.5)
        hts, mwrf_hail_vert = rdat.hid_vertical_fraction(hidhail,z_resolution =0.5)
        hts, mwrf_snow_vert = rdat.hid_vertical_fraction(hidsnow,z_resolution =0.5)

        plt.plot(mwrf_water_vert,hts,color='b')
        plt.plot(mwrf_graup_vert,hts,color='g')
        plt.plot(mwrf_hail_vert,hts,color='r')
        plt.plot(mwrf_snow_vert,hts,color = 'yellow')
        plt.savefig('{d}{p}_HID_prof_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty),dpi=300)
    if flags['all_cappi']== True:
        for t in rdat.date:
            #z=2.0
            #print xlim
            rdat.cappi_multiplot(ts=t,xlim=xlim,ylim=ylim,z=2.0)
            plt.savefig('{d}{p}_polcappi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{z}km.{t}'.format(d=dir,p=rdat.expr,s=t,r=rdat.radar_name,x=extra,t=ty,z=z),dpi=300)
    if flags['all_xsec']== True:
        for t in rdat.date:
            #y=-12.5
            rdat.xsec_multiplot(ts=t,y=y,vectors=True,res = [15,2],xlim=xlim)    
            plt.savefig('{d}{p}_polrhi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.expr,s=t,r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)
    if flags['up_width'] == True:
        tmp, m_warea_wrf = rdat.updraft_width_profile(thresh_dz=True)
        #print m_warea_wrf
        plt.plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
        plt.ylim(20,-60)
        plt.xlabel('Updraft Width (km$^2$)')
        plt.ylabel('Temperature (deg C)')
        plt.savefig('{d}{p}_upwidth_{s:%Y%m%d%H%M}-{e:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.expr,s=ts,e=te,r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)
    if flags['qr_cappi'] == True:
        mx_vars = ['qc','qr','qg','qi','qh']
        for t in rdat.date:
            rdat.cappi_multiplot(z=z,ts=t,xlim=xlim,ylim=ylim,varlist=mx_vars)
            plt.savefig('{d}{p}_qcappi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{z}km.{t}'.format(d=dir,p=rdat.expr,s=t,r=rdat.radar_name,x=extra,t=ty,z=z),dpi=300)
    if flags['qr_rhi'] == True:
        mx_vars = ['qc','qr','qg','qi','qh']
        for t in rdat.date:
            rdat.xsec_multiplot(ts=t,y=y,xlim=xlim,varlist=mx_vars)
            plt.savefig('{d}{p}_qrhi_6panel_{s:%Y%m%d%H%M}_{r}_{x}_{y}.{t}'.format(d=dir,p=rdat.expr,s=t,r=rdat.radar_name,x=extra,t=ty,y=y),dpi=300)

def make_compare_pplots(rdat1,rdat2,flags,dir,exp1='TWPICE',exp2='MC3E',ty='png',extra=''):
#    ts = tspan[0]
#    te = tspan[1]
    if flags['w_area_comp'] == True:
        fig, ax = plt.subplots(1,2,figsize=(18,9))
        axf = ax.flatten()

        tmp, m_warea_wrf = rdat1.updraft_width_profile(thresh_dz=True)
        #print m_warea_wrf
        axf[0].plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
        axf[0].set_ylim(20,-60)
        axf[0].set_xlabel('Updraft Width (km$^2$)')
        axf[0].set_ylabel('Temperature (deg C)')
        axf[0].set_title('{e} updraft width'.format(e=exp1))

        tmp, m_warea_wrf = rdat2.updraft_width_profile(thresh_dz=True)
        #print m_warea_wrf
        axf[1].plot(m_warea_wrf,tmp,color='k',label='Obs',lw=5)
        axf[1].set_ylim(20,-60)
        axf[1].set_xlabel('Updraft Width (km$^2$)')
        axf[1].set_ylabel('Temperature (deg C)')
        axf[1].set_title('{e} updraft width'.format(e=exp2))
        plt.savefig('{d}{e1}_{e2}_w_area_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)

    if flags['hid_vert_comp'] == True:
        fig, ax = plt.subplots(1,3,figsize=(18,9))
        axf = ax.flatten()

        hidwater = [1,2,10]
        hidgraup = [7,8]
        hidhail = [9]
        hidsnow =[3,4,5,6]

        hts, tw_water_vert1 = rdat1.hid_vertical_fraction(hidwater,z_resolution =1.0)
        hts,tw_graup_vert1 = rdat1.hid_vertical_fraction(hidgraup,z_resolution =1.0)
        hts, tw_hail_vert1 = rdat1.hid_vertical_fraction(hidhail,z_resolution =1.0)
        hts, tw_snow_vert1 = rdat1.hid_vertical_fraction(hidsnow,z_resolution =1.0)

        hts, tw_water_vert2 = rdat2.hid_vertical_fraction(hidwater,z_resolution =1.0)
        hts, tw_graup_vert2 = rdat2.hid_vertical_fraction(hidgraup,z_resolution =1.0)
        hts, tw_hail_vert2 = rdat2.hid_vertical_fraction(hidhail,z_resolution =1.0)
        hts, tw_snow_vert2 = rdat2.hid_vertical_fraction(hidsnow,z_resolution =1.0)

        lw=3
        axf[0].plot(tw_water_vert1,hts,color='blue',label='Water',lw=lw)
        axf[0].plot(tw_graup_vert1,hts,color='green',label='Graupel',lw=lw)
        axf[0].plot(tw_hail_vert1,hts,color='red',label='Hail',lw=lw)
        axf[0].plot(tw_snow_vert1,hts,color='goldenrod',label='Snow',lw=lw)
        axf[0].set_title('{e1} Hydromeor Freq.'.format(e1=exp1),fontsize=20)
        axf[0].set_xlabel('Frequency',fontsize=18)
        axf[0].legend(loc='best',fontsize=18)
        axf[0].set_ylabel('Height (km)',fontsize=18)
        axf[0].set_ylim(0,20)

        axf[1].plot(tw_water_vert2,hts,color='blue',label='Water',lw=lw)
        axf[1].plot(tw_graup_vert2,hts,color='green',label='Graupel',lw=lw)
        axf[1].plot(tw_hail_vert2,hts,color='red',label='Hail',lw=lw)
        axf[1].plot(tw_snow_vert2,hts,color='goldenrod',label='Snow',lw=lw)
        axf[1].set_title('{e1} Hydromeor Freq.'.format(e1=exp2),fontsize=20)
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
        axf[2].set_title('{e1}-{e2} Hydromeor Freq.'.format(e1=exp1,e2=exp2),fontsize=20)
        axf[2].set_xlabel('Frequency',fontsize=18)
        axf[2].legend(loc='best',fontsize=18)
        axf[2].set_ylabel('Height (km)',fontsize=18)
        axf[2].set_ylim(0,20)

        plt.savefig('{d}{e1}_{e2}_hid_vert_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)

    if flags['cfad_comp'] == True:
        fig, ax = plt.subplots(1,3,figsize=(18,9))
        axf = ax.flatten()

        cbins = np.arange(-25,26,0.5)
        dzbins = np.arange(-10,70,1)
        drbins = np.arange(-2,6,0.1)
        kdbins = np.arange(-2,4,0.05)
        rrbins = np.logspace(0.01,100.01,30)

        #W

        rdat1.cfad_plot(rdat1.w_name,ax = axf[0],bins=cbins,z_resolution=1.0,levels='levs')
        rdat2.cfad_plot(rdat2.w_name,ax = axf[1],bins=cbins,z_resolution=1.0,levels='levs')

        rdat1_cfad,radz= rdat1.cfad(rdat1.w_name,value_bins=cbins,z_resolution=1.0,ret_z=1)
        rdat2_cfad= rdat2.cfad(rdat2.w_name,value_bins=cbins,z_resolution=1.0)
        
        diff_cfad = rdat1_cfad - rdat2_cfad
        cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
        #print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))
    
        #levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
        levels=np.arange(-8,8.1,0.1)
        cb=axf[2].contourf(cbins[0:-1],radz,cfad_ma,levels,cmap='bwr',extend='both')

    #    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('Height (km MSL)',fontsize=18)
        axf[2].set_xlabel(rdat1.w_name,fontsize = 18)

        axf[2].set_title('{e1}-{e2} {v}'.format(e1=exp1,e2=exp2,v=rdat1.w_name))

        plt.tight_layout()
        plt.savefig('{d}{e1}_{e2}_W_CFAD_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)

        #dBZ
        fig, ax = plt.subplots(1,3,figsize=(18,9))
        axf = ax.flatten()

        rdat1.cfad_plot(rdat1.dz_name,ax = axf[0],bins=dzbins,z_resolution=1.0,levels='levs')
        rdat2.cfad_plot(rdat2.dz_name,ax = axf[1],bins=dzbins,z_resolution=1.0,levels='levs')

        rdat1_cfad,radz= rdat1.cfad(rdat1.dz_name,value_bins=dzbins,z_resolution=1.0,ret_z=1)
        rdat2_cfad= rdat2.cfad(rdat2.dz_name,value_bins=dzbins,z_resolution=1.0)
        
        diff_cfad = rdat1_cfad - rdat2_cfad
        cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
        #print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))
    
        #levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
        levels=np.arange(-8,8.1,0.1)
        cb=axf[2].contourf(dzbins[0:-1],radz,cfad_ma,levels,cmap='bwr',extend='both')

    #    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('Height (km MSL)',fontsize=18)
        axf[2].set_xlabel('%s %s' %(rdat1.names[rdat1.dz_name], rdat1.units[rdat1.dz_name]))

        axf[2].set_title('{e1}-{e2} {v}'.format(e1=exp1,e2=exp2,v=rdat1.names[rdat1.dz_name]))

        plt.tight_layout()
        plt.savefig('{d}{e1}_{e2}_dz_CFAD_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)


        #KDP
        fig, ax = plt.subplots(1,3,figsize=(18,9))
        axf = ax.flatten()

        rdat1.cfad_plot(rdat1.kdp_name,ax = axf[0],bins=kdbins,z_resolution=1.0,levels='levs')
        rdat2.cfad_plot(rdat2.kdp_name,ax = axf[1],bins=kdbins,z_resolution=1.0,levels='levs')

        rdat1_cfad,radz= rdat1.cfad(rdat1.kdp_name,value_bins=kdbins,z_resolution=1.0,ret_z=1)
        rdat2_cfad= rdat2.cfad(rdat2.kdp_name,value_bins=kdbins,z_resolution=1.0)
        
        diff_cfad = rdat1_cfad - rdat2_cfad
        cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
        #print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))
    
        #levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
        levels=np.arange(-8,8.1,0.1)
        cb=axf[2].contourf(kdbins[0:-1],radz,cfad_ma,levels,cmap='bwr',extend='both')

    #    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('Height (km MSL)',fontsize=18)
        axf[2].set_xlabel('%s %s' %(rdat1.names[rdat1.kdp_name], rdat1.units[rdat1.kdp_name]))

        axf[2].set_title('{e1}-{e2} {v}'.format(e1=exp1,e2=exp2,v=rdat1.names[rdat1.kdp_name]))

        plt.tight_layout()
        plt.savefig('{d}{e1}_{e2}_KDP_CFAD_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)


        #Zdr
        fig, ax = plt.subplots(1,3,figsize=(18,9))
        axf = ax.flatten()

        rdat1.cfad_plot(rdat1.zdr_name,ax = axf[0],bins=drbins,z_resolution=1.0,levels='levs')
        rdat2.cfad_plot(rdat2.zdr_name,ax = axf[1],bins=drbins,z_resolution=1.0,levels='levs')

        rdat1_cfad,radz= rdat1.cfad(rdat1.zdr_name,value_bins=drbins,z_resolution=1.0,ret_z=1)
        rdat2_cfad= rdat2.cfad(rdat2.zdr_name,value_bins=drbins,z_resolution=1.0)
        
        diff_cfad = rdat1_cfad - rdat2_cfad
        cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
        #print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))
    
        #levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
        levels=np.arange(-8,8.1,0.1)
        cb=axf[2].contourf(drbins[0:-1],radz,cfad_ma,levels,cmap='bwr',extend='both')

    #    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
        plt.colorbar(cb,ax=axf[2])
        axf[2].set_ylabel('Height (km MSL)',fontsize=18)
        axf[2].set_xlabel('%s %s' %(rdat1.names[rdat1.zdr_name], rdat1.units[rdat1.zdr_name]))

        axf[2].set_title('{e1}-{e2} {v}'.format(e1=exp1,e2=exp2,v=rdat1.names[rdat1.zdr_name]))

        plt.tight_layout()
        plt.savefig('{d}{e1}_{e2}_zdr_CFAD_compare_{x}.{t}'.format(d=dir,e2=exp2,e1=exp1,x=extra,t=ty),dpi=300)

