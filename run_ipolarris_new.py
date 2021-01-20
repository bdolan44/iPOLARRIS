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
#from polarris_config import run_exper
#from polarris_config import get_data
import warnings
warnings.filterwarnings('ignore')
import GeneralFunctions as GF
from skewPy import SkewT
from collections import OrderedDict
from polarris_driver_new import polarris_driver
import os
import sys


configfile = sys.argv[1:]
#print sys.argv[1:]

rdata, config = polarris_driver(configfile)
#config['image_dir'] ='./'
print(config['extrax'],'EXTRA 1 is')
#########################################

if sys.argv[2:]:
    configfile1 = sys.argv[2:]
    rdata2, config2 = polarris_driver(configfile1)

    print('calculating CFAD differences')

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=False,cscfad=False)
    ax[0].set_title(rdata.exper)
    ax[1].set_title(rdata2.exper)
    ax[2].set_title("{e} - {v}".format(e=rdata.exper,v=rdata2.exper))    
    plt.suptitle("Reflectivity")
    plt.savefig('{d}CFAD_diff_{e1}_{e2}_{c}{l}_{x}.{p}'.format(p=config['ptype'],d=config['image_dir'],c='ALL',x=config['extrax'],e1=rdata.exper,e2=rdata2.exper,l='reflectivity'),dpi=400,bbox_inches='tight')
    plt.close()


    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,'Z$_{dr}$',config,config2,bins=np.arange(-2,8,0.2),savefig=True,cscfad=False)
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,'K$_{dp}$',config,config2,bins=np.arange(-2,6,0.2),savefig=True,cscfad=False)
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad=False)


    fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad=None)
    ##Convective

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=True,cscfad='convective')
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,'Z$_{dr}$',config,config2,bins=np.arange(-2,8,0.2),savefig=True,cscfad='convective')
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,'K$_{dp}$',config,config2,bins=np.arange(-2,6,0.2),savefig=True,cscfad='convective')
    
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad='convective')


    fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad='convective',savefig=True)

################################################################################
##################Now you can just start plotting!##############################
################################################################################

### To see the variables that are available to plot, type:

#rdata.data.keys()

else:
    ################################################################################
    ##Plot a composite reflectivity at a given time.


    #tdate = datetime.datetime(2011,5,23,22,00)
    # tdate = datetime.datetime(2006,1,23,18,0,0)
    # whdate = np.where(np.abs(tdate-np.array(rdata.date)) == np.min(np.abs(tdate-np.array(rdata.date))))
    print('In run_ipolarris...running the COMPOSITE figs.')
    for i,d in enumerate(np.array(rdata.date)):
        print('plotting composites by time....')
        fig, ax = plot_driver.plot_composite(rdata,rdata.dz_name,i,cs_over=True)
        print('made composite')
        rtimematch = d
        ax.set_title('{e} {r} composite {d:%Y%m%d %H%M}'.format(d=rtimematch,e=rdata.exper,r=rdata.radar_name))
        minlat = config['ylim'][0]
        maxlat = config['ylim'][1]
        minlon = config['xlim'][0]
        maxlon = config['xlim'][1]
        ax.set_extent([minlon, maxlon, minlat,maxlat])

        plt.tight_layout()
        plt.savefig('{i}Composite_{v}_{t:%Y%m%d%H%M}_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],v=rdata.dz_name,t=rtimematch,e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
        plt.close()

        print('plotting cappis at 1 km by time...')
        fig, ax = plt.subplots(1,1,figsize=(8,8))
        if 'd' in rdata.data[rdata.z_name].dims:
            try:
                whz = np.where(rdata.data[rdata.z_name].sel(d=i).values==config['z'])[0][0]
            except IndexError as ie:
                #print('checking z...',rdata.data[rdata.z_name].sel(d=i).values)
                zdiffs = np.median(np.diff(rdata.data[rdata.z_name].values))
                whz = np.where(np.isclose(rdata.data[rdata.z_name].sel(d=i).values,config['z'],rtol=zdiffs))[0][0]
                
        else:
            whz = np.where(rdata.data[rdata.z_name].values==config['z'])[0][0]
        print('whz in run 122',whz)
        rdata.cappi(rdata.dz_name,z=whz,ts=d,contour='CS',ax=ax)
        ax.set_title('CAPPI DZ {t:%Y%m%d_%M%H%S} {h} km'.format(t=d,h=rdata.data['z'].values[whz]))
        ax.set_xlim(config['xlim'][0],config['xlim'][1])
        ax.set_ylim(config['ylim'][0],config['ylim'][1])
#        ax.set_extent([minlon, maxlon, minlat,maxlat])
        plt.savefig('{i}DZ_CAPPI_{h}_{v}_{t:%Y%m%d%H%M}_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],h=config['z'],v=rdata.dz_name,t=rtimematch,e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
        plt.close()

        fig, ax = plt.subplots(1,1,figsize=(8,8))
    #    whz = np.where(rdata.data[rdata.z_name].values==config['z'])[0][0]
        rdata.cappi(rdata.rr_name,z=whz,ts=d,contour='CS',ax=ax)
        ax.set_xlim(config['xlim'][0],config['xlim'][1])
        ax.set_ylim(config['ylim'][0],config['ylim'][1])
        ax.set_title('CAPPI RR {t:%Y%m%d_%M%D%S} {h} km'.format(t=d,h=rdata.data['z'].values[2]))
        plt.savefig('{i}RR_CAPPI_{h}_{v}_{t:%Y%m%d%H%M}_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],h=config['z'],v=rdata.dz_name,t=rtimematch,e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
        plt.close()
  

    # tdate = datetime.datetime(2006,1,23,18,00)
    # whdate = np.where(np.abs(tdate-np.array(rdata.date)) == np.min(np.abs(tdate-np.array(rdata.date))))
    # fig, ax = plot_driver.plot_composite(rdata,rdata.cs_name,whdate[0][0])
    # rtimematch = rdata.date[whdate[0][0]]
    # ax.set_title('C/S composite {d:%Y%m%d %H%M}'.format(d=rtimematch))
    # plt.tight_layout()
    # plt.savefig('{i}Composite_{v}_{t:%Y%m%d%H%M}_{e}_{m}_{x}.{p}'.format(i=config['image_dir'],v=rdata.cs_name,t=rtimematch,e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
    # plt.clf()

    ################################################################################
    ##Calculate a timeseries for writing out
    rrstratu,rrconvu,rrallu = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=-0.1)
    rrstrat,rrconv,rrall = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=0.0)

    import csv
    tformat = '%Y%m%d-%H%M%S'
    with open('{i}{e}_rr_uncondmean_stats.txt'.format(i=config['image_dir'],e=config['exper']), mode='w') as csv_file:
        v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
        v_writer.writerow(['Date', 'Unc_Conv_RR', 'Unc_Strat_RR', 'Unc_Tot_RR'])
        for i,v in enumerate(rdata.date):
            print( v)
            tim = v.strftime(tformat)
            dum =[tim,rrconvu[i].values,rrstratu[i].values,rrallu[i].values]
            v_writer.writerow(dum)

    tformat = '%Y%m%d-%H%M%S'
    with open('{i}{e}_rr_condmean_stats.txt'.format(i=config['image_dir'],e=config['exper']), mode='w') as csv_file:
        v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
        v_writer.writerow(['Date', 'Conv_RR', 'Strat_RR', 'Tot_RR'])
        for i,v in enumerate(rdata.date):
            print (v)
            tim = v.strftime(tformat)
            dum =[tim,rrconv[i].values,rrstrat[i].values,rrall[i].values]
            v_writer.writerow(dum)

    conv = np.where(rdata.data[rdata.cs_name].values == 2)
    strat = np.where(rdata.data[rdata.cs_name].values == 1)
    hist, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values)),bins=np.logspace(-1,2.4,40))
    histc, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values[conv])),bins=np.logspace(-1,2.4,40))
    hists, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values[strat])),bins=np.logspace(-1,2.4,40))


    tformat = '%Y%m%d-%H%M%S'
    with open('{i}{e}_rr_histgram_{m}.txt'.format(i=config['image_dir'],e=config['exper'],m=config['mphys']), mode='w') as csv_file:
        v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
        v_writer.writerow(['Date', 'Con', 'Strat', 'Tot'])
        for i,v in enumerate(eg[:-1]):
            dum =[v,histc[i],hists[i],hist[i]]
            v_writer.writerow(dum)


    ###Areas

    rrstratu_area,rrconvu_area,rrallu_area = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=-0.1,areas=True)
    rrstrat_area,rrconv_area,rrall_area = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=0.0,areas=True)

    #grid_area=rdata.radar_area()
    rain_area = rdata.radar_area
    import csv
    tformat = '%Y%m%d-%H%M%S'
    with open('{i}{e}_domain_area_stats.txt'.format(i=config['image_dir'],e=config['exper']), mode='w') as csv_file:
        v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
        v_writer.writerow(['Date', 'Unc_Con', 'Unc_Strat', 'Unc_Tot'])
        for i,v in enumerate(rdata.date):
            print (v)
            tim = v.strftime(tformat)
            dum =[tim,rrconvu_area[i].values.astype(float)*rdata.dx*rdata.dy,rrstratu_area[i].values.astype(float)*rdata.dx*rdata.dy,rrallu_area[i].values.astype(float)*rdata.dx*rdata.dy]
            v_writer.writerow(dum)

    tformat = '%Y%m%d-%H%M%S'
    with open('{i}{e}_rel_frequency_stats.txt'.format(i=config['image_dir'],e=config['exper']), mode='w') as csv_file:
        v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
        v_writer.writerow(['Date', 'Conv', 'Strat', 'Tot'])
        for i,v in enumerate(rdata.date):
            print( v)
            tim = v.strftime(tformat)
            dum =[tim,rrconv[i].values*rdata.dx*rdata.dy/rain_area*100.,rrstrat[i].values*rdata.dx*rdata.dy/rain_area*100.,rrall[i].values*rdata.dx*rdata.dy/rain_area*100.]
            v_writer.writerow(dum)




    ################################################################################
    ##First make a timeseries of rain rate, unconditional and conditional. This puts strat, conv, and total on the same plot but you can split the out by putting cs==False.
    ## The conditional rain rate is achieved by sending threshold = 0.
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,make_zeros=False)
    ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,ls='--',typ='uncond',make_zeros=True)#,zlev=0)

    ax.set_ylabel('Rain Rate (mm/hr)')
    ax.set_title('Precipitation Timeseries ')
    plt.tight_layout()
    plt.savefig('{i}Precip_timeseries_convstrat_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
    plt.close()


    ############################################################################

    ################################################################################
    ##Next let's make quantile (50,90,99) plots of the vertical velocity. This splits it by up and down, but you can turn split_updn == False
    if rdata.w_name is not None:
        fig,ax = plt.subplots(1,1,figsize=(10,10))
        ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=True)
        ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=False)
        ax.set_xlabel('Vertical velocity m/s')
        ax.set_title('Vertical velocity profiles')
        plt.tight_layout()
        plt.savefig('{i}Quantile_vvel_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
        plt.close()

        p99u,p90u,p50u,ht = rdata.percentile(wup=True)
        p99d,p90d,p50d,ht = rdata.percentile(wdown=True)
        p99a,p90a,p50a,ht = rdata.percentile(wdown=False)

        file = open('{i}{e}_{m}_updown_percentiles.txt'.format(i=config['image_dir'],e=rdata.exper,m=rdata.mphys),'w') 
 
        file.write("Updraft\n") 
        file.write("Height (km).    P99.    P90.     P50\n") 
        for i,h in enumerate(ht):
            file.write("{h}   {p1}   {p2}   {p3}\n".format(h=h,p1=p99u[i],p2=p90u[i],p3=p50u[i])) 

        file.write("Downdraft\n") 
        file.write("Height (km).    P99.    P90.     P50\n") 
        for i,h in enumerate(ht):
            file.write("{h}   {p1}   {p2}   {p3}\n".format(h=h,p1=p99d[i],p2=p90d[i],p3=p50d[i])) 

        file.write("ALL\n") 
        file.write("Height (km).    P99.    P90.     P50\n") 
        for i,h in enumerate(ht):
            file.write("{h}   {p1}   {p2}   {p3}\n".format(h=h,p1=p99a[i],p2=p90a[i],p3=p50a[i])) 


        file.close()
    else:
        print("No vertical velocity data.")
    ################################################################################

    ################################################################################
    ##Next let's make mean vertical profile of reflectivity
    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax = plot_driver.plot_verprof(rdata.data[rdata.dz_name],rdata.data[rdata.z_name],ax,split_updn=False,lab='dz',thresh=-50)
    ax.set_title('Vertical profile of reflectivity')
    ax.set_xlabel('Reflectivity')
    plt.tight_layout()
    plt.savefig('{i}MeanProfile_refl_{e}_{m}_{x}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)

    plt.close()
    ################################################################################
    ##Next let's make a reflectivity CFAD

#    cfaddat,vbins = plot_driver.cfad(rdata.data[rdata.dz_name],rdata,rdata.data[rdata.z_name],var=rdata.dz_name,nbins=40)
    cfaddat,vbins,r1ht = rdata.cfad(rdata.dz_name,ret_z=1,z_resolution=1.0,value_bins=np.arange(0,82,2),cscfad=False)

    fig,ax = plt.subplots(1,1,figsize=(10,10))
    ax = plot_driver.plot_cfad(cfaddat, hts = r1ht,  vbins = vbins,ax=ax,cfad_on = 0,tspan = config['date'],maxval=20,cont=True,levels = True)

    ax.set_xlabel('Reflectivity')
    ax.set_ylabel('Height (km)')
    ax.set_title('{c} CFAD'.format(c=rdata.exper))
    plt.tight_layout()
    plt.savefig('{i}CFAD_refl_{e}_{m}_{x}_new.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys,x=config['extrax']),dpi=400)
    plt.close()
    
    
    flags = {}
    for k in eval(config['ks']):
        flags[k]=config[k]

    if any(flags.values()) == True:
        plot_driver.make_single_pplots(rdata,flags,config)

    
