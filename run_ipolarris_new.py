#===================================================
#============== RUN_IPOLARRIS_NEW.PY ===============
#===================================================

# Import core Python packages
from collections import OrderedDict
import csv
import datetime
import glob
import matplotlib
matplotlib.use('Agg')
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import sys
import time
import warnings
warnings.filterwarnings('ignore')
import xarray as xr

# Import iPOLARRIS functions
import GeneralFunctions as GF
from polarris_driver_new import polarris_driver
import plot_driver
import RadarData
import RadarConfig
from skewPy import SkewT

#--------------- Main Program ----------------

print('\niPOLARRIS INITIATING...')
time.sleep(3)

print('\n#############################################')
print('####### Starting run_ipolarris_new.py #######')
print('#############################################')

configfile = sys.argv[1:] # Feed config file name as arg
#print sys.argv[1:]

print('\n#################################################')
print('############ Calling polarris_driver_new.py ####')
print('#################################################')
time.sleep(3)

rdata, config, config['uname'], config['vname'], config['wname'] = polarris_driver(configfile)

print('\n#################################################')
print('####### Returning to run_ipolarris_new.py #######')
print('#################################################')

config['image_dir'] = config['image_dir']+\
        config['exper']+'_'+config['sdatetime']+'_'+config['edatetime']+'/'+\
        config['type']+'/'
#print(,'EXTRA 1 is')

# If a second argument is passed for WRF config file, produce a bunch of comparison plots!
# More comments in this section TBD!
if sys.argv[2:]:
    configfile1 = sys.argv[2:]
    rdata2, config2 = polarris_driver(configfile1)

    print('calculating CFAD differences')

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=False,cscfad=False)
    ax[0].set_title(rdata.exper)
    ax[1].set_title(rdata2.exper)
    ax[2].set_title("{e} - {v}".format(e=rdata.exper,v=rdata2.exper))    
    plt.suptitle("Reflectivity")
    plt.savefig('{d}CFAD_diff_{e1}_{e2}_{c}{l}.{p}'.format(p=config['ptype'],d=config['image_dir'],c='ALL',e1=rdata.exper,e2=rdata2.exper,l='reflectivity'),dpi=400,bbox_inches='tight')
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
################## Now you can just start plotting! ############################
################################################################################

### To see the variables that are available to plot, type:

#rdata.data.keys()

else:
    ################################################################################
    ##Plot a composite reflectivity at a given time.

    # tdate = datetime.datetime(2011,5,23,22,00)
    # tdate = datetime.datetime(2006,1,23,18,0,0)
    # whdate = np.where(np.abs(tdate-np.array(rdata.date)) == np.min(np.abs(tdate-np.array(rdata.date))))

    if (config['compo_ref'] | config['all1']):
    
        print('IN RUN_IPOLARRIS_NEW... creating COMPOSITE figures.')
        print('\nPlotting composites by time for variable '+rdata.dz_name+'...')
            
        outdir = config['image_dir']+'composite_'+rdata.dz_name+'/'
        os.makedirs(outdir,exist_ok=True)
       
        for i,rtimematch in enumerate(np.array(rdata.date)):

            fig, ax = plot_driver.plot_composite(rdata,rdata.dz_name,i,cs_over=False,statpt=True)
            #ax.set_title('{e} {r} Composite {d:%Y-%m-%d %H%M} UTC'.format(d=rtimematch,e=rdata.exper,r=rdata.radar_name), fontsize=18)
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=rtimematch), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            
            #minlat = config['ylim'][0]
            #maxlat = config['ylim'][1]
            #minlon = config['xlim'][0]
            #maxlon = config['xlim'][1]
            #ax.set_extent([minlon, maxlon, minlat,maxlat])

            #plt.tight_layout()
            #plt.savefig('{i}composite_{v}_{d:%Y-%m-%d_%H%M%S}_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir']+'composite_'+rdata.dz_name+'/',v=rdata.dz_name,d=rtimematch,e=rdata.exper,m=rdata.mphys),dpi=400, bbox_inches='tight') 
            plt.savefig('{i}{e}_{v}_{d:%Y-%m-%d_%H%M%S}.{p}'.format(p=config['ptype'],e=rdata.exper,i=outdir,d=rtimematch,v=rdata.dz_name),dpi=400,bbox_inches='tight')
            plt.close()
            print(rtimematch)
        
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
        
    if (config['cappi_ref'] | config['all1']):
        
        print('\nIN RUN_IPOLARRIS_NEW... creating CAPPI figures.')
        print('\nPlotting CAPPIs at height z = '+str(config['z'])+'km by time for variable '+rdata.dz_name+'...')
 
        outdir = config['image_dir']+'cappi_'+rdata.dz_name+'/'
        os.makedirs(outdir,exist_ok=True)
       
        for i,rtimematch in enumerate(np.array(rdata.date)):

            fig, ax = plot_driver.plot_cappi(rdata,rdata.dz_name,rdata.z_name,config['z'],i,rtimematch,cs_over=False,statpt=True)
                    
            #ax.set_xlim(config['xlim'][0],config['xlim'][1])
            #ax.set_ylim(config['ylim'][0],config['ylim'][1])           
            #ax.set_extent([minlon, maxlon, minlat,maxlat])
            
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=rtimematch), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(0.99, 0.99, 'z = {a} km'.format(a=config['z']), horizontalalignment='right',verticalalignment='top', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes)

            #plt.savefig('{i}dz_cappi_{h}_{v}_{t:%Y-%m-%d_%H%M%S}_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir']+'cappi_'+rdata.dz_name+'/',h=config['z'],v=rdata.dz_name,t=rtimematch,e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight') 
            plt.savefig('{i}{e}_{v}_cappi_{h}_{t:%Y-%m-%d_%H%M%S}.{p}'.format(p=config['ptype'],i=outdir,e=rdata.exper,h=config['z'],v=rdata.dz_name,t=rtimematch),dpi=400,bbox_inches='tight')
            plt.close()

            print(rtimematch)

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
    
    if (config['cappi_rr'] | config['all1']):

        print('\nIN RUN_IPOLARRIS_NEW... creating CAPPI figures.')
        print('\nPlotting CAPPIs at height z = '+str(config['z'])+'km by time for variable '+rdata.dz_name+'...')
 
        print('Plotting CAPPIs at height z = '+str(config['z'])+'km by time for variable '+rdata.rr_name+'...')
 
        outdir = config['image_dir']+'cappi_'+rdata.rr_name+'/'
        os.makedirs(outdir,exist_ok=True)

        for i,rtimematch in enumerate(np.array(rdata.date)):

            fig, ax = plot_driver.plot_cappi(rdata,rdata.rr_name,rdata.z_name,config['z'],i,rtimematch,cs_over=False,statpt=True)
            #fig, ax = plt.subplots(1,1,figsize=(10,8))
        #    whz = np.where(rdata.data[rdata.z_name].values==config['z'])[0][0]
            #rdata.cappi(rdata.rr_name,z=whz,ts=rtimematch,contour='CS',ax=ax)
            #ax.set_xlim(config['xlim'][0],config['xlim'][1])
            #ax.set_ylim(config['ylim'][0],config['ylim'][1])
            #ax.set_title('CAPPI RR {t:%Y%m%d_%M%D%S} {h} km'.format(t=rtimematch,h=rdata.data['z'].values[2]))
            
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=rtimematch), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(0.99, 0.99, 'z = {a} km'.format(a=config['z']), horizontalalignment='right',verticalalignment='top', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes)

            #plt.savefig('{i}rr_cappi_{h}_{v}_{t:%Y-%m-%d_%H%M%S}_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir']+'cappi_'+rdata.rr_name+'/',h=config['z'],v=rdata.dz_name,t=rtimematch,e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight') 
            plt.savefig('{i}{e}_{v}_cappi_{h}_{t:%Y-%m-%d_%H%M%S}.{p}'.format(p=config['ptype'],i=outdir,e=rdata.exper,h=config['z'],v=rdata.rr_name,t=rtimematch),dpi=400,bbox_inches='tight')
            plt.close()

            print(rtimematch)

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')

    # tdate = datetime.datetime(2006,1,23,18,00)
    # whdate = np.where(np.abs(tdate-np.array(rdata.date)) == np.min(np.abs(tdate-np.array(rdata.date))))
    # fig, ax = plot_driver.plot_composite(rdata,rdata.cs_name,whdate[0][0])
    # rtimematch = rdata.date[whdate[0][0]]
    # ax.set_title('C/S composite {d:%Y%m%d %H%M}'.format(d=rtimematch))
    # plt.tight_layout()
    # plt.savefig('{i}Composite_{v}_{t:%Y%m%d%H%M}_{e}_{m}_{x}.{p}'.format(i=config['image_dir'],v=rdata.cs_name,t=rtimematch,e=rdata.exper,m=rdata.mphys,x=),dpi=400)
    # plt.clf()

    ################################################################################
    
    if (config['rrstats_txt'] | config['all2']):
        
        print('\nIN RUN_IPOLARRIS_NEW... creating text files.')
        print('Printing unconditional-mean statistics for variable '+rdata.rr_name+'...')
        
        ##Calculate a timeseries for writing out
        rrstratu,rrconvu,rrallu = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2.5,cs_flag=True,thresh=-0.1)
        rrstrat,rrconv,rrall = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2.5,cs_flag=True,thresh=0.0)

        tformat = '%Y%m%d-%H%M%S'
        outdir = config['image_dir']+'txtfiles/'
        os.makedirs(outdir,exist_ok=True)
        with open('{i}{e}_{v}_uncondmean_stats.txt'.format(i=outdir,v=rdata.rr_name,e=rdata.exper),mode='w') as csv_file:
            v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
            v_writer.writerow(['Date', 'Unc_Conv_RR', 'Unc_Strat_RR', 'Unc_Tot_RR'])
            for i,v in enumerate(rdata.date):
                #print( v)
                tim = v.strftime(tformat)
                dum =[tim,rrconvu[i].values,rrstratu[i].values,rrallu[i].values]
                v_writer.writerow(dum)

        print('\nDone! Saved to '+outdir)
        print('Printing conditional-mean statistics for variable '+rdata.rr_name+'...')

        with open('{i}{e}_{v}_condmean_stats.txt'.format(i=outdir,v=rdata.rr_name,e=rdata.exper),mode='w') as csv_file:
            v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
            v_writer.writerow(['Date', 'Conv_RR', 'Strat_RR', 'Tot_RR'])
            for i,v in enumerate(rdata.date):
                #print (v)
                tim = v.strftime(tformat)
                dum =[tim,rrconv[i].values,rrstrat[i].values,rrall[i].values]
                v_writer.writerow(dum)

        print('\nDone! Saved to '+outdir)
        print('Printing relative frequency statistics for variable '+rdata.rr_name+'...')

        rain_area = rdata.radar_area
        with open('{i}{e}_{v}_rel_frequency_stats.txt'.format(i=outdir,v=rdata.rr_name,e=rdata.exper), mode='w') as csv_file:
            v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
            v_writer.writerow(['Date', 'Conv', 'Strat', 'Tot'])
            for i,v in enumerate(rdata.date):
                #print( v)
                tim = v.strftime(tformat)
                dum =[tim,rrconv[i].values*rdata.dx*rdata.dy/rain_area*100.,rrstrat[i].values*rdata.dx*rdata.dy/rain_area*100.,rrall[i].values*rdata.dx*rdata.dy/rain_area*100.]
                v_writer.writerow(dum)

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')

    if (config['rrhist_txt'] | config['all2']):

        print('\nIN RUN_IPOLARRIS_NEW... creating text files.')
        print('Printing histogram data for variable '+rdata.rr_name+'...')
 
        conv = np.where(rdata.data[rdata.cs_name].values == 2)
        strat = np.where(rdata.data[rdata.cs_name].values == 1)
        hist, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values)),bins=np.logspace(-1,2.4,40))
        histc, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values[conv])),bins=np.logspace(-1,2.4,40))
        hists, eg = np.histogram(np.ravel((rdata.data[rdata.rr_name].values[strat])),bins=np.logspace(-1,2.4,40))

        #tformat = '%Y%m%d-%H%M%S'
        outdir = config['image_dir']+'txtfiles/'
        os.makedirs(outdir,exist_ok=True)
        with open('{i}{e}_{v}_rr_histgram_{m}.txt'.format(i=outdir,v=rdata.rr_name,e=rdata.exper,m=rdata.mphys), mode='w') as csv_file:
            v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
            v_writer.writerow(['Date', 'Con', 'Strat', 'Tot'])
            for i,v in enumerate(eg[:-1]):
                dum =[v,histc[i],hists[i],hist[i]]
                v_writer.writerow(dum)
 
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
   
    if (config['rrstats_areas_txt'] | config['all2']):
        ###Areas
        print('\nIN RUN_IPOLARRIS_NEW... creating text files.')
        print('Printing domain area statistics for '+rdata.rr_name+'...')
 
        rrstratu_area,rrconvu_area,rrallu_area = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=-0.1,areas=True)
        rrstrat_area,rrconv_area,rrall_area = rdata.calc_timeseries_stats(rdata.rr_name,ht_lev=2,cs_flag=True,thresh=0.0,areas=True)

        #grid_area=rdata.radar_area()
        rain_area = rdata.radar_area
        tformat = '%Y%m%d-%H%M%S'
        outdir = config['image_dir']+'txtfiles/' 
        os.makedirs(outdir,exist_ok=True)
        with open('{i}{e}_domain_area_stats.txt'.format(i=outdir,e=rdata.exper), mode='w') as csv_file:
            v_writer = csv.writer(csv_file, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_NONNUMERIC)
            v_writer.writerow(['Date', 'Unc_Con', 'Unc_Strat', 'Unc_Tot'])
            for i,v in enumerate(rdata.date):
                print (v)
                tim = v.strftime(tformat)
                dum =[tim,rrconvu_area[i].values.astype(float)*rdata.dx*rdata.dy,rrstratu_area[i].values.astype(float)*rdata.dx*rdata.dy,rrallu_area[i].values.astype(float)*rdata.dx*rdata.dy]
                v_writer.writerow(dum)

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')

    if (config['rr_timeseries'] | config['all1']):

        print('\nIN RUN_IPOLARRIS_NEW... creating timeseries.')
        print('Plotting timeseries for variable '+rdata.rr_name+'...')
 
        ################################################################################
        ##First make a timeseries of rain rate, unconditional and conditional. This puts strat, conv, and total on the same plot but you can split the out by putting cs==False.
        ## The conditional rain rate is achieved by sending threshold = 0.
        fig,ax = plt.subplots(1,1,figsize=(12,8))
        ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,make_zeros=False)
        ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,ls='--',typ='uncond',make_zeros=True)#,zlev=0)

        ax.set_ylabel('Rain Rate (mm/hr)',fontsize=16)
        #ax.set_title('Precipitation Timeseries ')
        ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        #plt.tight_layout()
        #plt.savefig('{i}precip_timeseries_convstrat_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
        plt.savefig('{i}{e}_{v}_timeseries_convstrat.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,v=rdata.rr_name),dpi=400,bbox_inches='tight')
        plt.close()
 
        print('\nDone! Saved to '+config['image_dir'])
        print('Moving on.\n')
       

    ############################################################################

    ################################################################################
    ##Next let's make quantile (50,90,99) plots of the vertical velocity. This splits it by up and down, but you can turn split_updn == False
    if rdata.w_name is not None:

        if (config['vv_profiles'] | config['all1']):
 
            print('\nIN RUN_IPOLARRIS_NEW... creating vertical profile figure.')
            print('Plotting vertical profile for variable '+rdata.w_name+'...')
 
            outdir = config['image_dir']+'vertical_profile/'
            os.makedirs(outdir,exist_ok=True)
 
            fig,ax = plt.subplots(1,1,figsize=(12,8))
            ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=True)
            ax = plot_driver.plot_quartiles(rdata.data[rdata.w_name],0.9,0.5,0.99,rdata.data[rdata.z_name],ax,split_updn=False)
            ax.set_xlabel('Vertical Velocity (m/s)',fontsize=16)
            #ax.set_title('Vertical velocity profiles')
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
 #plt.tight_layout()
            #plt.savefig('{i}quantile_vvel_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
            plt.savefig('{i}{e}_{v}_vertprof.{p}'.format(p=config['ptype'],i=outdir,e=rdata.exper,v=rdata.w_name),dpi=400,bbox_inches='tight')
            plt.close()

            print('\nDone! Saved to '+config['image_dir'])
            print('Moving on.\n')
 
        if (config['percentiles_txt'] | config['all2']):

            print('\nIN RUN_IPOLARRIS_NEW... creating percentile text file.')
            print('Printing percentile data for variable '+rdata.w_name+'...')
 
            p99u,p90u,p50u,ht = rdata.percentile(wup=True)
            p99d,p90d,p50d,ht = rdata.percentile(wdown=True)
            p99a,p90a,p50a,ht = rdata.percentile(wdown=False)

            outdir = config['image_dir']+'txtfiles/'
            os.makedirs(outdir,exist_ok=True)
            file = open('{i}{e}_{v}_updown_percentiles.txt'.format(i=outdir,v=rdata.w_name,e=rdata.exper),'w') 
     
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

            print('\nDone! Saved to '+outdir)
            print('Moving on.\n')
 
    else:
        print("\nNo vertical velocity data.")
        print('Moving on.\n')
    
    ################################################################################

    ################################################################################
    
    if (config['vert_ref'] | config['all1']):
 
        print('\nIN RUN_IPOLARRIS_NEW... creating vertical profile figure.')
        print('Plotting vertical profile for variable '+rdata.dz_name+'...')
       
        outdir = config['image_dir']+'vertical_profile/'
        os.makedirs(outdir,exist_ok=True)
        
        ##Next let's make mean vertical profile of reflectivity
        fig,ax = plt.subplots(1,1,figsize=(12,8))
        ax = plot_driver.plot_verprof(rdata.data[rdata.dz_name],rdata.data[rdata.z_name],ax,split_updn=False,lab='dz',thresh=-50)
        #ax.set_title('Vertical profile of reflectivity')
        ax.set_xlabel('Reflectivity (dBZ)',fontsize=16)
        ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        #plt.tight_layout()
        #plt.savefig('{i}meanprofile_refl_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
        plt.savefig('{i}{e}_{v}_vertprof.{p}'.format(p=config['ptype'],i=outdir,e=rdata.exper,v=rdata.dz_name),dpi=400,bbox_inches='tight')
        plt.close()
 
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
  
    '''
    if config['refcfad']:

        print('\nIN RUN_IPOLARRIS_NEW... creating CFAD figure.')
        print('Plotting CFAD for variable '+rdata.dz_name+'...')
 
    ################################################################################
    ##Next let's make a reflectivity CFAD

#    cfaddat,vbins = plot_driver.cfad(rdata.data[rdata.dz_name],rdata,rdata.data[rdata.z_name],var=rdata.dz_name,nbins=40)
        cfaddat,vbins,r1ht = rdata.cfad(rdata.dz_name,ret_z=1,z_resolution=1.0,value_bins=np.arange(0,82,2),cscfad=False)

        fig,ax = plt.subplots(1,1,figsize=(12,8))
        ax = plot_driver.plot_cfad(fig,cfaddat, hts = r1ht,  vbins = vbins,ax=ax,cfad_on = 0,tspan = config['date'],maxval=20,cont=True,levels = True)

        ax.set_xlabel('Reflectivity (dBZ)',fontsize=16)
        #ax.set_title('{c} CFAD'.format(c=rdata.exper))
        ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.radar_name), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        #plt.tight_layout()
        #plt.savefig('{i}CFAD_refl_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
        plt.savefig('{i}{e}_{v}_CFAD.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,v=rdata.dz_name),dpi=400,bbox_inches='tight')
        plt.close()

        print('\nDone! Saved to '+config['image_dir'])
        print('Moving on.\n')
    '''

    #################################################################################

    print('\n########################################')
    print('############ Calling plot_driver.py ####')
    print('#########################################\n')
    time.sleep(3)

    plot_driver.make_single_pplots(rdata,config)

    print('\n#################################################')
    print('####### Returning to run_ipolarris_new.py #######')
    print('#################################################\n')

    print('\niPOLARRIS RUN COMPLETE FOR '+config['exper']+' '+config['sdatetime']+' - '+config['edatetime']+'\n')
