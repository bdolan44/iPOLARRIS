#===================================================
#================ RUN_IPOLARRIS.PY =================
#===================================================

import sys
import time
'''
# WARNING TO USERS to activate conda environment #
print("\n####################################")
print("### Welcome, user, to iPOLARRIS! ###")
print("####################################\n")

print("WARNING: Before proceeding, ensure that you have: \n\n (a) installed the Anaconda package manager. Latest versions and instructions can be found here: \n https://conda.io/projects/conda/en/latest/user-guide/install/index.html \n\n (b) installed the Conda environment required to run iPOLARRIS with the command `conda env create -f env.yml` \n\n (c) activated the new environment with the command `conda activate pol` \n\n (d) Put /usr/bin ahead of your 'default' executable directories (i.e. before /opt/local/bin) in $PATH, if it is not already. It does not need to be ahead of your custom executable directories (i.e. ~/../anaconda3/bin) \n\n (e) Run: `f2py -c calc_kdp_ray_fir.f -m calc_kdp_ray_fir`. This will allow you to use the Fortran compiler in /usr/bin to convert your .f file into a readable .so file for your MAC or Linux OS. \n")

print("If you have NOT performed the required setup above, click x and Enter to exit. Otherwise, press any other key. \n")

usersays=input()
if usersays.lower().startswith('x'): 
    print('\nExiting gracefully.\n')
    import sys
    sys.exit()
else: 
    print('\niPOLARRIS INITIATING... If Conda env activated, no import errors...')
    import time
    time.sleep(3)
'''
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
from polarris_driver import polarris_driver
import plot_driver
import RadarData
import RadarConfig
from skewPy import SkewT

#--------------- Main Program ----------------

time.sleep(3)

print('\nSUCCESS! Requisite packages loaded.')

if len(sys.argv) > 2:
    print('\n***Entering SIMULATION MODE: you are about to compare radar observations with simulated radar observables created from wrfout files by POLARRIS-f!***')
else:
    print('\n***Entering OBSERVATION MODE: you are about to analyze radar observations recorded by a station!***')

time.sleep(3)

print('\n#############################################')
print('########## Starting run_ipolarris.py ########')
print('#############################################')

configfile = sys.argv[1:] # Feed config file name as arg
#print sys.argv[1:]

print('\n##########################################################')
print('############ Calling polarris_driver.py to read in obs ###')
print('##########################################################')

time.sleep(3)

rdata, config, config['uname'], config['vname'], config['wname'] = polarris_driver(configfile)

print('\n#################################################')
print('########## Returning to run_ipolarris.py ########')
print('#################################################')

#print(,'EXTRA 1 is')

# If a second argument is passed for WRF config file, produce a bunch of comparison plots!
# More comments in this section TBD!
if sys.argv[2:]:
   
    configfile1 = sys.argv[2:]

    print('\n###############################################################')
    print('############ Calling polarris_driver.py to read in sim data ###')
    print('###############################################################')
    time.sleep(3)

    rdata2, config2, config2['uname'], config2['vname'], config2['wname'] = polarris_driver(configfile1)

    print('\n#################################################')
    print('########## Returning to run_ipolarris.py ########')
    print('#################################################')

    if (config2['cfad_compare'] | config2['all3']):
 
        print('\nIN RUN_IPOLARRIS... creating CFAD COMPARISON figures.')
        outdir = config['image_dir']+'cfad_diff_individ/'
        os.makedirs(outdir,exist_ok=True)
 
        zmax = config['zmax']
        st = rdata.date[0].strftime('%Y%m%d_%H%M%S')
        en = rdata.date[-1].strftime('%Y%m%d_%H%M%S')

        if st.startswith(en): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+' UTC'
        elif st[0:8].startswith(en[0:8]): dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+'-'+en[9:11]+':'+en[11:13]+' UTC'
        else: dtlab = st[0:4]+'-'+st[4:6]+'-'+st[6:8]+' '+st[9:11]+':'+st[11:13]+'\n- '+en[0:4]+'-'+en[4:6]+'-'+en[6:8]+' '+en[9:11]+':'+en[11:13]+' UTC'

        for i,v in enumerate(eval(config['cfad_compare_vars'])):
            
            if v is None:
                continue 
            else:
                
                if v.startswith('HID'):
 
                    print(v)

                    fig, ax = plt.subplots(1,2,figsize=(14,8),gridspec_kw={'wspace': 0.08, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.}) 
                    if not isinstance(ax, np.ndarray) or not isinstance(ax, list): 
                        ax = np.array([ax])
                    axf = ax.flatten()

                    if not zmax == '':
                        rdata.plot_hid_cdf(ax=axf[0],cbar=False,z_resolution=config['z_resolution'],zmax=zmax)
                        rdata2.plot_hid_cdf(ax=axf[1],ylab=False,cbar=False,z_resolution=config['z_resolution'],zmax=zmax)
                    else:
                        rdata.plot_hid_cdf(ax=axf[0],cbar=False,z_resolution=config['z_resolution'])
                        rdata2.plot_hid_cdf(ax=axf[1],ylab=False,cbar=False,z_resolution=config['z_resolution'])
                    
                    lur,bur,wur,hur = axf[0].get_position().bounds
                    lur2,bur2,wur2,hur2 = axf[1].get_position().bounds
                    cbar_ax_dims = [lur,bur-0.15,lur2+wur2,0.05]
                    rdata.HID_barplot_colorbar(fig,cbar_ax_dims,orientation='horizontal',names='longnames')

                    axf[0].text(0,1,'{e} {r}'.format(e=rdata.exper,r=rdata.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=axf[0].transAxes)
                    axf[1].text(0,1,'{e} {r}'.format(e=rdata2.exper,r=rdata2.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=axf[1].transAxes)
                    axf[1].text(1,1, dtlab, horizontalalignment='right', verticalalignment='bottom', size=18, color='k', zorder=10, weight='bold', transform=axf[1].transAxes) # (a) Top-left
                   
                    if config['ptype'].startswith('mp4'):
                        plt.savefig('{d}{p}_HID_CFAD_{t1}-{t2}.png'.format(d=outdir,p=rdata.exper,t1=st,t2=en),dpi=400,bbox_inches='tight')
                    else: 
                        plt.savefig('{d}{p}_HID_CFAD_{t1}-{t2}.{t}'.format(d=outdir,p=rdata.exper,t=config['ptype'],t1=st,t2=en),dpi=400,bbox_inches='tight')

                    plt.close()
                
                else:

                    if not rdata.cfbins[config[v]] == '' and config[v] in rdata.data.variables.keys():
                        
                        print(v)

                        fig, ax = plt.subplots(1,3,figsize=(16,8),gridspec_kw={'wspace': 0.1, 'top': 1., 'bottom': 0., 'left': 0., 'right': 1.}) 
                        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): 
                            ax = np.array([ax])
                        axf = ax.flatten()

                        if not zmax == '':
                            ocfad, hts, pc, fig0, ax0 = rdata.cfad_plot(config[v],ax=axf[0],cbar=False,bins=rdata.cfbins[config[v]],z_resolution=config['z_resolution'],levels=1,zmax=zmax)
                            scfad, hts2, pc2, fig1, ax1 = rdata2.cfad_plot(config2[v],ax=axf[1],ylab=False,cbar=False,bins=rdata2.cfbins[config2[v]],z_resolution=config['z_resolution'],levels=1,zmax=zmax)
                            dcfad, hts3, pc3, fig2, ax2 = rdata.cfad_plot(config[v],cfad=ocfad-scfad,hts=hts,ax=axf[2],ylab=False,cbar=False,bins=rdata.cfbins[config[v]],z_resolution=config['z_resolution'],levels=1,zmax=zmax,diff=1)
                        else:
                            ocfad, hts, pc, fig0, ax0 = rdata.cfad_plot(config[v],ax=axf[0],cbar=False,bins=rdata.cfbins[config[v]],z_resolution=config['z_resolution'],levels=1)
                            scfad, hts2, pc2, fig1, ax1 = rdata2.cfad_plot(config2[v],ax=axf[1],ylab=False,cbar=False,bins=rdata2.cfbins[config2[v]],z_resolution=config['z_resolution'],levels=1)
                            dcfad, hts3, pc3, fig2, ax2 = rdata.cfad_plot(config[v],cfad=ocfad-scfad,hts=hts,ax=axf[2],ylab=False,cbar=False,bins=rdata.cfbins[config[v]],z_resolution=config['z_resolution'],levels=1,diff=1)

                        lur,bur,wur,hur = axf[0].get_position().bounds
                        lur2,bur2,wur2,hur2 = axf[1].get_position().bounds
                        cbar_ax_dims = [lur,bur-0.13,lur2+wur2,0.03]
                        cbar_ax = fig.add_axes(cbar_ax_dims)
                        cbt = plt.colorbar(pc,cax=cbar_ax,orientation='horizontal')
                        cbt.ax.tick_params(labelsize=16)
                        cbt.set_ticks(rdata.cfad_levs)
                        cbt.set_label('Frequency (%)', fontsize=16, labelpad=10)

                        lur3,bur3,wur3,hur3 = axf[2].get_position().bounds
                        cbar_ax_dims3 = [lur3,bur3-0.13,wur3,0.03]
                        cbar_ax3 = fig.add_axes(cbar_ax_dims3)
                        cbt3 = plt.colorbar(pc3,cax=cbar_ax3,orientation='horizontal')
                        cbt3.ax.tick_params(labelsize=16)
                        cbt.set_ticks(rdata.cfad_levs)
                        cbt3.set_label('Frequency Difference (%)', fontsize=16, labelpad=10)

                        axf[0].text(0,1,'{e} {r}'.format(e=rdata.exper,r=rdata.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=axf[0].transAxes)
                        axf[1].text(0,1,'{e} {r}'.format(e=rdata2.exper,r=rdata2.band+'-band'),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=axf[1].transAxes)
                        axf[2].text(0,1,'({e1} - {e2})'.format(e1=rdata.exper,e2=rdata2.exper),horizontalalignment='left',verticalalignment='bottom',size=18,color='k',zorder=10,weight='bold',transform=axf[2].transAxes)
                        axf[2].text(0.99,0.99, dtlab, horizontalalignment='right', verticalalignment='top', size=18, color='k', zorder=10, weight='bold', transform=axf[2].transAxes, bbox=dict(facecolor='w', edgecolor='none', pad=0.0)) # (a) Top-left
                      
                        if config['ptype'].startswith('mp4'):
                            plt.savefig('{d}{p}_{v}_CFAD_{t1}-{t2}.png'.format(d=outdir,p=rdata.exper,v=rdata.names_uc[config[v]],t1=st,t2=en),dpi=400,bbox_inches='tight')
                        else: 
                            plt.savefig('{d}{p}_{v}_CFAD_{t1}-{t2}.{t}'.format(d=outdir,p=rdata.exper,v=rdata.names_uc[config[v]],t=config['ptype'],t1=st,t2=en),dpi=400,bbox_inches='tight')
        
                        plt.close()

                    else:
                        
                        continue

                print('\nDone! Saved to '+outdir)
                print('Moving on.\n')
  
        '''
        fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,config,bins=config['dzbins'],xlab=rdata.longnames[rdata.dz_name]+' '+rdata.units[rdata.dz_name],cscfad=False,xlim=[min(rdata.cfbins[rdata.dz_name]),max(rdata.cfbins[rdata.dz_name])],ylim=config['zmax'],nor=10)
        ax[0].set_title(rdata.exper+' '+rdata.exper,fontsize=16,fontweight='bold')
        ax[0].set_ylabel('Height (km MSL)',fontsize=16)
        ax[1].set_title(rdata2.mphys.upper(),fontsize=16,fontweight='bold')
        ax[3].set_title('({e} - {v})'.format(e=rdata.exper,v=rdata2.mphys.upper()),fontsize=16,fontweight='bold')    
        plt.savefig('{i}{e}_{m}_{v}_CFAD_diff.{p}'.format(p=config2['ptype'],i=outdir,e=rdata2.exper,m=rdata2.mphys.upper(),v=rdata.dz_name),dpi=400,bbox_inches='tight')
        plt.close(fig)

        print('\nDone! Saved to '+outdir)
        print('Moving on.')
        print('\nPlotting composites by time for variable '+rdata.zdr_name+'...')
        
        fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,config,bins=config['drbins'],xlab=rdata.longnames[rdata.zdr_name]+' '+rdata.units[rdata.zdr_name],cscfad=False,xlim=[min(config['drbins']),max(config['drbins'])+1],ylim=config['zlim'],nor=3)
        ax[0].set_title(rdata.exper+' '+rdata.exper,fontsize=16,fontweight='bold')
        ax[0].set_ylabel('Height (km MSL)',fontsize=16)
        ax[1].set_title(rdata2.mphys.upper(),fontsize=16,fontweight='bold')
        ax[3].set_title('({e} - {v})'.format(e=rdata.exper,v=rdata2.mphys.upper()),fontsize=16,fontweight='bold')    
        plt.savefig('{i}{e}_{m}_{v}_CFAD_diff.{p}'.format(p=config2['ptype'],i=outdir,e=rdata2.exper,m=rdata2.mphys.upper(),v=rdata.zdr_name),dpi=400,bbox_inches='tight')
        plt.close(fig)

        print('\nDone! Saved to '+outdir)
        print('Moving on.')
        print('\nPlotting composites by time for variable '+rdata.kdp_name+'...')

        fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,config,bins=config['kdbins'],xlab=rdata.longnames[rdata.kdp_name]+' '+rdata.units[rdata.kdp_name],cscfad=False,xlim=[min(config['kdbins']),max(config['kdbins'])+1],ylim=config['zlim'],nor=3)
        ax[0].set_title(rdata.exper+' '+rdata.exper,fontsize=16,fontweight='bold')
        ax[0].set_ylabel('Height (km MSL)',fontsize=16)
        ax[1].set_title(rdata2.mphys.upper(),fontsize=16,fontweight='bold')
        ax[3].set_title('({e} - {v})'.format(e=rdata.exper,v=rdata2.mphys.upper()),fontsize=16,fontweight='bold')    
        plt.savefig('{i}{e}_{m}_{v}_CFAD_diff.{p}'.format(p=config2['ptype'],i=outdir,e=rdata2.exper,m=rdata2.mphys.upper(),v=rdata.kdp_name),dpi=400,bbox_inches='tight')
        plt.close(fig)

        print('\nDone! Saved to '+outdir)
        print('Moving on.')
        print('\nPlotting composites by time for variable '+rdata.hid_name+'...')
        #fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad=False)
        fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad=None)
        ax[0].set_title(rdata.exper+' '+rdata.exper,fontsize=16,fontweight='bold')
        ax[1].set_title(rdata2.mphys.upper(),fontsize=16,fontweight='bold')
        plt.setp(ax, ylim=[0,10])
        plt.savefig('{i}{e}_{m}_{v}_CFAD_diff.{p}'.format(p=config2['ptype'],i=outdir,e=rdata2.exper,m=rdata2.mphys.upper(),v=rdata.hid_name),dpi=400,bbox_inches='tight')
        plt.close(fig)

        print('\nDone! Saved to '+outdir)
        print('Moving on.')

    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.dz_name,rdata2.dz_name,'Reflectivity',config,config2,bins=np.arange(0,82,2),savefig=True,cscfad='convective')
        
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.zdr_name,rdata2.zdr_name,'Z$_{dr}$',config,config2,bins=np.arange(-2,8,0.2),savefig=True,cscfad='convective')
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.kdp_name,rdata2.kdp_name,'K$_{dp}$',config,config2,bins=np.arange(-2,6,0.2),savefig=True,cscfad='convective')
    
    fig,ax = plot_driver.plot_difference_cfad(rdata,rdata2,rdata.w_name,rdata2.w_name,'Vertical Velocity',config,config2,bins=np.arange(-20,21,1),savefig=True,cscfad='convective')

    fig,ax = plot_driver.plot_hid_comparison_cfad(rdata,rdata2,config=config,cscfad='convective',savefig=True)
    '''
    print('\niPOLARRIS RUN COMPLETE FOR '+config2['mphys'].upper()+' '+config['sdatetime']+' - '+config['edatetime']+'\n')

################################################################################
################## Now you can just start plotting! ############################
################################################################################

### To see the variables that are available to plot, type:

#rdata.data.keys()

else:
    
    #############################################################################
   
    if (config['compo_ref'] | config['all1']):
    
        print('\nIN RUN_IPOLARRIS_NEW... creating COMPOSITE figures.')
        print('\nPlotting composites by time for variable '+rdata.dz_name+'...')
            
        outdir = config['image_dir']+'composite_'+rdata.dz_name+'/'
        os.makedirs(outdir,exist_ok=True)
      
        for i,rtimematch in enumerate(np.array(rdata.date)):

            fig, ax = rdata.plot_composite(rdata.dz_name,i,statpt=True)
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
            ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=rtimematch), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
           
            if not config['ptype'].startswith('mp4'):
                plt.savefig('{i}{e}_{v}_{d:%Y%m%d_%H%M%S}.{p}'.format(p=config['ptype'],e=rdata.exper,i=outdir,d=rtimematch,v=rdata.dz_name),dpi=400,bbox_inches='tight')
            else: 
                if len(rdata.date) < 6:
                    plt.savefig('{i}{e}_{v}_{d:%Y%m%d_%H%M%S}.png'.format(e=rdata.exper,i=outdir,d=rtimematch,v=rdata.dz_name),dpi=400,bbox_inches='tight')
                else:
                    plt.savefig(outdir+'/fig'+str(i).zfill(3)+'.png',dpi=400,bbox_inches='tight')

            plt.close()
            print(rtimematch)
        
        if config['ptype'].startswith('mp4') and len(rdata.date) >= 6:
            
            st = rdata.date[0].strftime('%Y%m%d_%H%M%S')
            en = rdata.date[-1].strftime('%Y%m%d_%H%M%S')
           
            os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdata.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_{v}_{t1}-{t2}.mp4'.format(p=config['ptype'],e=rdata.exper,i=outdir,v=rdata.dz_name,t1=st,t2=en))

        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
    
    #############################################################################

    if (config['cappi_rr'] | config['all1']):

        print('\nIN RUN_IPOLARRIS_NEW... creating CAPPI figures.')
        print('Plotting CAPPIs for all heights by time for variable '+rdata.rr_name+'...')
 
        outdir = config['image_dir']+'cappi_'+rdata.rr_name+'/'
        os.makedirs(outdir,exist_ok=True)

        if not config['z'] == '': zspan = list(eval(str([config['z']])))
        else: zspan = rdata.data[rdata.z_name].values

        for z in zspan:

            print('\nz = '+str(z))

            for i,rtimematch in enumerate(np.array(rdata.date)):
                
                print(rtimematch)

                dummy, ax = rdata.cappi(rdata.rr_name,z=z,xlim=config['xlim'],ylim=config['ylim'],ts=rtimematch,latlon=config['latlon'],statpt=True,xlab=True,ylab=True,dattype=config['type'])

                ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                ax.text(1, 1, '{d:%Y-%m-%d %H:%M:%S} UTC'.format(d=rtimematch), horizontalalignment='right', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
                ax.text(0.99, 0.99, 'z = {a} km'.format(a=z), horizontalalignment='right',verticalalignment='top', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes, bbox=dict(facecolor='w', edgecolor='none', pad=0.0))

                if not config['ptype'].startswith('mp4'):
                    plt.savefig('{i}{e}_{v}_cappi_{d:%Y%m%d_%H%M%S}_{h}.png'.format(i=outdir,e=rdata.exper,h=z,v=rdata.rr_name,d=rtimematch),dpi=400,bbox_inches='tight')
                else: 
                    if len(rdata.date) < 6:
                        plt.savefig('{i}{e}_{v}_cappi_{d:%Y%m%d_%H%M%S}_{h}.png'.format(i=outdir,e=rdata.exper,h=z,v=rdata.rr_name,d=rtimematch),dpi=400,bbox_inches='tight')
                    else:
                        plt.savefig(outdir+'/fig'+str(i).zfill(3)+'.png',dpi=400,bbox_inches='tight')

                plt.close()

            if config['ptype'].startswith('mp4') and len(rdata.date) >= 6:

                st = rdata.date[0].strftime('%Y%m%d_%H%M%S')
                en = rdata.date[-1].strftime('%Y%m%d_%H%M%S')

                os.system('ffmpeg -nostdin -y -r 1 -i '+outdir+'/fig%03d.png -c:v libx264 -r '+str(len(np.array(rdata.date)))+' -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" '+'{i}{e}_{v}_{t1}-{t2}_{h}.mp4'.format(p=config['ptype'],e=rdata.exper,i=outdir,v=rdata.rr_name,t1=st,t2=en,h=z))
            
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')

    #############################################################################
    
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
        
        ##First make a timeseries of rain rate, unconditional and conditional. This puts strat, conv, and total on the same plot but you can split the out by putting cs==False.
        ## The conditional rain rate is achieved by sending threshold = 0.
        fig,ax = plt.subplots(1,1,figsize=(12,8))
        ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,make_zeros=False)
        ax = plot_driver.plot_timeseries(rdata.data[rdata.rr_name],rdata.date,ax,cs=True,rdata=rdata,thresh=0,zlev=1,ls='--',typ='uncond',make_zeros=True)#,zlev=0)

        ax.set_ylabel('Rain Rate (mm/hr)',fontsize=16)
        #ax.set_title('Precipitation Timeseries ')
        ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.exper), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        #plt.tight_layout()
        #plt.savefig('{i}precip_timeseries_convstrat_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
        plt.savefig('{i}{e}_{v}_timeseries_convstrat.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,v=rdata.rr_name),dpi=400,bbox_inches='tight')
        plt.close()
 
        print('\nDone! Saved to '+config['image_dir'])
        print('Moving on.\n')
       

    ############################################################################

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
            ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.exper), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
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
    
    #############################################################################
 
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
        ax.text(0, 1, '{e} {r}'.format(e=rdata.exper,r=rdata.band+'-band'), horizontalalignment='left', verticalalignment='bottom', size=16, color='k', zorder=10, weight='bold', transform=ax.transAxes) # (a) Top-left
        #plt.tight_layout()
        #plt.savefig('{i}meanprofile_refl_{e}_{m}.{p}'.format(p=config['ptype'],i=config['image_dir'],e=rdata.exper,m=rdata.mphys),dpi=400,bbox_inches='tight')
        plt.savefig('{i}{e}_{v}_vertprof.{p}'.format(p=config['ptype'],i=outdir,e=rdata.exper,v=rdata.dz_name),dpi=400,bbox_inches='tight')
        plt.close()
 
        print('\nDone! Saved to '+outdir)
        print('Moving on.\n')
  
    #############################################################################

    print('\n########################################')
    print('############ Calling plot_driver.py ####')
    print('#########################################\n')
    time.sleep(3)

    plot_driver.make_single_pplots(rdata,config)

    print('\n#################################################')
    print('####### Returning to run_ipolarris_new.py #######')
    print('#################################################\n')

    print('\niPOLARRIS RUN COMPLETE FOR '+config['exper']+' '+config['sdatetime']+' - '+config['edatetime']+'\n')
