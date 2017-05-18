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
from polarris_config import run_exper
from polarris_config import get_data
import warnings
warnings.filterwarnings('ignore')
import GeneralFunctions as GF


def foo(s1):
    return '{}'.format(s1.rstrip())
#######################################################
#######################################################
def get_time(time_parse,filename,dformat):
    tp = time_parse[0]
    te = time_parse[1]
#    print filename
    base = os.path.basename(filename)
    radar_name = base[:3]
    radcdate=np.str(base[tp:tp+te])
    date=datetime.datetime.strptime(radcdate,dformat)
    
    return date

#mc3e_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_mc3e.txt'
mc3e_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt'
mc3e_yp = 'wrf'
mc3e_exper = 'MC3E'
mc3e_mphys = None
mc3e_date = '20110523'
mc3e_time_parse=[11,19]
#mc3e_wdate_format = '%Y%m%d_%H%M%S'
mc3e_wdate_format = '%Y-%m-%d_%H-%M-%S'


#twpice_radar_files = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/iPOLARRIS/wrf_sbm_twpice.txt'
twpice_radar_files = r'/Users/bdolan/scratch/POLARRIS_2/wrf_twp_files.txt'
twpice_yp = 'wrf'
twpice_exper = 'TWPICE'
twpice_mphys = None
twpice_date = '2006123'
twpice_time_parse=[11,19]
twpice_wdate_format = '%Y-%m-%d_%H-%M-%S'
#twpice_wdate_format = '%Y%m%d_%H%M%S'



ptype = 'png'
#image_dir = r'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/'
image_dir = r'/Users/bdolan/scratch/iPOLARRIS_images/'

mc3e_dat = run_exper(mc3e_radar_files,mc3e_exper,mc3e_mphys,mc3e_date,mc3e_time_parse,mc3e_wdate_format,mc3e_yp)
twpice_dat = run_exper(twpice_radar_files,twpice_exper,twpice_mphys,twpice_date,twpice_time_parse,twpice_wdate_format,twpice_yp)


st_mc3e = 'MC3E {m} 21-02'.format(m=mc3e_mphys)#mc3e_dat['rconf'].sav_title()
st_twpice = 'TWPICE {m} 16-21'.format(m=twpice_mphys)#twpice_dat['rconf'].sav_title()
st_diff = 'MC3E-TWPICE_{m}'.format(m=mc3e_mphys)

mc3ecnt = np.shape(mc3e_dat['hts'])[0]
twpicecnt = np.shape(twpice_dat['hts'])[0]



#############1st dBZ##############


fig, ax = plt.subplots(1,3,figsize=(18,8))
axf = ax.flatten()

m_dzcfad_all = np.sum(mc3e_dat['dzcfad'],axis=0)/mc3ecnt
t_dzcfad_all = np.sum(twpice_dat['dzcfad'],axis=0)/twpicecnt

fig, ax = GF.cfad_plot('DZ',cfad = m_dzcfad_all, hts = mc3e_dat['hts'][0],  bins = mc3e_dat['dzbins'],ax=axf[0],cfad_on = 0,rconf = mc3e_dat['rconf'],tspan = mc3e_dat['time'],maxval=20,cont=True,levels = True)

fig, ax = GF.cfad_plot('DZ',cfad = t_dzcfad_all, hts = twpice_dat['hts'][0],  bins = twpice_dat['dzbins'],ax=axf[1],cfad_on = 0,rconf = twpice_dat['rconf'],tspan = twpice_dat['time'],maxval=20,cont=True,levels = True)


diff_cfad = m_dzcfad_all - t_dzcfad_all
cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
#print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))

#levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
levels=np.arange(-2,2.1,0.1)
cb=axf[2].contourf(mc3e_dat['dzbins'][:-1],mc3e_dat['hts'][0],cfad_ma,levels,cmap='bwr',extend='both')

#    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
plt.colorbar(cb,ax=axf[2])
axf[2].set_ylabel('Height (km MSL)',fontsize=18)
axf[2].set_xlabel(mc3e_dat['rconf'].names['DZ'],fontsize = 18)

axf[2].set_title('{d} {v}'.format(d=st_diff,v=mc3e_dat['rconf'].dz_name))

plt.tight_layout()

plt.savefig('{id}CFADDZ{s}.{t}'.format(id=image_dir,s=st_diff,t=ptype),dpi=200)

#############Now DR##############

fig, ax = plt.subplots(1,3,figsize=(18,8))
axf = ax.flatten()

m_drcfad_all = np.sum(mc3e_dat['drcfad'],axis=0)/mc3ecnt
t_drcfad_all = np.sum(twpice_dat['drcfad'],axis=0)/twpicecnt

fig, ax = GF.cfad_plot('DR',cfad = m_drcfad_all, hts = mc3e_dat['hts'][0],  bins = mc3e_dat['drbins'],ax=axf[0],cfad_on = 0,rconf = mc3e_dat['rconf'],tspan = mc3e_dat['time'],maxval=20,cont=True,levels = True)

fig, ax = GF.cfad_plot('DR',cfad = t_drcfad_all, hts = twpice_dat['hts'][0],  bins = twpice_dat['drbins'],ax=axf[1],cfad_on = 0,rconf = twpice_dat['rconf'],tspan = twpice_dat['time'],maxval=20,cont=True,levels = True)



diff_cfad = m_drcfad_all - t_drcfad_all
cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
#print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))

#levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
levels=np.arange(-8,8.1,0.1)
cb=axf[2].contourf(mc3e_dat['drbins'][:-1],mc3e_dat['hts'][0],cfad_ma,levels,cmap='bwr',extend='both')

#    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
plt.colorbar(cb,ax=axf[2])
axf[2].set_ylabel('Height (km MSL)',fontsize=18)
axf[2].set_xlabel(mc3e_dat['rconf'].names['DR'],fontsize = 18)

axf[2].set_title('{d} {v}'.format(d=st_diff,v=mc3e_dat['rconf'].dr_name))

plt.tight_layout()

plt.savefig('{id}CFADDR{s}.{t}'.format(id=image_dir,s=st_diff,t=ptype),dpi=200)

#############Now W##############

fig, ax = plt.subplots(1,3,figsize=(18,8))
axf = ax.flatten()

m_kdcfad_all = np.sum(mc3e_dat['kdcfad'],axis=0)/mc3ecnt
t_kdcfad_all = np.sum(twpice_dat['kdcfad'],axis=0)twpicecnt
fig, ax = GF.cfad_plot('KD',cfad = m_kdcfad_all, hts = mc3e_dat['hts'][0],  bins = mc3e_dat['kdbins'],ax=axf[0],cfad_on = 0,rconf = mc3e_dat['rconf'],tspan = mc3e_dat['time'],maxval=20,cont=True,levels = True)

fig, ax = GF.cfad_plot('KD',cfad = t_kdcfad_all, hts = twpice_dat['hts'][0],  bins = twpice_dat['kdbins'],ax=axf[1],cfad_on = 0,rconf = twpice_dat['rconf'],tspan = twpice_dat['time'],maxval=20,cont=True,levels = True)


diff_cfad = m_kdcfad_all - t_kdcfad_all
cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
#print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))

#levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
levels=np.arange(-8,8.1,0.1)
cb=axf[2].contourf(mc3e_dat['kdbins'][:-1],mc3e_dat['hts'][0],cfad_ma,levels,cmap='bwr',extend='both')

#    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
plt.colorbar(cb,ax=axf[2])
axf[2].set_ylabel('Height (km MSL)',fontsize=18)
axf[2].set_xlabel(mc3e_dat['rconf'].names['KD'],fontsize = 18)

axf[2].set_title('{d} {v}'.format(d=st_diff,v=mc3e_dat['rconf'].kd_name))

plt.tight_layout()

plt.savefig('{id}CFADKD{s}.{t}'.format(id=image_dir,s=st_diff,t=ptype),dpi=200)

#############Now W##############

fig, ax = plt.subplots(1,3,figsize=(18,8))
axf = ax.flatten()

m_wcfad_all = np.sum(mc3e_dat['wcfad'],axis=0)/mc3ecnt
t_wcfad_all = np.sum(twpice_dat['wcfad'],axis=0)/twpicecnt
fig, ax = GF.cfad_plot('KD',cfad = m_wcfad_all, hts = mc3e_dat['hts'][0],  bins = mc3e_dat['wbins'],ax=axf[0],cfad_on = 0,rconf = mc3e_dat['rconf'],tspan = mc3e_dat['time'],maxval=20,cont=True,levels = True)

fig, ax = GF.cfad_plot('KD',cfad = t_wcfad_all, hts = twpice_dat['hts'][0],  bins = twpice_dat['wbins'],ax=axf[1],cfad_on = 0,rconf = twpice_dat['rconf'],tspan = twpice_dat['time'],maxval=20,cont=True,levels = True)


diff_cfad = m_wcfad_all - t_wcfad_all
cfad_ma = np.ma.masked_where(diff_cfad == 0, diff_cfad)
#print np.shape(bins[0:-1]), np.shape(radz), np.shape(np.transpose(cfad_ma))

#levels=[1.,5.,10.,20.,25.,30.,35.,40.,45.,50.,55.,60.]
levels=np.arange(-8,8.1,0.1)
cb=axf[2].contourf(mc3e_dat['wbins'][:-1],mc3e_dat['hts'][0],cfad_ma,levels,cmap='bwr',extend='both')

#    pc = plt.contour([bins[0:-1]],radz, np.transpose(cfad_ma))#, vmin = 0.02, vmax = 15, cmap = cmap)
plt.colorbar(cb,ax=axf[2])
axf[2].set_ylabel('Height (km MSL)',fontsize=18)
axf[2].set_xlabel(mc3e_dat['rconf'].names['Wvar'],fontsize = 18)

axf[2].set_title('{d} {v}'.format(d=st_diff,v=mc3e_dat['rconf'].dz_name))

plt.tight_layout()

plt.savefig('{id}CFADW{s}.{t}'.format(id=image_dir,s=st_diff,t=ptype),dpi=200)

