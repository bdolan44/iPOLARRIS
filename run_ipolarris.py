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
from plot_driver import make_single_pplots
from plot_driver import make_compare_pplots

#######################################################


# mc3e_wrf = get_data(exper='MC3E',radar_files=r'/Users/bdolan/scratch/POLARRIS_2/wrf_mc3e_files.txt')
twp_wrf = get_data(exper='TWPICE',radar_files=r'./wrf_twp_files.txt')
flags = {'cfad_4panel_flag':False,
         'hid_cfad_flag':False,
         'joint_flag':False,
         'cfad_individ_flag':False,
         'hid_prof':False,
         'all_cappi':True,
         'all_xsec':False,
         'up_width':False,
         'qr_cappi':False,
         'qr_rhi':False,}
#         
#         
 make_single_pplots(twp_wrf,[twp_wrf.date[0],twp_wrf.date[-1]],flags,'/gpfsm/dnb32/bcabell/GSDSU_MASTER_V4Beta/POLARRIS_images/',exp='TWPICE_WRF',extra=twp_wrf.mphys)
# make_single_pplots(mc3e_wrf,[mc3e_wrf.date[0],mc3e_wrf.date[-1]],flags,'/Users/bdolan/scratch/gitlab/iPOLARRIS/',exp='MC3E_WRF',extra=mc3e_wrf.mphys)
#         
# flags={'cfad_comp':False,
#         'hid_vert_comp':True,
#       'w_area_comp':False,}
# 
# make_compare_pplots(twp_wrf,mc3e_wrf,flags,'./',extra=mc3e_wrf.mphys)
# 

#######################################################

