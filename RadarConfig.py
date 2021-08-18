# this houses the RadarConfig object which is basically just carrying
# along some options about names, colorbars, plotting limits, etc.

from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
import datetime as dt


class RadarConfig(object):
    
    def __init__(self, dz='DZ', zdr='DR', kdp='KD', ldr='LH', rho='RH', hid = 'HID',conv='Con',
            temp='T', x='x', y='y', z='z', u='U', v='V',rr='RR', w='Wvar',vr='VR',mphys='None',exper = 'Case',
            band = 'C',lat_0 = 0,lon_0=90.0,lat_r=None,lon_r=None,lat=None,lon=None,tm = None,radar_name = None,
            color_blind = False):
        # ******** first the polarimetric stuff *************
        self.dz_name = dz
        self.zdr_name = zdr
        self.kdp_name = kdp
        self.ldr_name = ldr
        self.rho_name = rho
        self.rr_name = rr
        if self.rr_name == None:
            self.rr_name = 'RR'
        self.temp_name = temp
        self.hid_name = hid
        self.vr_name = vr
        self.conv_name=conv
        self.cs_name = 'CSS'

        # ********** next the cartesian coordinates **********
        self.x_name = x
        self.y_name = y
        self.z_name = z
        # ********** dual doppler names **********************
        self.u_name = u
        self.v_name = v
        self.w_name = w

        self.lat_name=lat
        self.lon_name=lon
        print('in config, lat name is ',self.lat_name)
        
        ########set up some experiment parameters ############
        self.mphys=mphys
        self.exper=exper
        self.radar_lat = lat_r
        self.radar_lon = lon_r
        self.lat_0 = lat_0
        self.lon_0 = lon_0
        
        self.band = band
        self.expr =exper
        #print tm
        self.date = tm
        self.radar_name = radar_name



        self.species = np.array(['DZ','RN','CR','AG','WS','VI','LDG','HDG','HA','BD'])
        #self.hid_colors = ['White','LightBlue','MediumBlue','Darkorange','LightPink','Cyan','DarkGray',\
        #    'Lime','Yellow','Red','Fuchsia']
        self.hid_colors = ['LightBlue','MediumBlue','Darkorange','LightPink','Cyan','DarkGray',\
            'Lime','Yellow','Red','Fuchsia']       
        self.pol_vars = np.array([self.dz_name, self.zdr_name, self.kdp_name, self.ldr_name, self.rho_name, self.hid_name])

        self.cs_colors = ['#FFFFFF', 'DodgerBlue', 'Red', 'Khaki']
        self.cs_labels = ['', 'Strat', 'Conv', 'Mixed']

        self.set_dbz_colorbar(color_blind=color_blind)
        self.set_hid_colorbar()
        self.set_cs_colorbar()

        # Now just set some defaults
        self.lims = {dz: [0,80], zdr: [-1, 3], kdp: [-0.5, 3], ldr: [-35, -20], rho: [0.95, 1.00], hid: [0, len(self.species)+1],w:[-25,25],vr:[-25,25],self.cs_name:[0,4],self.rr_name:[0.01,150]}
        self.delta = {dz: 10, zdr: 1, kdp: 1, ldr: 5, rho: 0.005, hid: 1,w:5,vr:5,self.cs_name:1,self.rr_name:10}
        self.units = {dz: '(dBZ)', zdr: '(dB)', kdp: '($^{\circ}$/km)', ldr: '(dB)', rho: '', hid: '',w:'(m s$^{-1}$)',vr:'(m s$^{-1}$)',self.cs_name:'',self.rr_name:'(mm hr$^{-1}$)'}
        self.names = {dz: 'Z', zdr: 'Z$_{DR}$', kdp: 'K$_{dp}$', ldr: 'LDR', rho: r'$\rho_{hv}$', hid: '',w:'',vr:'V$_r$',self.cs_name:'',self.rr_name:'RR'}
        self.longnames = {dz: 'Reflectivity', zdr: 'Differntial reflectivity', kdp: 'Specific differential phase',\
                ldr: 'Linear depolarization ratio', rho: 'Correlation coefficient', hid: 'Hydrometeor identification',w:'Vertical Velocity',vr:'Radial Velocity',\
                self.cs_name: 'Convective/Stratiform',self.rr_name:'Rain Rate'}
        self.cmaps = {dz: self.temp_cmap, zdr: plt.cm.Spectral_r, kdp: plt.cm.gist_heat_r, \
                ldr: plt.cm.gist_rainbow_r, rho: plt.cm.jet, hid: self.hid_cmap,w:plt.cm.seismic,vr:plt.cm.bwr,self.cs_name: self.cs_cmap,self.rr_name:plt.cm.Spectral_r}
        self.ticklabels = {dz: np.arange(0, 90, 10), zdr: np.arange(-1, 4, 1), kdp: np.arange(-0.5, 4.5, 1), 
                ldr: np.arange(-35, -15, 5), rho: np.arange(0.95, 1.01, 0.005), hid: np.append('', self.species),w:np.arange(-25,30.0,5.0),vr:np.arange(-25,30.0,5.0),
                self.cs_name: self.cs_labels,self.rr_name:[0.1,1,10,30,50,70,100,130,150]}
#############################################################################################################

    def print_date(self,tm=None, fmt='%Y-%m-%d %H:%M:%S %Z'):
#        print tm
        if tm is not None:
            #print tm
            if len(tm) > 1:
                date1 = tm[0]
                date2 = tm[-1]
                tms = dt.datetime.strftime(date1,fmt)
                tme = dt.datetime.strftime(date2,fmt)
                tmf = '{s}-{e}'.format(s=tms,e=tme)
            else:
                tmf = dt.datetime.strftime(tm[0],fmt)
        else:
            tmf = dt.datetime.strftime(np.array(self.date)[0],fmt)
        #print tmf
        return tmf

#############################################################################################################
    def print_title(self,tm = None):
        extra = ''
        if self.exper is not None:
            extra = '{e} {x}'.format(e=extra,x = self.exper)
        if self.radar_name is not None:
            extra = '{e} {r}'.format(e=extra,r=self.radar_name)
        if self.mphys is not None:
            extra = '{e} {m}'.format(e=extra,m=self.mphys)
        
        extra = '{e}, {t}'.format(e=extra, t=self.print_date(tm))
        return extra
#############################################################################################################

    def sav_title(self,tm = None):
        extra = ''
        if self.exper is not None:
            extra = '{e}_{x}'.format(e=extra,x = self.exper)
#        if self.radar_name is not None:
#            extra = '{e}_{r}'.format(e=extra,r=self.radar_name)
        if self.mphys is not None:
            extra = '{e}_{m}'.format(e=extra,m=self.mphys)
        
#        extra = '{e}_{t}'.format(e=extra, t=self.print_date(tm,fmt = '%Y%m%d%_H%M%S'))
        return extra
#############################################################################################################

    def set_dbz_colorbar(self, color_list=None, color_blind=False):
        if color_list is None:
            # just use the default here
            if color_blind is not True:
                radarcbar = ['PeachPuff','Aqua','DodgerBlue','Blue','Lime', \
                    'LimeGreen','Green','Yellow','Orange','DarkOrange','Red', \
                    'Crimson','Fuchsia','Purple','Indigo','MidnightBlue'] 
            #radarcbar = ['PeachPuff','Aqua','DodgerBlue','MediumBlue','Lime', \
            #    'LimeGreen','Green','Yellow','Orange','OrangeRed','Red', \
            #    'Crimson','Fuchsia','Indigo','DarkCyan','White'] 
            else:
                radarcbar = ['Lavender', 'Thistle', 'Plum', 'MediumPurple', 'CornFlowerBlue', 'SkyBlue', 'PaleTurquoise', 'LightCyan', 'Yellow', 'Gold', 'Orange', 'DarkOrange', 'Chocolate', 'IndianRed', 'FireBrick', 'Maroon']
        else: 
            radarcbar = deepcopy(color_list)

        temp_cmap = colors.ListedColormap(radarcbar)
        self.temp_cmap = temp_cmap # export for use with trajectory object


#############################################################################################################


    def set_hid_colorbar(self, color_list=None):

        if color_list is None:
         hidcbar = deepcopy(self.hid_colors)

        else:
            hidcbar = deepcopy(color_list)
        self.hid_cmap = colors.ListedColormap(hidcbar)

        #self.boundshid = np.arange(0,12)
        self.boundshid = np.arange(self.hid_cmap.N+1)
        self.normhid = colors.BoundaryNorm(self.boundshid, self.hid_cmap.N)
#############################################################################################################
    def set_cs_colorbar(self, color_list=None):
        if color_list is None:
            cscbar = deepcopy(self.cs_colors)

        else:
            cscbar = deepcopy(color_list)

        self.cs_cmap = colors.ListedColormap(cscbar)
        self.cs_bounds = np.arange(0,5)
        self.cs_norm = colors.BoundaryNorm(self.cs_bounds, self.cs_cmap.N)

#############################################################################################################

    def set_lims(self, var_key, varlims):
        "Reset limits of a variable for plotting, feed it a 2 element list"
        self.lims[var_key] = varlims
#############################################################################################################

    def get_lims(self, var_key):
        return self.lims[var_key]
#############################################################################################################

    def set_delta(self, var_key, new_delta):
        self.delta[var_key] = new_delta
#############################################################################################################

    def get_delta(self, var_key):
        return self.delta[var_key]
#####################################
