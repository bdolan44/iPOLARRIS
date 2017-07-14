# Brody Fuchs, CSU, October 2016
# brfuchs@atmos.colostate.edu

# Classes and function that read in polarimetric/DD files
# and plot them and do things with them

# Also included here is a Cell object to do some cell calculations
# which will be used for my PhD stuff

# 10/13/16: This is hooked into my old version of things like HID, will need to get 
#             csuradartools into this at some point

# 2/22/2016: Trying to modularize to allow for integration with larger datasets such as 
# WRF output and other functions.

from __future__ import division
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
#from pyhid import beta_functions, cdf_fhc, radar_calculations
from matplotlib import colors
import matplotlib.ticker as plticker
from general_tools import make_cmap
import scipy.ndimage as ND
import matplotlib as mpl
from copy import deepcopy
import datetime
import scipy.ndimage as ndi
import xarray as xr
import raintype as rt



#import analysis_tools as AT
#import lightning_tools as LT
import csu_fhc
import general_tools as gentools
import RadarConfig

# Up here are just some general functions


class RadarData(RadarConfig.RadarConfig): 


    def __init__(self, data,times, ddata = None,dz='DZ', zdr='DR', kdp='KD', ldr='LH', rho='RH', hid='HID',
            temp='T', x='x', y='y', z='z', u='U', v='V', w='Wvar', rr='RR',lat=None, lon=None, band='C',exper='CASE',
            radar_name= None,mphys=None,dd_data = None,z_thresh=-10.0,cs_z = 2.0,zconv = 41.,zdr_offset=0, 
            conv_types = ['CONVECTIVE'],strat_types = ['STRATIFORM'],mixed_types = ['UNCERTAIN']): 

        super(RadarData, self).__init__(dz=dz, zdr=zdr, kdp=kdp, ldr=ldr, rho=rho, hid=hid, temp=temp, x=x, y=y,
                      z=z, u=u, v=v, w=w,mphys=mphys,exper=exper,lat=lat,lon=lon,tm = times,radar_name = radar_name)

        # ********** initialize the data *********************
#        self.data = {} 
        #self.date = None
#        self.radar_file = radar_file
#        self.dd_file = dd_file

        # now if radar_file and dd_file have been set, go and try to read them in, else do nothing
#        if self.radar_file is not None:
#            self.time_parse = time_parse
#            self.dformat = dformat
        self.conv_types = conv_types
        self.strat_types = strat_types
        self.mixed_types = mixed_types
        self.band = band
        if ddata is not None:
            self.data= data.combine_first(ddata)
        else:
            self.data = data
        self.zdr_offset = zdr_offset
        if zdr_offset != 0:
            print 'fixing Zdr!'
            self.zdr_correct()
        self.cs_z =cs_z
        self.zconv =zconv
        self.t_name = temp
        self.z_thresh=z_thresh
        self.zind = 1
        self.yind = 2
        self.xind = 3
        self.ntimes =1
#         print self.data.variables['z']
#         print np.shape(self.data['z'].data)
        try:
            self.nhgts = np.shape(self.data[self.z_name].data)[self.zind][0]
        except:
            self.nhgts = np.shape(self.data[self.z_name].data)
#            self.read_data_from_nc(self.radar_file)
        self.calc_deltas()
        self.rr_name = rr
        self.mask_dat()
        self.raintype()
#        else:
#            pass

        # down here goes some zdr checking
        self.zdr_offset = 0 # initialize as 0
        try:
            self.top_index = np.where(self.data[self.z_name]== np.max(self.data[self.z_name]))[self.zind]
        except:
            self.top_index = np.where(self.data[self.z_name]== np.max(self.data[self.z_name]))[0]
        try:
            self.bot_index = np.where(self.data[self.z_name]== np.min(self.data[self.z_name]))[self.zind]
        except:
            self.bot_index = np.where(self.data[self.z_name]== np.min(self.data[self.z_name]))[0]

#         if hasattr(self.nc,'Latitude_deg') == True:
#             print 'Getting attribute!'
#             try:
#                 self.radar_lat = getattr(self.nc,'Latitude_deg')[0]
#             except:
#                 self.radar_lat = getattr(self.nc,'Latitude_deg')
#             
#         if hasattr(self.nc,'Longitude_deg') == True:
#             print 'Getting attribute!'
#             try:
#                 self.radar_lat = getattr(self.nc,'Longitude_deg')[0]
#             except:
#                 self.radar_lat = getattr(self.nc,'Longitude_deg')
        try:
            self.T = self.data[temp]
        except:
            print ''
#        self.set_masks()
#        self.remove_nc()
    # Have some options of how to pass it files? 
    # Nc gridded files
    # pyart?
    # some other format??

    def read_data_from_nc(self, filename, format = 'grid', date = None,dd_flag =0):
    # will need to change time_parse if filename is different, also assumes 4 digits for year
    # this will read in a netcdf file calling functions below assuming a certain filename format
    # RADyyyymmdd_hhmmss*.ext
        tp = self.time_parse[0]
        te = self.time_parse[1]
        self.filename = filename
        self.base = os.path.basename(filename)
        self.radar_name = self.base[:3]
        if date is None:
            try:
                radcdate=np.str(self.base[tp:tp+te])
                #print radcdate
                self.date=datetime.datetime.strptime(radcdate,self.dformat)
            except ValueError:
                pass

        # run some basic functions so not needed to be typed out
        self.nc4_read()
        self.nc4_varnames()
        
#        self.nc4_dict()

#############################################################################################################

    def nc4_read(self):
        "Read netcdf4 file with netCDF4"
        self.nc = Dataset(self.filename)
    #############################################################################################################

    def nc4_varnames(self):
        "Returns all of the variables of the netcdf4 file"
        self.ncvariables = self.nc.variables.keys()

    #############################################################################################################

    def raintype(self):
        print 'Setting up default C/S parameters. Can change in RadarData.py'
        minZdiff = 20 
        deepcoszero = 40
        shallowconvmin = 28
        #truncZconvthres = 43;
        #dBZformaxconvradius = 46;
        truncZconvthres = self.zconv
        dBZformaxconvradius = self.zconv+3;
        weakechothres = 7
        backgrndradius = 5       #(in km)
        maxConvRadius = 10       #(in km)
        minsize = 8              #(in km^2)
        startslope = 50          #(in km^2)
        maxsize = 2000           #(in km^2)

        if self.x_name == 'longitude':
            dx = self.dx*110.
        else:
            dx = self.dx
#        print dx
        zlev = self.get_ind(self.cs_z,self.data[self.z_name].data)
#        print zlev
        
        refl = np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev)).data)
        if len(np.shape(refl)) > 2:
            refl = np.squeeze(self.data[self.dz_name].sel(z=slice(zlev,zlev+1)).data)

        raintype,self.rtypes = rt.raintype(refl, refl_missing_val=self.data[self.dz_name].data.min(), 
                                   refl_dx=dx, minZdiff=minZdiff, deepcoszero=deepcoszero,
                                   shallowconvmin=shallowconvmin,truncZconvthres=truncZconvthres,
                                   dBZformaxconvradius=dBZformaxconvradius,
                                   weakechothres=weakechothres, backgrndradius=backgrndradius,
                                   maxConvRadius=maxConvRadius,minsize=minsize,
                                   startslope=startslope, maxsize=maxsize)
        nlevs = np.shape(self.data[self.z_name].data)[0]
#        print nlevs
        rpt = np.tile(raintype,(nlevs,1,1))
        self.raintype= rpt        
        self.def_convstrat()
    #############################################################################################################

    def check_size(self,ddata):
        if self.x_name in self.data:
            ochksz = np.shape(self.data[self.x_name])
#            print self.nc.variables[self.x_name]
            ddsz = np.shape(ddata.variables[self.x_name])
            
            if ochksz != ddsz:
                "Create a new array of the dimensions for the pol file to put DD data in"
                xvals = ddata.variables[self.x_name][:]
                yvals = ddata.variables[self.y_name][:]
                zvals = ddata.variables[self.z_name][:]
                
                whx = np.nonzero(np.in1d(self.data[self.x_name], xvals))[0]
                why = np.nonzero(np.in1d(self.data[self.y_name], yvals))[0]
                whz = np.nonzero(np.in1d(self.data[self.z_name], zvals))[0]

                
                for key in ddata.variables.keys():
                    if ddata.variables[key].ndim == 3:
                        test =np.zeros_like(self.data[self.dz_name])
                        test = np.ma.masked_where(test==0,test)
                        test[whz[0]:whz[-1]+1,why[0]:why[-1]+1,whx[0]:whx[-1]+1] = ddata.variables[key][:]
    #                    for i,v in enumerate(whz):
    #                        for j,u in enumerate(why):
    #                            for k,w in enumerate(whx):
    #                                test[v,u,w]=self.nc.variables[key][i,j,k]
#                        print np.shape(test)
                        if not key in self.data.variables.keys():
#                            print key
                            self.data[key]=((self.data.dims['z'],self.data.dims['y'],self.data.dims['x']),test)
                
                
    #############################################################################################################

    def nc4_dict(self):
        "Goes into the netcdf file and grabs all of the data and puts into a dictionary"
        #out = {}
        for key in self.nc.variables.keys():
            self.data[key] = self.nc.variables[key][:]

    #self.data[data_type] = out
    #############################################################################################################

    def get_data_from_dict(self, dict):    
        # if want to pass a dictionary already
        self.data = dict
#############################################################################################################

    def mask_dat(self):    
        # if want to pass a dictionary already
        for k in self.data.variables.keys():
#            print k, np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data)
            if k != 'TIME' and k != self.t_name:
#                print k
#                print np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data)
                try:
                    whbad = np.where(np.less_equal(self.data[self.dz_name].data, self.z_thresh))
                    self.data[k].data[whbad]= np.nan
                    whbad2 = np.where(np.isnan(self.data[self.dz_name].data))
                    self.data[k].data[whbad2]= np.nan
#                    self.data[k].data = np.ma.masked_where(self.data[self.dz_name].data < self.z_thresh,self.data[k].data)
                except:
                    pass
#                    print np.shape(self.data[self.dz_name].data), np.shape(self.data[k].data), 'wrong shapes!'
#############################################################################################################

    def convert_t(self):    
        # if want to pass a dictionary already
        self.T=np.ma.masked_less(self.T,0)
        self.T = self.T-273.15
        
#############################################################################################################
    def set_masks(self):
        for k in self.data.keys():
            self.data[k] = np.ma.masked_where(self.data[k] < -998.0,self.data[k])

    def valid_vars(self):
        return np.intersect1d(self.pol_vars, self.data.keys())
#############################################################################################################


    def calc_deltas(self): # get grid sizes for x, y, z
        try:
            self.dx = np.average(np.abs(np.diff(self.data[self.x_name].data[0,:])))
        except:
            self.dx = np.average(np.abs(np.diff(self.data[self.x_name].data)))
        try:
            self.dy = np.average(np.abs(np.diff(self.data[self.y_name].data[:,0])))
        except:
            self.dy = np.average(np.abs(np.diff(self.data[self.y_name].data)))
        try:
            self.dz = np.average(np.abs(np.diff(self.data[self.z_name].data[:])))
        except:
            self.dz = np.average(np.abs(np.diff(self.data[self.z_name].data)))
            


#############################################################################################################


    def add_field(self, in_array, name, pol_var=False):
        "Add some field to the radar object and give it a name"
        self.data[name] = in_array
        if pol_var:
            self.add_pol_var(name)
#############################################################################################################

    def add_pol_var(self, name):
        self.pol_vars = np.append(self.pol_vars, name)
#############################################################################################################

    def get_names(self):
        pass

#############################################################################################################

    ##### some private worker bee functions #####

    def _get_ab_incides(self, above=None, below=None):

        if above is not None:
            bot_index = np.argsort(np.abs(self.data[self.z_name] - above))[0]
        else:
            bot_index = deepcopy(self.bot_index)
        if below is not None:
            top_index = np.argsort(np.abs(self.data[self.z_name] - below))[0]
        else:
            top_index = deepcopy(self.top_index)

        return int(bot_index), int(top_index)


    def _pick_data(self, data, pick):

        if pick is None:
            return data
        else:
            # want to mask the values that are false, hence the logical not
            return np.ma.masked_where(np.logical_not(pick), data) 


    def convert_xy_to_ll(self, x_data, y_data):
    # code to use the radar lat/lon to convert x and y coords to lat lon
        lat = self.radar_lat + y_data/111.0
        lon = self.radar_lon + x_data/(111.0*np.cos(self.radar_lat*np.pi/180.0))

        pass


    def convert_ll_to_xy(self, lat_data, lon_data):
        # code to use the radar lat/lon to convert lat/lon data to x/y radar relative coords
        # just going to assume a flat plane for now
        if (self.radar_lat is not None) and (self.radar_lon is not None):
            dlat = lat_data - self.radar_lat
            dlon = lon_data - self.radar_lon    

            x = dlon*111.0*np.cos(self.radar_lat*np.pi/180.0)
            y = dlat*111.0

            return x, y

        else:
            print 'Need to set the radar lat/lon'
            return None, None
    


#############################################################################################################

    ################### Some quality control stuff ###############################################

    def zdr_check(self, bins=np.arange(-3, 3, 0.15), thresh=0.3):
            # This will do a binned histogram of the ZDR everywhere and makes sure that it peaks around 0
        hist, edges = np.histogram(np.isfinite(self.data[self.zdr_name].data), bins=bins)
        dbin = np.diff(bins)[0]/2.0
        maxarg = np.argmax(hist)

        offset = edges[maxarg]+dbin
#        print hist, edges
        self.zdr_offset = deepcopy(offset)
        if np.abs(offset) >= thresh:
            self.zdr_correct()

    def zdr_correct(self):
        self.data[self.zdr_name].data -= self.zdr_offset




    
    def find_melting_level(self):
        # this sets the levels that are above and below the melting level
        mlevel = self.get_T_height(0)

        # now need to 
        self.above_melting = np.zeros((self.data[self.dz_name].shape), dtype=bool)
        self.below_melting = np.zeros((self.data[self.dz_name].shape), dtype=bool)

        # alright, now need to loop thru the vertical
        for iz in range(len(self.data[self.z_name])):
            if self.data[self.z_name][iz] < mlevel:
                self.below_melting[iz, :, :] = True
            elif self.data[self.z_name][iz] > mlevel:
                self.above_melting[iz, :, :] = True
            else:
                pass # this only happens if they're equal I think


#############################################################################################################

    def get_max_loc(self, variable, smooth=0, above=None, below=None):
        # this will find the maximum updraft and return the x, y, z position
        
        #print smooth
        bot_index, top_index = self._get_ab_incides(above=above, below=below)

        if smooth > 0:
            arr = ndi.filters.gaussian_filter(self.data[variable], smooth)[bot_index:top_index]
        else:
            arr = self.data[variable][bot_index:top_index]

        max_ind = np.where(arr == arr.max())
        zloc = self.data[self.z_name][bot_index:top_index][max_ind[0]][0]
        yloc = self.data[self.y_name][max_ind[1]][0]
        xloc = self.data[self.x_name][max_ind[2]][0]

        return zloc, yloc, xloc

    def remove_nc(self):
        del self.nc


    ##################### Here is some sounding stuff ####################################

    def add_sounding_object(self, sobj):
        self.sobj = sobj
        # add the object to be used later or for other purposes
        # now call add_sounding from the sounding object info
        self.add_sounding(height=sobj.data['hght']/1000.0, temp=sobj.data['temp']) # have to flip it to make right side up




    def add_sounding(self, height=None, temp=None, interp=False):
        """This will add sounding data and interp to the radar heights, if you leave data=None, it will put in a basic
        default sounding"""
        if (height is not None) and (temp is not None):
            self.snd_temp = deepcopy(temp)
            self.snd_height = deepcopy(height)
        else: # do a default temp structure
            t0 = 25.0 # base temp
            z0 = 0.5 # base height
            gamma = 6.5 # lapse rate
            nz = 20

            self.snd_height = np.arange(20)+z0
            self.snd_temp = -1*gamma*np.arange(20)+t0
            self.snd_temp[self.snd_temp <= -60] = -60 # make it -60 as the coldest, don't think this will matter too much

#############################################################################################################


    def interp_sounding(self):

        """This will take the sounding data and interpolate it to the radar coordinates"""
        self.gridded_height = np.zeros(self.data[self.dz_name].shape)
        for i in range(self.data[self.z_name].shape[0]):
                 self.gridded_height[i,:,:] = self.data[self.z_name][i]

        self.T = np.interp(self.gridded_height, self.snd_height, self.snd_temp)

#############################################################################################################


    def get_T_height(self, temp, interp=False):
        temp_index = np.argmin(np.abs(self.T[:,0,0] - temp))
        return self.gridded_height[:,0,0][temp_index]


#############################################################################################################
   ############ Here is calling CSU radartools for HID, RR, etc... ############################
#############################################################################################################
    def calc_pol_analysis(self):
        self.set_hid(use_temp = 'True',band=self.band,zthresh = self.z_thresh)
        self.calc_qr_pol()
        self.calc_rr_pol()


#############################################################################################################
   ############ Here is calling CSU radartools for HID, RR, etc... ############################
#############################################################################################################


   # Just a wrapper on the CSU radartools HID function
    def set_hid(self, band=None, use_temp=False, name='HID',zthresh = -9999.0):
       bad = self.data[self.dz_name].data < self.z_thresh
#       print type(self.data[self.dz_name])

       
       if band is None:
           self.hid_band = self.band
       else:
           self.hid_band = band
       if use_temp and hasattr(self, 'T'):
            #print 'using T!'
            #print self.dz_name
            #print self.zdr_name
            #print self.kdp_name
            #print self.rho_name
            #print np.ma.max(self.T)
#             try:
#                 #msk = self.data[self.dz_name].data.where(self.data[self.dz_name].data < -10)
# 
#                 self.data[self.dz_name].data = np.ma.masked_where(self.T.mask,self.data[self.dz_name].data)
#                 self.data[self.zdr_name].data = np.ma.masked_where(self.T.mask,self.data[self.zdr_name].data)
#                 self.data[self.kdp_name].data = np.ma.masked_where(self.T.mask,self.data[self.kdp_name].data)
#                 self.data[self.rho_name].data = np.ma.masked_where(self.T.mask,self.data[self.rho_name].data)
#             except:
#                 self.data[self.zdr_name].data = np.ma.masked_where(bad,self.data[self.zdr_name].data)
#                 self.data[self.kdp_name].data = np.ma.masked_where(bad,self.data[self.kdp_name].data)
#                 self.data[self.rho_name].data = np.ma.masked_where(bad,self.data[self.rho_name].data)
# 
            #print self.hid_band
#             print self.data[self.dz_name]
#             print self.data[self.zdr_name]
#             print self.data[self.kdp_name]
#             print self.data[self.rho_name]
#           
            self.scores = csu_fhc.csu_fhc_summer(dz=self.data[self.dz_name].data, zdr=self.data[self.zdr_name].data, rho=self.data[self.rho_name].data, 
                                kdp=self.data[self.kdp_name].data, band=self.hid_band, use_temp=True, T=self.T)
       else:
           self.scores = csu_fhc.csu_fhc_summer(dz=self.data[self.dz_name].data, zdr=self.data[self.zdr_name].data, rho=self.data[self.rho_name].data, 
                                kdp=self.data[self.kdp_name].data, band=self.hid_band, use_temp=False) 


       self.data[self.dz_name].data[bad] = np.nan
       dzmask = np.isnan(self.data[self.dz_name].data)

           # set the hid
       self.hid = np.argmax(self.scores, axis=0)+1
       try:
#           print 'Trying to mask HID!'
           self.hid[dzmask] = -1
       except:
           print 'Problem trying to threshold HID'
       try:
            scmask = np.isnan(self.scores[0,...])
#            print 'trying to mask via scores!'
            self.hid[scmask] =-1
       except:
            print 'Cant threshold on scores!'
#           self.hid=np.ma.masked_where(self.data[self.dz_name].data < self.z_thresh,self.hid)
       self.add_field((self.data[self.dz_name].dims,self.hid,), name)

#############################################################################################################
        ###Calculate the mixing ratios and add to the radar object

    def calc_qr_pol(self):
    ### This is where I do mixing ratio calculations ### 
        #print self.band
        #print self.expr
#        print self.band, self.exper
        if self.band == 'C' or self.band == 'c':
 #           print 'got C!'
#            print (self.exper)
            if self.exper == 'MC3E':
                k_c= 1.4195
                k_m=0.7489

                z_c=0.0014297
                z_m=0.6729
            if self.exper == 'TWPICE': 
                k_c= 0.89191
                k_m=0.6529

                z_c=0.00194
                z_m=0.5957
        elif self.band == 'S' or self.band == 's':
            if self.exper == 'MC3E':
                k_c= 1.4489
                k_m=0.6589

                z_c=0.0019393
                z_m=0.592032
            if self.exper == 'TWPICE': 
                k_c= 2.559
                k_m=0.76687

                z_c=0.0014377
                z_m=0.66743
        else:
            print 'Your wavelength has not been run yet! Please return to fundamentals.'
            return
            


        dbzz = np.ma.masked_less(self.data[self.dz_name].data,35)
        kdpp = np.ma.masked_where(dbzz.mask,self.data[self.kdp_name].data)

        kdppp = np.ma.masked_less(kdpp,0.2)
        dbzzz = np.ma.masked_where(kdppp.mask,dbzz)

        M_Kdp = k_c*(kdppp**k_m)

        kdp_notmet = np.ma.masked_greater_equal(self.data[self.kdp_name].data,0.2)
        dbz_notmet = np.ma.masked_where(kdp_notmet.mask,self.data[self.dz_name].data)

        dbz_notmett = np.ma.masked_greater_equal(dbz_notmet,35.0)
        kdp_notmett = np.ma.masked_where(dbz_notmett.mask,kdp_notmet)

        kdp_notmettt = np.ma.masked_less(kdp_notmett,-100.)
        dbz_notmettt = np.ma.masked_where(kdp_notmettt.mask,dbz_notmett)

        linz = 10.**(dbz_notmettt/10.)  # for method w/o HID
        M_Z = z_c*(linz**z_m)
    
        M_Z.set_fill_value(0.0)
        M_Kdp.set_fill_value(0.0)
        lwccc = M_Z.filled()+M_Kdp.filled()
        qrr = lwccc/1.225
        qrr = np.ma.masked_less_equal(qrr,0.0)

        self.add_field((self.data[self.dz_name].dims,qrr,), 'rqr')

    
    def calc_rr_pol(self):

#        import pydisdrometer as pyd
#        import pytmatrix as pyt
        import csu_blended_rain_julie

    ### This is where I do mixing ratio calculations ### 

        if self.band == 'C' or self.band == 'c':
            if self.expr == 'MC3E':
                k_c=20.59
                k_m=0.75

                z_c=0.0179
                z_m=0.6855
            
                azzdr_coeff=0.00842
                bzzdr_coeff=0.894
                czzdr_coeff=-3.4029
            
                a=azdrk_coeff=33.0097
                b=bzdrk_coeff=0.9148
                c=czdrk_coeff=-0.92778
            if self.expr == 'TWPICE': 
                k_c=28.0669
                k_m=0.8049

                z_c=0.0152889
                z_m=0.738
            
                azzdr_coeff=0.008097
                bzzdr_coeff=0.87966
                czzdr_coeff=-2.6350
            
                a=azdrk_coeff=41.0249
                b=bzdrk_coeff=0.882359
                c=czdrk_coeff=-1.30948
        elif self.band == 'S' or self. band == 's':
            if self.expr == 'MC3E':
                k_c=35.77975
                k_m=0.767

                z_c=0.017850
                z_m=0.6817
            
                azzdr_coeff=0.00463
                bzzdr_coeff=0.9539
                czzdr_coeff=-3.9243
            
                a=azdrk_coeff=96.8418
                b=bzdrk_coeff=0.965216
                c=czdrk_coeff=-1.9881
            if self.expr == 'TWPICE': 
                k_c=51.7402
                k_m=0.81788

                z_c=0.015364
                z_m=0.73256
            
                azzdr_coeff=0.006482
                bzzdr_coeff=0.92738
                czzdr_coeff=-4.0552
            
                a=azdrk_coeff=94.1082
                b=bzdrk_coeff=0.9337
                c=czdrk_coeff=-1.9350
        else:
            print 'Sorry, your wavelength has not been run yet! Return to first principles!'
            return

        rr,rm = csu_blended_rain_julie.csu_hidro_rain(self.data[self.dz_name].data,self.data[self.zdr_name].data,self.data[self.kdp_name].data,z_c,z_m,k_c,k_m,azdrk_coeff,bzdrk_coeff,
                                                      czdrk_coeff,azzdr_coeff,bzzdr_coeff,czzdr_coeff,band='S',fhc=self.hid)

        self.add_field((self.data[self.dz_name].dims,rr,), 'RRB')
        self.add_field((self.data[self.dz_name].dims,rm,), 'RRM')

        return
#############################################################################################################

    def score_xsec_plot(self, y=None, title_flag=False, *args, **kwargs):

        "Cross section showing scores for all species, used as HID analysis tool"
        # first, get the appropriate y index from the y that's wanted
        if y is None:
            y_ind = int(len(self.data[self.y_name])/2.0)
        else:
            y_ind = np.argmin(np.abs(y - self.data[self.y_name]))

        fig, ax = plt.subplots(5,2, figsize = (9,9))
        
        axf = ax.flatten()
        fig.subplots_adjust(top = 0.90, bottom = 0.12, left = 0.05, right = 0.92, hspace = 0.4)

        for c in range(10): # looping thru each number corresponding to each species

            dummy = axf[c].pcolormesh(self.data[self.x_name], self.data[self.z_name], self.scores[c,:,y_ind,:], cmap = plt.cm.gist_heat_r,
                     vmin = 0, vmax = 1.0)


        axf[c].set_title(self.species[c])

        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])
        cb = fig.colorbar(dummy, cax = cbar_ax, orientation = 'horizontal')
        cb.set_label('$\mu$ score')

        if title_flag:
            fig.suptitle('%s %s Cross Section HID scores' %(self.print_date(), self.radar_name), fontsize = 14)


        return fig, ax

#############################################################################################################

    def score_cappi_plot(self, z=1.0, title_flag=False, *args, **kwargs):

        "CAPPI plot showing scores for all species, used as HID analysis tool"

        # first, get the appropriate z index from the z that's wanted in altitude
        z_ind = np.argmin(np.abs(z - self.data[self.z_name]))


        fig, ax = plt.subplots(5,2, figsize = (9,9))
        
        axf = ax.flatten()
        fig.subplots_adjust(top = 0.90, bottom = 0.12, left = 0.05, right = 0.92, hspace = 0.4)

        for c in range(10): # looping thru each number corresponding to each species

            dummy = axf[c].pcolormesh(self.data[self.x_name], self.data[self.y_name], self.scores[c,z_ind,:,:], cmap = plt.cm.gist_heat_r,
                     vmin = 0, vmax = 1.0)

        axf[c].set_title(self.species[c])

        cbar_ax = fig.add_axes([0.15, 0.06, 0.7, 0.03])
        cb = fig.colorbar(dummy, cax = cbar_ax, orientation = 'horizontal')
        cb.set_label('$\mu$ score')

        if title_flag:
            fig.suptitle('%s %s CAPPI HID scores' %(self.print_date(), self.radar_name), fontsize = 14)


        return fig, ax




### NOW WE GET DOWN TO THE PLOTTING FUNCTIONS HERE


########### STARTING WITH CROSS SECTIONS ################

    def xsec(self, var, y=None, xlim=None, zlim=None, ts = None,varlist=None, ax=None, title_flag=False, 
                vectors=False, cblabel=None, res=2.0,cbpad=0.03, **kwargs):
        "Just one axis cross-section plot of a variable"
        # first, get the appropriate y index from the y that's wanted
        ts=self.date
        tsi = 0
        if y is None:
            y_ind = int(len(self.data[self.y_name].data)/2.0)
        else:
            if self.y_name == 'latitude':
                y_ind = self.get_ind(y,self.data[self.y_name].data)
            else:
                y_ind = y

        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
                xmin = xlim[0]
                xmax = xlim[1]
            else:
                xmini, xmaxi = xlim
                xmin , xmax = xlim
                
        if zlim is None:
            zmin, zmax = self.data[self.z_name].data.min(), self.data[self.z_name].data.max()
            zmini = self.get_ind(zmin,self.data[self.z_name].data)
            zmaxi = self.get_ind(zmax,self.data[self.z_name].data)
        else:
            zmin, zmax = zlim
                    # If ax is not given, open a fig and ax object. This is not advisable
        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()

        # now actually doing the plotting

        #print self.data[self.x_name].shape, self.data[self.z_name].shape, self.data[var][:,y_ind,:].shape
        #print self.data[self.z_name]
    
        #print xmini,xmaxi,zmini,zmaxi
        # if this variable is already included in the defaults, then this is straightforward
       #print tsi, tsi, zmini,zmaxi,xmini,xmaxi,y_ind,var
#        print zmini,zmaxi,y_ind,xmini,xmaxi
        if self.y_name == 'latitude':
            data = np.squeeze(self.data[var].sel(z=slice(zmini,zmaxi),y=slice(y_ind,y_ind+1),x=slice(xmini,xmaxi)).data)
        else:
            data = np.squeeze(self.data[var].sel(z=slice(zmini,zmaxi),y=slice(y_ind,y_ind),x=slice(xmini,xmaxi)).data)

#         if np.shape(data) > 2:
#             data = np.squeeze(self.data[var].sel(z=slice(zmini,zmaxi),x=slice(xmini,xmaxi)).data)
            

        try:
            xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi),y=slice(y_ind,y_ind+1)))
        except:
            xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi)))
#        print zmini,zmaxi
        
        zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi)))
        
        data = np.ma.masked_less(data,-900.0)
        data = np.ma.masked_where(~np.isfinite(data),data)
#        print np.shape(data)
        #print np.shape(data),np.shape(xdat),np.shape(zdat)
        #print np.shape(xdat),np.shape(zdat)
#        print 'data',np.shape(data),'zdat',np.shape(zdat),'xdat',np.shape(xdat)
        if var in self.lims.keys():
            range_lim = self.lims[var][1] - self.lims[var][0]
            dummy = ax.pcolormesh(xdat,zdat, data,
                vmin = self.lims[var][0], vmax = self.lims[var][1], cmap = self.cmaps[var], **kwargs)
        else:
            dat = self.data[var].data
            data[data<-900.0]=np.nan
            range_lim  = np.nanmax(dat) - np.nanmin(dat)
            dummy = ax.pcolormesh(xdat,zdat, data,
                vmin = np.nanmin(dat), vmax = np.nanmax(dat), **kwargs)
        if range_lim < 1:
            cb_format = '%.2f'
        if range_lim >= 1:
            cb_format = '%.1f'
        if range_lim >= 10:
            cb_format = '%d'



        cb = fig.colorbar(dummy, ax=ax, fraction=0.03, format=cb_format, pad=cbpad)
        if var in self.lims.keys():
            cb.set_label(' '.join([self.names[var], self.units[var]]).strip())
            cb.set_ticks(np.arange(self.lims[var][0], self.lims[var][1]+self.delta[var], self.delta[var]))
            cb.set_ticklabels(self.ticklabels[var])
        else:
            cb.set_label(var)


        
        ###### this sets the limits #######
#        print zmin, zmax
        if self.x_name == 'longitude':
            ax.axis([xmin, xmax, zmin, zmax])
        else:
            ax.axis([xmin, xmax, zmin, zmax])
        ax.set_xlabel('Distance E of radar (km)')
        ax.set_ylabel('Altitude (km MSL)')


        if vectors:
            try:
                #print zlim
                self.xsec_vector(ax=ax, y=y,zlim=zlim,xlim=xlim,ts=ts,res=res)
            except Exception, e:
                print 'Error trying to plot xsec vectors: {}'.format(e)

        if title_flag:
            ax.set_title('%s %s Cross Section' %(ts, self.radar_name), fontsize = 14)

        return dummy

#############################################################################################################

    def xsec_multiplot(self, y=None, xlim=None, zlim=None, ts=None,varlist=None, vectors=False,res=2.0, **kwargs):
        "multipanel cross-section plot showing all available polarimetric variables and HID, if available"

    # first, get the appropriate y index from the y that's wanted

        ts=self.date
        tsi = 0
        if y is None:
            y_ind = int(len(self.data[self.y_name].data)/2.0)
            
        else:
            y_ind = self.get_ind(y,self.data[self.y_name].data)
        try:
            yval = self.data[self.y_name].data[y_ind][0]
        except:
            yval = self.data[self.y_name].data[y_ind]
        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
            else:
                xmini, xmaxi = xlim
        if zlim is None:
            zmini, zmaxi = self.data[self.z_name].data.min(), self.data[self.z_name].data.max()
        else:
            zmini = self.get_ind(zlim[0],self.data[self.z_name].data)
            zmaxi = self.get_ind(zlim[1],self.data[self.z_name].data)

        # first get how many varialbes there are?
        if varlist is not None:
            good_vars = varlist
        else:
            good_vars = self.valid_vars()
        nvars = len(good_vars)
        if hasattr(self, 'scores'):
            nvars += 1
        if nvars <= 3:
            ncols = 1
            nrows = deepcopy(nvars)
            figx = 7
            figy = 3*nrows
        else:
            ncols = 2
            nrows = int(np.ceil(nvars/2))
            figx = 10
            figy = 3*nrows
    
        fig, ax = plt.subplots(nrows, ncols, figsize = (figx, figy), sharex = True, sharey = True)
        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): ax = np.array([ax])
        axf = ax.flatten()


        # BF 3/30/16: TAKING OUT IMSHOW AND PUTTING IN PCOLORMESH
        for i, var in enumerate(sorted(good_vars)):
            dummy = self.xsec(var, ts=ts, y=y, vectors=vectors, xlim=xlim, zlim=zlim, ax=axf[i],res=res, **kwargs)
        # now do the HID plot, call previously defined functions

            fig.tight_layout()
            fig.subplots_adjust(top = 0.9)

        fig.suptitle('%s %s Cross Section y = %s' %(ts, self.radar_name,yval), fontsize = 18)

        return fig, ax

    def get_ind(self,val,dat):
        dum = np.abs(val - dat)
        wh_t = np.where(dum == np.min(dum))
        return wh_t[0][0]
        
######################### Here is the CAPPI stuff ##############################

    def cappi(self, var, z=1.0, xlim=None, ylim=None, ax=None,ts = None, title_flag=False, vectors=False, cblabel=None, labels=True, res = 2.0, thresh_dz=False,**kwargs):
        "Just make a Constant Altitude Plan Position Indicator plot of a given variable"

        # first, get the appropriate z index from the z that's wanted in altitude
        #z_ind = np.argmin(np.abs(z - self.data[self.z_name].data))
        z_ind = self.get_ind(z,self.data[self.z_name].data)


        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
                xmin = xlim[0]
                xmax = xlim[1]
            else:
                xmini, xmaxi = xlim
        if ylim is None:
            ymini, ymaxi = self.data[self.y_name].data.min(), self.data[self.y_name].data.max()
    
        else:
            if self.y_name == 'latitude':
                ymini = self.get_ind(ylim[0],self.data[self.y_name].data[:,0])
                ymaxi = self.get_ind(ylim[1],self.data[self.y_name].data[:,0])
                ymin = ylim[0]
                ymax = ylim[1]

            else:
                ymini, ymaxi = ylim
            
        ts=self.date
        tsi = 0

        # If ax is not given, open a fig and ax object. This is not advisable
        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()

        #print np.shape(self.data[var])
        try:
            data = np.squeeze(self.data[var].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi),y=slice(ymini,ymaxi)).data)
        except:
            data = np.squeeze(self.data[var].data[z_ind,ymini:ymaxi,xmini:xmaxi])
        if len(np.shape(data)) > 2:
#            print 'data shape is wrong!',np.shape(data)
            data = data[0,...]
        try:
            xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi),y=slice(ymini,ymaxi)).data)
        except:
            xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi)))
        try:
            ydat = np.squeeze(self.data[self.y_name].sel(x=slice(xmini,xmaxi),y=slice(ymini,ymaxi)).data)
        except:
            ydat = np.squeeze(self.data[self.y_name].sel(y=slice(ymini,ymaxi)))
        
        
#        print 'xmini, xmaxi, xmin,xmax',xmini,xmaxi,xmin,xmax,ymini,ymaxi
#        print xdat[xmax]
#        data[dzmask] =np.nan
        data = np.ma.masked_where(~np.isfinite(data),data)


        if var in self.lims.keys():
            range_lim = self.lims[var][1] - self.lims[var][0]
#            print np.shape(data), np.shape(xdat),np.shape(ydat)
#            print 'in var',var
            dummy = ax.pcolormesh(xdat,ydat, data,
                vmin = self.lims[var][0], vmax = self.lims[var][1], cmap = self.cmaps[var], **kwargs)
        else:
            dat = self.data[var].data
            dat[dat<-900.0]=np.nan
            range_lim  = np.nanmax(dat) - np.nanmin(dat)
            dummy = ax.pcolormesh(xdat,ydat, data,
                vmin = np.nanmin(dat), vmax = np.nanmax(dat), **kwargs)

        if range_lim < 1:
            cb_format = '%.2f'
        if range_lim >= 1:
            cb_format = '%.1f'
        if range_lim >= 10:
            cb_format = '%d'


        if labels:
            cb = fig.colorbar(dummy, ax=ax, fraction=0.03, pad=0.03, format=cb_format)
            if var in self.lims.keys():
                cb.set_label(' '.join([self.names[var], self.units[var]]).strip())
                cb.set_ticks(np.arange(self.lims[var][0], self.lims[var][1]+self.delta[var], self.delta[var]))
                cb.set_ticklabels(self.ticklabels[var])

            else:
                cb.set_label(var)
#                cb.set_ticks(np.arange(self.lims[var][0], self.lims[var][1]+self.delta[var], self.delta[var]))
#                cb.set_ticklabels(self.ticklabels[var])

            
        else: # if this variable is not included in the defaults, have a lot less customization
            # can get around this with the **kwargs
            dummy = ax.pcolormesh(self.data[self.x_name], self.data[self.y_name], self.data[var][z_ind,:,:], **kwargs)

            cb = fig.colorbar(dummy, ax=ax, fraction=0.03, pad=0.03)
            if cblabel is not None:
                cb.set_label(cblabel)



        ####### plotting limits getting set here ######
        if self.x_name == 'longitude':
            ax.axis([xmin, xmax, ymin, ymax])
        else:
            ax.axis([xmini, xmaxi, ymini, ymaxi])
        if labels:
            ax.set_xlabel('Distance E of radar (km)')
            ax.set_ylabel('Distance N of radar (km)')


        # Now check for the vectors flag, if it's there then plot it over the radar stuff
        if vectors:
#            try:
                self.plan_vector(ax=ax, z=z,res=res,thresh_dz=thresh_dz,xlim=xlim,ylim=ylim)
#            except Exception, e:
#                print 'Error trying to plot vectors: {}'.format(e)


        if title_flag:
            ax.set_title('%s %s CAPPI %.1f km MSL' %(ts, self.radar_name, \
                    self.data[self.z_name].data[0][z_ind]), fontsize = 14)

        return dummy

#############################################################################################################

    def cappi_multiplot(self, z=1.0, xlim=None, ylim=None, ts=None,res = 2, varlist=None, vectors=False, **kwargs):
        "6 panel CAPPI plot showing all the polarimetric variables and HID"
        
        # first, get the appropriate z index from the z that's wanted in altitude
        z_ind = self.get_ind(z,self.data[self.z_name].data)
        
        #print xlim, 'In cappi_multiplot'
        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
                xmin = xlim[0]
                xmax = xlim[1]
            else:
                xmini, xmaxi = xlim
#            print 'xlim', np.shape(self.data[self.x_name].data)
        if ylim is None:
            ymini, ymaxi = self.data[self.y_name].data.min(), self.data[self.y_name].data.max()
        else:
            if self.y_name == 'latitude':
                ymini = self.get_ind(ylim[0],self.data[self.y_name].data[:,0])
                ymaxi = self.get_ind(ylim[1],self.data[self.y_name].data[:,0])
                ymin = ylim[0]
                ymax = ylim[1]
            else:
                ymini, ymaxi = ylim
            
        ts=self.date
        tsi = 0



        if varlist is not None:
            good_vars = varlist
        else:
            good_vars = self.valid_vars()
        nvars = len(good_vars)
        if hasattr(self, 'scores'):
            nvars += 1
        if nvars <= 3:
            ncols = 1
            nrows = deepcopy(nvars)
            figx = 5
            figy = 3*nrows
        else:
            ncols = 2
            nrows = int(np.ceil(nvars/2))
        figx = 9
        figy = 3*nrows

        fig, ax = plt.subplots(nrows, ncols, figsize = (figx, figy), sharex = True, sharey = True)
        if not isinstance(ax, np.ndarray) or not isinstance(ax, list): ax = np.array([ax], **kwargs)
        axf = ax.flatten()
        
        for i, var in enumerate((good_vars)):
            #print var    
            dummy = self.cappi(var, z=z, ax=axf[i], xlim=xlim, ylim=ylim,ts = ts, vectors=vectors,res=res)
        # now do the HID plot, call previously defined functions
        # try:
        #     dummy_hid = self.HID_plot(self.HID_from_scores(self.scores, rank = 1)[z_ind,:,:], 
        #             axis = axf[-1],extent=ext)
        #     self.HID_colorbar(dummy_hid, axis = axf[-1], figure = fig, fraction = 0.03, pad = 0.03)
        # except AttributeError:
        #     print 'No HID scores, not plotting'
        #     pass

        fig.tight_layout()
        fig.subplots_adjust(top=0.93)
        fig.suptitle('%s %s CAPPI %.1f km MSL' %(ts, self.radar_name, \
                    self.data[self.z_name][z_ind]), fontsize = 18)


        return fig, ax


#############################################################################################################
    # Down here is dual doppler plotting stuff

    def xsec_vector(self, y=None, xlim=None,zlim=None,ts=None,ax=None, res=2.0, ht_offset=0.2, **kwargs):
        if ts is None:
            ts=self.date
            tsi = 0
#         else:
#             print 'Check your dates!', ts
#        print 'ts:',ts

        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
            else:
                xmini, xmaxi = xlim
    
        if zlim is None:
            zmini = self.get_ind(self.data[self.z_name].data.min(), self.data[self.z_name].data)
            zmaxi = self.get_ind(self.data[self.z_name].data.max(), self.data[self.z_name].data)

        else:
            zmini, zmaxi = zlim

        if np.size(res) == 1:
            resx = res
            resz = res
        else:
            resx = res[0]
            resz = res[1]
        
        if self.data[self.x_name].units == "[deg]":
            skip = int(np.round(resx/self.dx))
            xskip = int(np.round(skip/110.))
            
        else:
            skip = int(np.round(resx/self.dx))
            xskip = skip
        if self.data[self.z_name].units == "[deg]":
            skip = int(np.round(resz/self.dz))
            zskip = skip/110.
        else:
            skip = int(np.round(resz/self.dz))
            zskip = 1
            
        #print skip,xskip, zskip
        #print skip

    # first, get the appropriate y index from the y that's wanted
        if y is None:
            y_ind = int(len(self.data[self.y_name].data[0])/2.0)
        else:
            if self.y_name == 'latitude':
                y_ind = self.get_ind(y,self.data[self.y_name].data[:,0])
            else:
                y_ind = y
#        print 'y in xsec', y

        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()

#        print xmini,xmaxi,y_ind,zmini,zmaxi
        
        try:
            if self.y_name == 'latitude':
#                print y_ind
                xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind+1)).data)
                zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).data)
                udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind+1)).data)
                wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind+1)).data)
            else:
#                print np.shape(xdat), np.shape(zdat)
                xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind)).data)
                zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).data)
                udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind)).data)
                wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind)).data)

        except:
#            print 'uh-oh, exception'
            xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1)).data)
            zdat = np.squeeze(self.data[self.z_name].sel(z=slice(zmini,zmaxi+1)).data)
#            print np.shape(xdat),np.shape(zdat),np.shape(self.data[self.u_name].data)
            udat = np.squeeze(self.data[self.u_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind)).data)
            wdat = np.squeeze(self.data[self.w_name].sel(z=slice(zmini,zmaxi+1),x=slice(xmini,xmaxi+1),y=slice(y_ind,y_ind)).data)
        
     #   print 'vect shp',np.shape(udat),np.shape(wdat),np.shape(xdat),np.shape(zdat)
#        print xskip
#        print np.shape(xdat),np.shape(zdat),np.shape(udat),np.shape(wdat)
#        print ht_offset, xskip
        q_handle = ax.quiver(xdat[::xskip], zdat[::zskip]+ht_offset, \
            udat[::zskip, ::xskip], wdat[::zskip, ::xskip], \
                scale=30, scale_units='inches', pivot='middle', width=0.0025, headwidth=6, **kwargs)
        qk = ax.quiverkey(q_handle, 0.85, 0.85, 20, r'20 $\frac{m}{s}$',coordinates='axes', \
                    fontproperties={'weight':'bold','size':14})        

        return q_handle

#############################################################################################################

    def plan_vector(self, z=1.0, ax=None, xlim=None,ylim = None,ts=None,res=2.0, ht_offset=0.2,thresh_dz=False, **kwargs):


        if np.size(res) == 1:
            resx = res
            resy = res
        else:
            resx = res[0]
            resy = res[1]


        if self.data[self.x_name].units == "[deg]":
            xskip = int(np.round(resx/(self.dx*110.)))
            yskip = int(np.round(resy/(self.dy*110.)))
            #yskip = xskip
            #yskip = int(np.round(res/(self.dy*110.)))
            #skip = int(np.round(skip/110.))
        else:
            xskip = int(np.round(resx/self.dx))
            yskip = int(np.round(resy/self.dy))

        if ts is None:
            ts=self.date
            tsi = 0
        else:
            tsi = self.get_ind(ts,np.array(self.date))

#        print 'xlim',xlim
        if xlim is None:
            xmini, xmaxi = self.data[self.x_name].data.min(), self.data[self.x_name].data.max()
        else:
            if self.x_name == 'longitude':
#                print xlim
                xmini = self.get_ind(xlim[0],self.data[self.x_name].data[0,:])
                xmaxi = self.get_ind(xlim[1],self.data[self.x_name].data[0,:])
            else:
                xmini, xmaxi = xlim
    
        if ylim is None:
            ymini, ymaxi = self.data[self.y_name].data.min(), self.data[self.y_name].data.max()
        else:
            if self.y_name == 'latitude':
                ymini = self.get_ind(ylim[0],self.data[self.y_name].data[:,0])
                ymaxi = self.get_ind(ylim[1],self.data[self.y_name].data[:,0])
            else:
                ymini, ymaxi = ylim

#        print xskip,yskip

        if z is None:
            z_ind = int(len(self.data[self.z_name].data[0])/2.0)
        else:
            z_ind = z

        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()
#        print xmini,xmaxi,ymini,ymaxi

        xdat = np.squeeze(self.data[self.x_name].sel(x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data)
        ydat = np.squeeze(self.data[self.y_name].sel(x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data)
        udat = np.squeeze(self.data[self.u_name].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data)
        vdat = np.squeeze(self.data[self.v_name].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data)
        if thresh_dz == True:
            dzdat = np.squeeze(self.data[self.dz_name].sel(z=slice(z_ind,z_ind+1),x=slice(xmini,xmaxi+1),y=slice(ymini,ymaxi+1)).data)
            #print 'trying to threshold...',np.shape(vdat),np.shape(dzdat)
            msk = np.less(dzdat,self.z_thresh)
            vdat= np.ma.masked_where(msk,vdat)
            udat= np.ma.masked_where(msk,udat)
            xdat= np.ma.masked_where(msk,xdat)
            ydat= np.ma.masked_where(msk,ydat)
            #print type(vdat)
        #print np.max(vdat)
        #print 'vect shp',np.shape(udat),np.shape(vdat),np.shape(xdat),np.min(ydat),np.max(ydat)
#        print np.shape(xdat),np.shape(ydat),np.shape(udat)
        q_handle = ax.quiver(xdat[::yskip][::xskip], ydat[::yskip][::xskip], \
            udat[::yskip][::xskip], vdat[::yskip][::xskip], \
                scale=100, scale_units='inches', pivot='middle', width=0.0025, headwidth=4, **kwargs)

#        q_handle = ax.quiver(self.data[self.x_name][::skip], self.data[self.y_name][::skip]+ht_offset, \
#            self.data[self.u_name][z_ind,:,:][::skip, ::skip], self.data[self.v_name][z_ind,:,:][::skip, ::skip], \
#                scale=100, scale_units='inches', pivot='middle', width=0.002, headwidth=4, **kwargs)
        qk = ax.quiverkey(q_handle, 0.8, 0.8, 20, r'20 $\frac{m}{s}$',coordinates='axes', \
                    fontproperties={'weight':'bold','size':14})        

        return q_handle

#############################################################################################################


    def cfad(self, var, value_bins=None, above=2.0, below=15.0,tspan=None, pick=None, ret_z=0,z_resolution=1.0,cscfad = False):
    # pick a variable and do a CFAD for the cell

        if value_bins is None: # set a default if nothing is there
            value_bins = np.linspace(self.lims[var][0], self.lims[var][1], 20)
        else:
            pass
        #print value_bins

        nbins = value_bins.shape[0]
        bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #print self.data.keys()
#         data = self._pick_data(self.data[var], pick)
#         if 't' in self.data.dims.keys():
#             data = np.mean(data,axis=0)
#             
#         print np.shape(data),type(data)

        if np.mod(z_resolution, self.dz) != 0:
                print 'Need even multiple of vertical resolution: %.1f'%self.dz
                return

        multiple = np.int(z_resolution/self.dz)
        sz=np.shape(self.data[self.z_name].data)[0]
        #print np.shape(sz),sz, multiple
        looped = np.arange(0, sz, multiple)
        cfad_out = np.zeros((sz//multiple, nbins-1))
        #print looped
        #print cfad_out.shape, multiple
#        if tspan == None:
        tsi = 0
        tei = self.ntimes-1
#         else:
#             ts=tspan[0]
#             te=tspan[1]
#             tsi = self.get_ind(ts,np.array(self.date))
#             tei = self.get_ind(te,np.array(self.date))
# 
#        print 'cscfad',cscfad
        if cscfad == 'convective':
            #mask = np.where(self.raintype != 2)
            mask= np.where(self.raintype != 2)
            holddat = deepcopy(self.data[var].data)
            self.data[var].data[mask] = np.nan
#            print 'in conv'
        elif cscfad == 'stratiform':
            mask = np.where(self.raintype != 1)
            holddat = deepcopy(self.data[var].data)
            self.data[var].data[mask] = np.nan
#            print 'in strat',type(self.data[var].data)
        else:
            mask = np.where(self.raintype > 100)

        # if left blank, check the whole thing
        for ivl, vl in enumerate(looped[:-1]):
            #print ivl, vl
#             try:
            #cfad, ed = (np.histogram(mc3e_wrf.data['w'].sel(z=slice(3,15))))            
            #print ts,te,'ts,te'
            v = self.data[self.z_name].data[vl]
            v2 = self.data[self.z_name].data[vl+multiple]
            dum = (self.data[var].sel(z=slice(v,v2)).data)
            #print np.max(dum)
#            dum2 = np.ma.masked_less(dum,-900.0)
            dum2 = np.where(np.isfinite(dum))
            #print dum[dum2]
            lev_hist, edges = np.histogram(dum[dum2], bins=value_bins, density=True) 
#             except:
#            lev_hist, edges = np.histogram(data[vl:vl+multiple].ravel(), bins=value_bins, density=True) 
            #print lev_hist, edges
                    # this keeps it general, can make more elaborate calls in other functions
            lev_hist = 100.0*lev_hist/np.sum(lev_hist)
            if np.max(lev_hist) > 0:
                cfad_out[ivl, :] = lev_hist
        #print np.shape(cfad_out)
        if cscfad == 'convective' or cscfad == 'stratiform':
#            print 'setting data back'
            self.data[var].data = holddat

        if ret_z == 1:
        
            return cfad_out,self.data[self.z_name].data[0][looped]
        else:
            return cfad_out

#############################################################################################################

    def cfad_plot(self, var, nbins=20, ax=None, maxval=10.0, above=2.0, below=15.0, bins=None, 
            log=False, pick=None, z_resolution=1.0,levels=None,tspan =None,cont = False,cscfad = False, **kwargs):

        from matplotlib.colors import from_levels_and_colors
        if bins is None:
            bins = np.linspace(self.lims[var][0], self.lims[var][1], nbins)
        else:
            pass


        multiple = np.int(z_resolution/self.dz)
#         print self.dz
#         print 'multiple: {}'.format(multiple)

        cfad = self.cfad(var, value_bins=bins, above=above, below=below, pick=pick, z_resolution=z_resolution,tspan=tspan,cscfad=cscfad)
    #print cfad.sum(axis=1)
        bot_index, top_index = self._get_ab_incides(above=above, below=below)

        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()

        if log:
            norm = colors.LogNorm(vmin=1e-5, vmax=1e2)
        else:
            norm = None


        # plot the CFAD
        cfad_ma = np.ma.masked_where(cfad==0, cfad)
#        print np.max(cfad_ma),var
        #print np.shape(cfad_ma)
#        print multiple, self.data[self.z_name].data[::multiple]
        if cont is True:
            levs = [0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.]
            cols = ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet']
            pc = ax.contourf(bins[0:-1],self.data[self.z_name].data[::multiple],(cfad_ma/np.sum(cfad_ma))*100.,levs,color=cols)
        else:

            if levels is not None:
                cmap, norm = from_levels_and_colors([0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.], ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet']) # mention levels and colors here
                #print cmap
                pc = ax.pcolormesh(bins, self.data[self.z_name].data[::multiple], cfad_ma, norm=norm, cmap=cmap)
            else:
                pc = ax.pcolormesh(bins, self.data[self.z_name].data[::multiple], cfad_ma, vmin=0, vmax=maxval, norm=norm, **kwargs)

#        print np.shape(cfad_ma)

        cb = fig.colorbar(pc, ax=ax)
        cb.set_label('Frequency (%)')
        ax.set_ylabel('Height (km MSL)')
#        try:
        ax.set_xlabel('%s %s' %(self.names[var], self.units[var]))
        ax.set_title("{d} {r} {v}".format(d=self.date,r=self.radar_name,v=self.longnames[var]))
#            ax.set_title('%s %s %s CFAD' % (self.print_date(), self.radar_name, self.longnames[var]))
#        except:
#            pass

        return fig, ax


#############################################################################################################

#############################################################################################################

    def hist2d(self, varx=None, vary=None, binsx=None, binsy=None, above=None, below=None, pick=None,xthr = -900.0,ythr=-900.0):
        # This will just look at the whole volume
#        if above is None:
            
 #       bot_index, top_index = self._get_ab_incides(above=above, below=below)
        if above is None:
            above = 0
        if below is None:
            below = self.nhgts[0]
        
        datax = self.data[varx].sel(z=slice(above,below)).data
        datay = self.data[vary].sel(z=slice(above,below)).data
                    
        dum = np.less(datax,xthr)
        dumy = np.less(datay,ythr)
        datax[dum]=np.nan
        datay[dumy]=np.nan

        hist, edges = gentools.hist2d(varx=datax, vary=datay, 
                binsx=binsx, binsy=binsy)

#         hist, edges = gentools.hist2d(varx=datax[bot_index:top_index], vary=datay[bot_index:top_index], 
#                 binsx=binsx, binsy=binsy)

        return hist, edges


#############################################################################################################
#############################################################################################################

    def plot_2dhist(self, hist,edge,ax=None,cbon = True):
        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()
    
        cb = ax.contourf(edge[0][:-1],edge[1][:-1],hist.T,norm=colors.Normalize(vmin=0, vmax=np.max(hist)),levels=np.arange(0.01,np.max(hist),0.01))
        if cbon == True:
            #print ' making colorbar'
            plt.colorbar(cb,ax=ax)


        # This will just look at the whole volume
#        if above is None:
        return fig, ax


#############################################################################################################



#               ADD some percentile plots   
#############################################################################################################

    def percentileplots(winds,height,exper,date,ptype,punit,pshort,extra,inst,use_frzht=False,frzhgt=None):
    #    print np.ma.min(np.ma.compressed(winds))
        # Time to make the percentile lists
        percentile_50 = np.zeros([np.size(height)])
        percentile_90 = np.zeros([np.size(height)])
        percentile_99 = np.zeros([np.size(height)])
        mean = np.zeros([np.size(height)])
    #    print np.shape(winds)
    #    print np.shape(height)
        for i,ht in enumerate(height):
            #summed_z = np.ma.sum(winds[:,i,:,:])
            try: 
                percentile_90[i]=np.percentile(np.ma.compressed((winds[i,:,:])),90) #use w/ masked arrays
                #percentile_90[i]=np.percentile(winds[:,i,:,:],90)
                #percentile_90[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_90[i]=np.nan
            try:
                percentile_99[i]=np.percentile(np.ma.compressed((winds[i,:,:])),99) #use w/ masked arrays
                #percentile_99[i]=np.percentile(winds[:,i,:,:],99)
                #percentile_99[i]=np.percentile(winds[i],99)
            except IndexError:
                percentile_99[i]=np.nan
            try:
                percentile_50[i]=np.percentile(np.ma.compressed((winds[i,:,:])),50) #use w/ masked arrays
                #percentile_50[i]=np.percentile(winds[:,i,:,:],50)
                #percentile_50[i]=np.percentile(winds[i],50)
            except IndexError:
                percentile_50[i]=np.nan
            #mean[i]=np.ma.mean(winds[:,i,:,:]) #use w/ masked arrays
            #mean[i]=np.mean(winds[:,i,:,:])

        ## This part is only necessary for the mixing ratio stuff ##
        if use_frzht==True:

            for i,hgt in enumerate(height):
    #            print np.float(hgt),frzhgt+np.float(1)
                if np.float(hgt) > frzhgt+np.float(1):
                    percentile_50[i] = 0
                    percentile_90[i] = 0
                    percentile_99[i] = 0
    #                print percentile_50[i]
    #    print percentile_50 

    #         Now doing the plotting
        p2 = plt.plot(percentile_50,height,color='red',linewidth=2,label='50th')
        p1 = plt.plot(percentile_90,height,color='orange',linewidth=2,label='90th')
        p4 = plt.plot(percentile_99,height,color='black',linewidth=2,label='99th')
    #         p5 = plt.plot(mean,z,color='black',linewidth=2,label='mean')
    #         p5 = plt.plot(mean,z,color='blue',linewidth=2,label='mean')
        lns = p2+p1+p4
    #   lns = p4
        labs = [l.get_label() for l in lns]
        lgd = plt.legend(lns,labs,bbox_to_anchor=(1.6,0.85),prop={'size':14})
    #         plt.title('MC3E 1 May 2011 90th Percentile, 99th Percentile, \nand Mean of Downward Winds',size=16,y=1.08,x=0.78)
        plt.title('{c} {dd} 50th, 90th, and 99th Percentile of {p} {e}'.format(c=exper,dd=date,p=ptype,e=extra),size=16,y=1.08,x=0.78)
        #plt.xlabel('Rainrate (mm hr$^-$$^1$)',size=14)
        plt.xlabel('{p} ({u})'.format(p=ptype,u=punit))
        plt.ylabel('Height (km)',size=14)
        plt.ylim([0,6]) # For qr plots
    #         plt.xlim([-1,1])
        return plt
#############################################################################################################

    def vertical_hid_volume(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None, cscfad = None):
        # This function only returns the height and volume arrays, no plotting
        if above is None:
            above = 0
        if below is None:
            #print self.nhgts[0]
            below = self.nhgts[0]

        #bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #data = self._pick_data(self.data[self.hid_name].data, pick)
        #print above,below,np.shape(self.data[self.hid_name].data)
        if cscfad == 'convective':
            mask= np.where(self.raintype != 2)
            holddat = deepcopy(self.data[self.hid_name].data)
            self.data[self.hid_name].data[mask] = -1
            print 'in conv'
        elif cscfad == 'stratiform':
            mask= np.where(self.raintype != 1)
            holddat = deepcopy(self.data[self.hid_name].data)
            self.data[self.hid_name].data[mask] = -1
            print 'in strat'
        else:
            mask = np.where(self.raintype < 100)
        
        
        data = self.data[self.hid_name].sel(z=slice(above,below)).data
        msk = (np.less(self.data[self.dz_name].data, -900))
        data[msk] = -1
        #print 'vhv dat',np.shape(data)
        # Not sure if I need this here.....
        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return

        multiple = np.int(z_resolution/self.dz)
        vol = np.zeros(int(self.data[self.z_name].data.shape[0]/multiple))
        hts = np.zeros(int(self.data[self.z_name].data.shape[0]/multiple))
        #print np.shape(vol)
        #print self.data[self.z_name].data.shape[1]
        looped = np.arange(0, int(self.data[self.z_name].data.shape[0]), multiple)
        #print looped,multiple
        for vi,vl in enumerate(looped):
            lev_hid = data[vl:vl+multiple,...] # go to vl+multiple cuz not inclusive
            #print 'lev_hid',np.shape(lev_hid)
#            print hid_nums, np.shape(lev_hid)
            where_this_hid = np.where(lev_hid == hid_nums)
#            print np.shape(where_this_hid)
            this_hid_vol = where_this_hid[0].shape[0]
            vol[vi] += this_hid_vol
            #print self.data[self.z_name].data[0][vl+multiple]
        hts = self.data[self.z_name].data[looped]

        if cscfad == 'convective' or cscfad == 'stratiform':
            self.data[self.hid_name].data = holddat

        #print np.shape(vol), np.max(vol)
        #print hts
#         print self.data[self.z_name].data[0]
#         print self.data[self.z_name].data[0][::looped]
#         hts = self.data[self.z_name].data[::multiple]
#        print 'hts in vert hid vol', np.shape(hts)
        return hts, vol


#############################################################################################################
#############################################################################################################

    def hid_vertical_fraction(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None,cscfad = None):


        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return None
        else:
            #hts = self.data[self.z_name].data[0][::int(np.round(z_resolution/self.dz))]
            
            hts2 = np.squeeze(self.data[self.z_name].sel(z=slice(above,below)).data)
            #hts2 = (self.data[self.z_name].sel(t=slice(0,-1),z=slice(above,below))).data
            #print np.shape(hts2)
            hts = hts2[::int(np.round(z_resolution/self.dz))]
            # This gives the percent of the storm that is taken up by the combo of HID values
        hid_nums = np.asarray(hid_nums)

        hidcdf = self.hid_cdf(z_resolution=z_resolution, cscfad = cscfad)
        hvf = np.zeros(hidcdf.shape[1])
#        print 'hvf in hvf', np.shape(hvf)
    # now loop thru each hid_num
        for hn in hid_nums:
            if hn == 1:
                hvf += hidcdf[hn-1, :]
            else:
                hvf += hidcdf[hn-1, :] - hidcdf[hn-2, :]

        return hts, hvf


#############################################################################################################

    def hid_cdf(self, z_resolution=1.0, pick=None,cscfad = None):
        # vertical HID_cdf with bar plots I think
        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return None

        # loop thru the species and just call the vertical hid volume
        all_vols = []
        for sp in range(len(self.species)):
            all_vols.append(self.vertical_hid_volume([sp+1], z_resolution=z_resolution, pick=pick,cscfad = cscfad)[1]) # need the +1


        all_vols = np.array(all_vols)
        all_cdf = np.zeros_like(all_vols)
#9        print np.shape(all_vols)
        # shape is 10,16, which is nspecies x nheights
        # need to do cdf on each level
        for iz in range(all_vols.shape[1]):
            # loop thru the vertical
            level_cum_vol = np.cumsum(all_vols[:, iz])
            all_cdf[:, iz] = 100.0*level_cum_vol/level_cum_vol[-1]

        return all_cdf


#############################################################################################################


    def plot_hid_cdf(self, data=None, z_resolution=1.0, ax=None, pick=None,cscfad = None):
        # this will just plot it

        if data is not None:
            pass
        else:
            pass # will call the hid_cdf function here
            data = self.hid_cdf(z_resolution=z_resolution, pick=pick,cscfad = cscfad)
        #print np.shape(data)
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()   

        fig.subplots_adjust(left = 0.07, top = 0.93, right = 0.87, bottom = 0.1)
        multiple = np.int(z_resolution/self.dz)

        for i, vl in enumerate(np.arange(0, self.data[self.z_name].data.shape[0], multiple)):
            #print vl,i
#            print self.data[self.z_name].data[vl]
            #print data[0,:]
            ax.barh(self.data[self.z_name].data[vl], data[0, i], left = 0., edgecolor = 'none', color = self.hid_colors[1]) 
            for spec in range(1, len(self.species)): # now looping thru the species to make bar plot
                #print spec, np.max(data[spec,i])
                ax.barh(self.data[self.z_name].data[vl], data[spec, i], left = data[spec-1, i], \
                color = self.hid_colors[spec+1], edgecolor = 'none')
        ax.set_xlim(0,100)
        ax.set_xlabel('Cumulative frequency (%)')
        ax.set_ylabel('Height (km MSL)')
        # now have to do a custom colorbar?
#        self.HID_barplot_colorbar(fig)  # call separate HID colorbar function for bar plots

            #fig.suptitle('%04d/%02d/%02d - %02d:%02d:%02d %s, cell %d, HID CDF' \
            #                %(self.year,self.month,self.date,self.hour,self.minute,self.second, \
            #                self.radar, self.cell_num), fontsize = 14)
        ax.set_title('%s %s HID CDF' % (self.print_date(), self.radar_name))

        return fig, ax 


#############################################################################################################
#############################################################################################################

    def vertical_profile(self, var, rep_func=np.average, above=None, below=None, pick=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #print bot_index,top_index
        vp = []
        heights = []

        data = self._pick_data(self.data[var], pick)

        for iz in range(bot_index, top_index):
            if isinstance(data, np.ma.core.MaskedArray):
                this_data = np.ma.compressed(data[iz,:,:])
            else:
                this_data = data[iz,:,:]
            if len(this_data):
                value = rep_func(this_data)
            else:
                value = np.nan
        vp.append(value)
        heights.append(self.data[self.z_name][iz])

        vp = np.ma.masked_invalid(np.array(vp))
        vp = np.ma.masked_where(vp < -100, vp)
        heights = np.array(heights)

        return heights, vp



#############################################################################################################
#############################################################################################################
#############################################################################################################

#############################################################################################################



    def find_updraft_downdraft(self, updraft_thresh=1.0, downdraft_thresh=-1.0):
        # This will make a mask for updraft and downdraft regions based on some thresholds
        # First, will make sure that the updraft thresh is > 0 and the downdraft is < 0
        if (updraft_thresh >= 0) and (downdraft_thresh < 0):
             self.updraft = self.data[self.w_name] >= updraft_thresh
             self.downdraft = (self.data[self.w_name] <= downdraft_thresh) & (self.data[self.w_name] >= -99.0)

        else:
            print 'Make sure your thresholds are valid'
            return



#############################################################################################################

    def updraft_percentile(self, perc):
        return np.percentile(self.data[self.w_name][self.updraft], perc)

#############################################################################################################

    def downdraft_percentile(self, perc):
        #print self.data[self.w_name][self.downdraft]
        try:
            return np.percentile(self.data[self.w_name][self.downdraft], perc)
        except IndexError:
            return None

#############################################################################################################


    def updraft_vol(self, up_val, above=None, below=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)
    
        # if left blank, check the whole thing
        uv = 0. # initialized the updraft volume
        for vl in range(bot_index, top_index): 
            lev_w = self.data[self.w_name][vl]
            in_vol = np.where(lev_w >= up_val)      
            uv += in_vol[0].shape[0]*self.dx*self.dy*self.dz

        return uv

#############################################################################################################

    def downdraft_vol(self, down_val, above=None, below=None):

        bot_index, top_index = self._get_ab_incides(above=above, below=below)

        # if left blank, check the whole thing
        dv = 0. # initialized the updraft volume
        for vl in range(bot_index, top_index): 
            lev_w = self.data[self.w_name][vl]
            in_vol = np.where( (lev_w <= down_val) & (lev_w >= -99.0) )
            dv += in_vol[0].shape[0]*self.dx*self.dy*self.dz

        return dv

#############################################################################################################

    def w_profile(self, rep_func=np.average, masked=True):
        # returns a vertical profile of vertical velocity with a given sign (up, down, all)
        # using a specified representative function. Default is average, but could do something like a median

        # first, set up the out, which will have a spot for each z and one for up/down draft
        vert_motion = deepcopy(self.data[self.w_name])
        out = {'up': np.zeros(self.data[self.z_name].shape), 'down': np.zeros(self.data[self.z_name].shape)}

        for il in range(self.data[self.z_name].shape[0]):
            good_up = np.where(self.updraft[il])[0]
            good_dn = np.where(self.downdraft[il])[0]
            if good_up.shape[0] > 0:
                out['up'][il] = rep_func(vert_motion[il][self.updraft[il]])
            else:
                out['up'][il] = np.nan

            if good_dn.shape[0] > 0:
                out['down'][il] = rep_func(vert_motion[il][self.downdraft[il]])
            else:
                out['down'][il] = np.nan


        if masked:
            out['up'] = np.ma.masked_invalid(out['up'])
            out['down'] = np.ma.masked_invalid(out['down'])


        return out

#############################################################################################################


    def plot_w_profile(self, rep_func=np.average, masked=True, ax=None):
        # this will plot the vertical profile
        wprof = self.w_profile(rep_func=rep_func, masked=masked)

        if ax is None:
            fig, ax = plt.subplots()
        else:
        # ax has been passed in, do nothing to ax, but need to get the parent fig
            fig = ax.get_figure()

        u = ax.plot(wprof['up'], self.data[self.z_name], color='red', linewidth=2, label='updrafts')
        d = ax.plot(wprof['down'], self.data[self.z_name], color='blue', linewidth=2, label='downdrafts')

        max_abs = np.abs(np.array(ax.get_xlim())).max()

        ax.set_xlim(-1*max_abs, max_abs)
        #ax.set_xlim()
        ax.axvline(x=0, color='black')
        ax.set_xlabel('Vertical motion (m/s)')
        ax.set_ylabel('Altitude (km MSL)')
        ax.grid(True)
        ax.legend(loc='best')
        ax.set_title('%s %s Vertical motion profile' % (self.print_date(), self.radar_name))

        return fig, ax


#############################################################################################################

    def updraft_width_profile(self, thresh=5.0, temps=np.arange(20,-60,-5),thresh_dz=False):
        import scipy.interpolate as sint
        # this gets the width of the updraft as a function of temperature
            
        uw = np.zeros(self.data[self.z_name].data[:].shape)
        #temps = self.data[self.z_name].data[0,:]
        # basically just loop thru the Z and get the associated temperature and area

        data = self.data[self.w_name].data
        if thresh_dz == True:
            data[self.data[self.dz_name].data < self.z_thresh]=np.nan
        for iz, z in enumerate(self.data[self.z_name].data):
            values_above = np.where(data[iz,...] >= thresh)[0]
            num_above = len(values_above)
            uw[iz] = num_above*self.dx*self.dy/self.ntimes
            if self.data[self.x_name].units == "[deg]":
                #print 'changing units'
                uw[iz]=uw[iz]*110.*110.
        #print np.shape(uw),np.max(uw),self.data[self.x_name].units
        #print np.shape(uw)
        #print np.shape(self.T[0,:,0,0])
        # now inerpolate this to the temps listed
        f_temp_u = sint.interp1d(self.T[:,0,0], uw, bounds_error=False)
        uwp_interp = f_temp_u(temps)
        #return temps, uw
        return temps,uwp_interp

#############################################################################################################

    def def_convstrat(self):
        
        for k in eval(self.conv_types):
            
            self.raintype[self.raintype == self.rtypes[k]] =2
        for k in eval(self.strat_types):
            self.raintype[self.raintype == self.rtypes[k]] =1

        for k in eval(self.mixed_types):
            self.raintype[self.raintype == self.rtypes[k]] =3

