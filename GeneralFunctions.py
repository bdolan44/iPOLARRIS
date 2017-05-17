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

#import analysis_tools as AT
#import lightning_tools as LT
import csu_fhc
import general_tools as gentools
import RadarConfig
from matplotlib.colors import from_levels_and_colors
import GeneralFunctions as GF

# Up here are just some general functions


########################
def _get_ab_incides(hts, above=None, below=None):

    if above is not None:
        bot_index = (np.argsort(np.abs(hts) - above))[0]
    else:
        bot_index = np.int(np.argmin(hts))
    if below is not None:
        top_index = (np.argsort(np.abs(hts) - below))[0]
        #print 't',top_index
    else:
        top_index = np.int(np.argmax(hts))

        print bot_index
        print 't',top_index

    return int(bot_index), int(top_index)

########################
def cfad(data = None,cfad =None,hts=None,value_bins=None, above=2.0, below=15.0,tspan=None, pick=None, norm=True,ret_z=0,ret_bin = 0,z_resolution=1.0,z_ind = 0,
        thresh=-900.0):
# pick a variable and do a CFAD for the cell

    if value_bins is None : # set a default if nothing is there
        value_bins = np.linspace(np.nanmin(data), np.nanmax(data), 20)
    else:
        pass

    if hts is None:
        print 'Please provide Nominal heights'
        return

    nbins = value_bins.shape[0]
    bot_index, top_index = _get_ab_incides(hts,above=above, below=below)

    delz = hts[1]-hts[0]
    if np.mod(z_resolution, delz) != 0:
            print 'Need even multiple of vertical resolution: {d.1f}'.format(d = delz)
            return

    multiple = np.int(z_resolution/delz)
    looped = np.arange(0, int(hts.shape[0]), multiple)
    cfad_out = np.zeros((int(np.shape(looped)[0]), nbins-1))

    for ivl, vl in enumerate(looped):

        dum1 = data.reshape(data.shape[z_ind],-1)
        dum2 = dum1[vl:vl+multiple,...]
        dum2[dum2<thresh]=np.nan
        dum3 = np.where(np.isfinite(dum2))
        
        lev_hist, edges = np.histogram(dum2[dum3], bins=value_bins, density=True) 
        if norm == True:
            lev_hist = 100.0*lev_hist/np.sum(lev_hist)
        if np.max(lev_hist) > 0:
            cfad_out[ivl, :] = lev_hist
    if ret_z == 1 and ret_bin == 1:    
        return cfad_out,hts[looped],value_bins
    elif ret_z == 1 and ret_bin == 0:
        return cfad_out,hts[looped]
    elif ret_z == 0 and ret_bin == 1:
        return cfad_out, value_bins
    else:
        return cfad_out

#############################################################################################################

def cfad_plot(var,data = None,cfad=None, hts=None, nbins=20, ax=None, maxval=10.0, above=2.0, below=15.0, bins=None, 
        log=False, pick=None, z_resolution=1.0,levels=None,tspan =None,cont = False, rconf = None,**kwargs):

    if hts is None:
        print 'please provide nominal heights to cfad_plot'
        return
    if data is not None:
        cfad, reshts, vbins =GF.cfad(data, hts, value_bins=bins, above=above, below=below, pick=pick, z_resolution=z_resolution,tspan=tspan,
                z_ind = 0,ret_z = 1,ret_bin=1)
    elif cfad is not None:
        if bins == None:
            vbins = np.arange(np.nanmin(data),np.nanmax(data),nbins)
        else:
            vbins = bins[:-1]
        reshts = hts
    else:
        print 'please specify data or cfad'
        return
 
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

    if cont is True:
        levs = [0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.]
        cols = ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet']
        pc = ax.contourf(bins[0:-1],reshts,cfad_ma,levs,colors=cols)
    else:

        if levels is not None:
            cmap, norm = from_levels_and_colors([0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,15.0,20.,25.], ['silver','darkgray','slategrey','dimgray','blue','mediumaquamarine','yellow','orange','red','fuchsia','violet']) # mention levels and colors here
            #print cmap
            pc = ax.pcolormesh(bins, reshts, cfad_ma, norm=norm, cmap=cmap)
        else:
            
            pc = ax.pcolormesh(bins, reshts, cfad_ma, vmin=0, vmax=maxval, norm=norm, **kwargs)

    cb = fig.colorbar(pc, ax=ax)
    cb.set_label('Frequency (%)')
    ax.set_ylabel('Height (km MSL)')
#        try:
    if rconf is not None:
        if var in rconf.names.keys():
            
            ax.set_xlabel('{n} {u}'.format(n=rconf.names[var], u=rconf.units[var]))
            #print rconf.print_title(tm=tspan)
            ax.set_title("{d}".format(d=rconf.print_title(tm=tspan)))
    #            ax.set_title('%s %s %s CFAD' % (self.print_date(), self.radar_name, self.longnames[var]))
#        except:
#            pass

    return fig, ax


#############################################################################################################

#############################################################################################################

def hist2d(datax, datay, varx=None, vary=None, binsx=None, binsy=None, above=None, below=None, pick=None,xthr = -900.0,ythr=-900.0):
#        if above is None:
    # This will just look at the whole volume
        
#       bot_index, top_index = self._get_ab_incides(above=above, below=below)
    if binsx is None:
        binsx = np.arange(np.nanmin(datax),np.nanmax(datax),20)
    if binsy is None:
        binsy = np.arange(np.nanmin(datay),np.nanmax(datay),20)
    
    
    
    dumx = np.less(datax,xthr)
    dumy = np.less(datay,ythr)
    datax[dumx]=np.nan
    datay[dumy]=np.nan

    hist, edges = gentools.hist2d(varx=datax, vary=datay, 
            binsx=binsx, binsy=binsy)

#         hist, edges = gentools.hist2d(varx=datax[bot_index:top_index], vary=datay[bot_index:top_index], 
#                 binsx=binsx, binsy=binsy)

    return hist, edges


#############################################################################################################
#############################################################################################################

def plot_2dhist(hist,edge,ax=None,cbon = True,rconf = None):
    if ax is None:
        fig, ax = plt.subplots()
    else:
    # ax has been passed in, do nothing to ax, but need to get the parent fig
        fig = ax.get_figure()

    cb = ax.contourf(edge[1][:-1],edge[0][:-1],hist,norm=colors.Normalize(vmin=0, vmax=np.max(hist)),levels=np.arange(0.01,np.max(hist),0.01))
    if cbon == True:
        print ' making colorbar'
        plt.colorbar(cb,ax=ax)
    if rconf is not None:
            
            #print rconf.print_title(tm=tspan)
        ax.set_title("{d}".format(d=rconf.print_title()))
    #            ax.set_title('%s %s %s CFAD' % (self.print_date(), self.radar_name, self.longnames[var]))
#        except:
#            pass


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

    def vertical_hid_volume(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None):
        # This function only returns the height and volume arrays, no plotting
        if above is None:
            above = 0
        if below is None:
            below = self.nhgts

        #bot_index, top_index = self._get_ab_incides(above=above, below=below)
        #data = self._pick_data(self.data[self.hid_name].data, pick)
        data = self.data[self.hid_name].sel(z=slice(above,below)).data
        msk = (np.less(self.data[self.dz_name].data, -900))
        data[msk] = -1
        #print 'vhv dat',np.shape(data)
        # Not sure if I need this here.....
        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return

        multiple = np.int(z_resolution/self.dz)
        vol = np.zeros(int(self.data[self.z_name].data.shape[1]/multiple))
        hts = np.zeros(int(self.data[self.z_name].data.shape[1]/multiple))
        #print np.shape(vol)
        #print self.data[self.z_name].data.shape[1]
        looped = np.arange(0, int(self.data[self.z_name].data[0].shape[0]), multiple)
        #print looped,multiple
        for vi,vl in enumerate(looped):
            lev_hid = data[:,vl:vl+multiple,...] # go to vl+multiple cuz not inclusive
            #print 'lev_hid',np.shape(lev_hid)
#            print hid_nums, np.shape(lev_hid)
            where_this_hid = np.where(lev_hid == hid_nums)
#            print np.shape(where_this_hid)
            this_hid_vol = where_this_hid[0].shape[0]
            vol[vi] += this_hid_vol
            #print self.data[self.z_name].data[0][vl+multiple]
        hts = self.data[self.z_name].data[0][looped]

        #print hts
#         print self.data[self.z_name].data[0]
#         print self.data[self.z_name].data[0][::looped]
#         hts = self.data[self.z_name].data[::multiple]
#        print 'hts in vert hid vol', np.shape(hts)
        return hts, vol


#############################################################################################################
#############################################################################################################

    def hid_vertical_fraction(self, hid_nums, z_resolution=1.0, above=None, below=None, pick=None):


        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return None
        else:
            #hts = self.data[self.z_name].data[0][::int(np.round(z_resolution/self.dz))]
            
            hts2 = np.squeeze(self.data[self.z_name].sel(t=slice(0,1),z=slice(above,below)).data)
            #hts2 = (self.data[self.z_name].sel(t=slice(0,-1),z=slice(above,below))).data
            #print np.shape(hts2)
            hts = hts2[::int(np.round(z_resolution/self.dz))]
            # This gives the percent of the storm that is taken up by the combo of HID values
        hid_nums = np.asarray(hid_nums)

        hidcdf = self.hid_cdf(z_resolution=z_resolution)
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

    def hid_cdf(self, z_resolution=1.0, pick=None):
        # vertical HID_cdf with bar plots I think
        if np.mod(z_resolution, self.dz) != 0:
            print 'Need even multiple of vertical resolution: %.1f'%self.dz
            return None

        # loop thru the species and just call the vertical hid volume
        all_vols = []
        for sp in range(len(self.species)):
            all_vols.append(self.vertical_hid_volume([sp+1], z_resolution=z_resolution, pick=pick)[1]) # need the +1


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


    def plot_hid_cdf(self, data=None, z_resolution=1.0, ax=None, pick=None):
        # this will just plot it

        if data is not None:
            pass
        else:
            pass # will call the hid_cdf function here
            data = self.hid_cdf(z_resolution=z_resolution, pick=pick)
        #print np.shape(data)
        if ax is None:
            fig, ax = plt.subplots(1,1)
        else:
            fig = ax.get_figure()   

        fig.subplots_adjust(left = 0.07, top = 0.93, right = 0.87, bottom = 0.1)
        multiple = np.int(z_resolution/self.dz)

        for i, vl in enumerate(np.arange(0, self.data[self.z_name].data.shape[1], multiple)):
            #print vl,i
#            print self.data[self.z_name].data[vl]
            #print data[0,:]
            ax.barh(self.data[self.z_name].data[0][vl], data[0, i], left = 0., edgecolor = 'none', color = self.hid_colors[1]) 
            for spec in range(1, len(self.species)): # now looping thru the species to make bar plot
                ax.barh(self.data[self.z_name].data[0][vl], data[spec, i], left = data[spec-1, i], \
                color = self.hid_colors[spec+1], edgecolor = 'none')
        ax.set_xlim(0,100)
        ax.set_xlabel('Cumulative frequency (%)')
        ax.set_ylabel('Height (km MSL)')
        # now have to do a custom colorbar?
#        self.HID_barplot_colorbar(fig)  # call separate HID colorbar function for bar plots

            #fig.suptitle('%04d/%02d/%02d - %02d:%02d:%02d %s, cell %d, HID CDF' \
            #                %(self.year,self.month,self.date,self.hour,self.minute,self.second, \
            #                self.radar, self.cell_num), fontsize = 14)
#        ax.set_title('%s %s Hydrometeor identification CDF' % (self.print_date(), self.radar_name))

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
            
        uw = np.zeros(self.data[self.z_name].data[0,...].shape)
        #temps = self.data[self.z_name].data[0,:]
        # basically just loop thru the Z and get the associated temperature and area

        data = self.data[self.w_name].data
        if thresh_dz == True:
            data[self.data[self.dz_name].data < -900.0]=np.nan
        for iz, z in enumerate(self.data[self.z_name].data[0]):
            values_above = np.where(data[:,iz,...] >= thresh)[0]
            num_above = len(values_above)
            uw[iz] = num_above*self.dx*self.dy/self.ntimes
            if self.data[self.x_name].units == "[deg]":
                #print 'changing units'
                uw[iz]=uw[iz]*110.*110.

        #print np.shape(uw)
        #print np.shape(self.T[0,:,0,0])
        # now inerpolate this to the temps listed
        f_temp_u = sint.interp1d(self.T[0,:,0,0], uw, bounds_error=False)
        uwp_interp = f_temp_u(temps)
        #return temps, uw
        return temps,uwp_interp
