# Brody Fuchs, May 2018
# brfuchs@atmos.colostate.edu

# Converting the Yuter and Houze (1998) and the Steiner and Houze (1995)
# Convective/stratiform partitioning algorithm to python
# Just gonna leave the IDL documentation here

import numpy as np
import scipy.ndimage as ndi
from copy import deepcopy




# ; Name:
# ;     CONV_STRAT_LATLON
# ; Purpose:
# ;     This program takes a 2D Cartesian-gridded reflectivity field 
# ;     and designates each point convective or stratiform based upon
# ;     either the Steiner, Houze and Yuter (SYH 1995 - JAMS) or  
# ;     Yuter and Houze (YH 1998 - QJMS) algorithm.
# ; Calling sequence:
# ;     CONV_STRAT_LATLON
# ; Input:
# ;     zfield      2-D radar reflectivity data (dBZ)
# ;     x         2-D X-axis coordinates
# ;     y         2-D Y-axis coordinates
# ;TJL
# ;NOTE: xy are assumed to be lat/lon
# ;     CoreThresh  Intensity threshold to identify convective cores, 
# ;                 represents a minimum (dBZ > thresh = core). YH uses 46.

# ; Optional Keywords:
# ;     method      Indicates whether to use the SYH or YH methodology for
# ;                 for convective core calculations. 
# ;     SYH         SYH (1995) used where radius around center following profile
# ;                     used in the paper
# ;     YH          YH (1998) used where  can be tuned for climate/region to 
# ;                     confirm stratiform region by evidence of uniformity in  
# ;                     horiz (e.g. laminar flow).  a and b to be found by 
# ;                     comparison with subjectively analyzed convective 
# ;                     core points and tuning the cosine function to fit. 
# ;     a           Tunig variable from YH [8. in YH]
# ;     b           Tunig variable from YH [64. in  YH]
# ;     missing     Value of missing data
# ; Output:
# ;     cs      2-D array (x,y) with the conv-strat definition produced by
# ;             this algorithm bkgnd is the mean background reflectivity field
# ;             =-99999. missing data
# ;             =0       stratiform
# ;             =1       convective (z < 25 dBZ)
# ;             =2       convective (25 < z < 30 dBZ)
# ;             =3       convective (30 < z < 35 dBZ)
# ;             =4       convective (35 < z < 40 dBZ)
# ;             =5       convective (z > 40 dBZ)
# ; Keywords:
# ;     None
# ; Author and history:
# ;     20 Apr 2011  Created by Nick Guy (CSU)
# ;                  Modified to speed up performance. 
# ;     27 Apr 2011  Corrected the mean reflectivity of background field
# ;                  to take mean of linear Z and not mean of dBZ.
# ;     22 Dec 2011  Arrays initialized to Nan.
# ;                  Missing value added and updated in procedure.
# ;===================================================================

re = 6371.0

def assign_radius(In1, outarray, latarray, lonarray, dist_thresh, box_size, assign_value):
#need distance threshold, box size and it will return the close_d_global

    for ii in range(0, len(In1[0])):
        #print 'latlon',np.shape(latarray)
        lat_ind = In1[0][ii]
        lon_ind = In1[1][ii]
        lat1 = latarray[lat_ind, lon_ind]
        lon1 = lonarray[lat_ind, lon_ind]

    
        i_lat_min = np.max([0, lat_ind-box_size])
        i_lat_max = np.min([latarray.shape[1], lat_ind+box_size])

        i_lon_min = np.max([0, lon_ind-box_size])
        i_lon_max = np.min([lonarray.shape[0], lon_ind+box_size])
    

    
        #print 'lat range: {} - {}, lon range: {} - {}'.format(i_lat_min, i_lat_max, i_lon_min, i_lon_max)

        lat2 = latarray[i_lat_min:i_lat_max+1, i_lon_min:i_lon_max+1]		
        lon2 = lonarray[i_lat_min:i_lat_max+1, i_lon_min:i_lon_max+1]		
        
    #	    print 'lats: {}'.format(lat2)
    #	print 'lat shape: {}'.format(lat2.shape)

        d = haversine(lat1, lon1, lat2, lon2)
    #	    print d

        close_d = np.where(d <= dist_thresh)
            #print close_d



        if len(close_d[0]):
    #             if assign_value == 5:
    #                print 'box_size',close_d[1]+lon_ind-box_size
    #             if all(close_d[1]+lon_ind-box_size) > 0:
#             if ((lon_ind-box_size)> 0 and (lon_ind-box_size)<i_lon_max):
            close_d_global = (close_d[0]+lat_ind-box_size, close_d[1]+lon_ind-box_size)
            outarray[close_d_global] = assign_value
#             else:
#                     close_d_global = (close_d[0]+lat_ind-box_size, close_d[1]+lon_ind-box_size)
#                     outarray[close_d_global] = assign_value
            

    return outarray


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    return c * re




def conv_strat_latlon(dbz, lat2d, lon2d, CoreThresh=42.0, method='SYH', a=8, b=64, tune_thresh=42.0, sm_rad=15, fill_dbz=25.0, bg_diff=10):

#    xdim = N_ELEMENTS(x[*,0])-1
#    ydim = N_ELEMENTS(y[0,*])-1

    assert method in ['SYH', 'YH'] # make sure it's either SYH or YH

#    delx = 0.01
#    dely = 0.01
    lon1d = lon2d[0, :]
    lat1d = lat2d[:, 0]

    #print x
    #print y

    #print lon
    #print lat

    dlat = np.average(np.diff(lat1d))
    dlon = np.average(np.diff(lon1d))
    #print 'dlat: %.2e, dlon: %.2e'%(dlat, dlon)

#    print lat[:12]
#    print lon[:12]
	
    # initialize the cs field with all masked values
    cs = np.full(dbz.shape, np.nan)
    cc = np.zeros_like(cs) - 1
    bkgnd = np.full(dbz.shape, np.nan)


    # if dbz is a masked array, then just grab the opposite of the mask
    #print 'steiner',type(dbz)
    if isinstance(dbz, np.ma.masked_array):
        bad = deepcopy(dbz.mask)
    else:
        bad = dbz == np.nan

    good = np.logical_not(bad)
	
#    zlin[bad] = np.nan
#    dbz[bad] = np.nan
    dbz[bad] = fill_dbz
#    zlin[bad] = 0.0

    # calculate linear reflectivity
    zlin = 10.**(dbz/10.)


    # now calculate a smoothed background field using a smoothing function
    # The IDL version of this used FILTER_IMAGE, which is homebrewed

    # I think we'll start with scipy.ndi.median_filter

    # not sure how this is gonna work with nans?
    #print 'shape sm_rad',sm_rad,np.shape(zlin)
    bkgnd_lin = ndi.median_filter(zlin, size=(sm_rad, sm_rad), mode='nearest')
#    bkgnd_lin = ndi.gaussian_filter(zlin, (sm_rad, sm_rad))

    bkgnd = 10.*np.log10(bkgnd_lin)

    bkgnd[bad] = np.nan
#    bkgnd[bad] = -10.0


    inThr = (bkgnd < 0.) & ((dbz - bkgnd) > bg_diff)
    cc[inThr] = 1

    #; Next line uses SYH 1995 radius algorithm
    if method == 'SYH':
        inCon = (bkgnd >= 0.) & (bkgnd < tune_thresh) & (dbz-bkgnd > (bg_diff-(bkgnd**2.)/180.))
#; This line uses YH (1998) climatological tuning algorithm
    else:
        inCon = (bkgnd >= 0.) & (bkgnd < tune_thresh) & (dbz-bkgnd > a*np.cos((np.pi*bkgnd)/(2.*b)))

    cc[inCon] = 2

    inbkCore = bkgnd >= tune_thresh

    cc[inbkCore] = 3

#; Test for convective core threshold
    inCore = dbz >= CoreThresh
    cc[inCore] = 4
#;-------------------------------------------------
#PRINT,'Beginning application of radius'

    # anywhere data is good, give it a 0 at least. 0 is stratiform
    cs[good] = 0

    # anything above 0 is "convective"
    # Check for 


    In1 = np.where((cc > 0) & (bkgnd <= 25.))
    if len(In1[0]): # if there are any points that satisfy this
        cs = assign_radius(In1, cs, lat2d, lon2d, 1.0, 10, 1)

    In2 = np.where((cc > 0) & (bkgnd > 25.) & (bkgnd <= 32.))
    if len(In2[0]):
        cs = assign_radius(In2, cs, lat2d, lon2d, 1.0, 10, 2)

    In3 = np.where((cc > 0) & (bkgnd > 32.) & (bkgnd <= 37.))
    if len(In3[0]):
        cs = assign_radius(In3, cs, lat2d, lon2d, 2.0, 10, 3)

    In4 = np.where((cc > 0) & (bkgnd > 37.) & (bkgnd <= 42.))
    if len(In4[0]):
        cs = assign_radius(In4, cs, lat2d, lon2d, 2.0, 10, 4)

    In5 = np.where((cc > 0) & (bkgnd > 42.))
    if len(In5[0]):
        cs = assign_radius(In5, cs, lat2d, lon2d, 2.0, 10, 5)


    return cs, cc, bkgnd


