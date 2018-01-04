# This just houses some general functions

import numpy as np
from copy import deepcopy
import scipy.ndimage as ND
import datetime
import matplotlib.dates as md
import sys
import linecache

def percentile95(inputdata, **kwargs):
    return np.percentile(inputdata, 95, **kwargs)

def percentile99(inputdata, **kwargs):
    return np.percentile(inputdata, 99, **kwargs)

def percentile90(inputdata, **kwargs):
    return np.percentile(inputdata, 90, **kwargs)

def percentile1(inputdata, **kwargs):
    return np.percentile(inputdata, 1, **kwargs)

def percentile5(inputdata, **kwargs):
    return np.percentile(inputdata, 5, **kwargs)

def percentile10(inputdata, **kwargs):
    return np.percentile(inputdata, 10, **kwargs)


def print_exception():
    exc_type, exc_obj, tb = sys.exc_info()
    f = tb.tb_frame
    lineno = tb.tb_lineno
    filename = f.f_code.co_filename
    linecache.checkcache(filename)
    line = linecache.getline(filename, lineno, f.f_globals)
    print 'EXCEPTION IN ({}, LINE {} "{}"): {}'.format(filename, lineno, line.strip(), exc_obj)


def size_filter(input_arr, filter, min_volume = 1, filter_type = '!=', \
			bad_val = -1, filter_out = True, filter_bad = -1, verbose = False):
    "Used in a number of other functions, will filter out spurious data based on some filter"

	# will do non-zero feature ID on some input_arr and return that feature output
	# unless filter_arr is not None, in which case it will use that feature output to filter
	# the filter_arr and return that
    small_volumes = 0. # keep track of number of volumes filtered out
    big_volumes = 0. # keep track of volumes that are big enough
    label_structure = np.ones((3,3,3)) # allows all 3d neighbors for ID
    out = deepcopy(input_arr) # copy of array to do filtering on and return
    dummy_arr = deepcopy(input_arr) # make a dummy array the same as input to do stuff with
#        bad = dummy_arr != filter # if not equal to filter, set it to 0 to not be ID'ed by ND label function
	# this one gonna try to use strings and evals
    bad = eval("dummy_arr %s filter"%filter_type) # if not equal to filter, set it to 0 to not be ID'ed by ND label function
    dummy_arr[bad] = 0 # now have only one certain HID identifier
    labelarray, nfeat = ND.label(dummy_arr, structure=label_structure)
	#print '%d features'%nfeat
        # now have the HIDs labelled
    for feat in range(1, nfeat+1): # dont do 0
    	points = np.where(labelarray == feat)
	npoints = points[0].shape[0]
        if npoints < min_volume:
#                print 'FILTERING OUT A VOLUME'
	    small_volumes += 1
            labelarray[points] = bad_val # resetting the feature array to 0 if too small
	    if filter_out: out.mask[points] = True
	else: big_volumes += 1
    if verbose:
	print '%d volumes large enough, %d filtered out'%(big_volumes, small_volumes)
    if filter_out:
    	if verbose: print 'Returning filtered array'
	return out
    else: # actually use the labelarray as a filter to modify some other array
	if verbose: print 'Returning feature array'
	return labelarray # return the labelarray, mainly for debugging purposes

#######################################################################################################################

def nder (y, n, smooth=0, delta = 1, normalize = 1): # nth derivative using gradient recursively
    out = np.zeros(y.shape) # out is same shape as input array
    if smooth > 0:
        for itime in range(y.shape[-1]): # loop thru time array

            out[:,itime] = ND.gaussian_filter1d(y[:,itime],smooth) # first smooth it for each time
            for i in range(n): # do n derivatives
                out[:,itime] = np.array(np.gradient(out[:,itime], delta)) # convert gradient output to array
            if normalize == 1:
                if out[:,itime].any(): # if any are not 0, how to normalize?
                    out[:,itime] = out[:,itime]/np.abs(out[:,itime]).max()
    return out


#######################################################################################################################


def pretty_time(ax):
    xfmt = md.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(xfmt)

#######################################################################################################################


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.itervalues():
        sp.set_visible(False)
#######################################################################################################################


def shadow_line_plot(x, y, ax, fgcolor='white', bgcolor='black', fgwidth=1.5, bgwidth=2.5, alpha=1.0, **kwargs):
    # kwargs can be passed to the foreground line plot
    ax.plot(x, y, bgcolor, linewidth=bgwidth, alpha=alpha)
    line_handle = ax.plot(x, y, fgcolor, linewidth=fgwidth, alpha=alpha, **kwargs)
    return line_handle
#######################################################################################################################



def get_cells_times(cell_st, track_st, track_num, year=2011):

    c = np.where(track_st.num == track_num)
    track_cell = track_st[c]
    cells = np.array([np.int(tc) for tc in track_cell.cells[0].split(',')[:-1]])
    
    times = []
    for i in range(cells.shape[0]):
        wh_cell = np.where(cell_st.num == cells[i])[0][0]
        this_cell = cell_st[wh_cell]
        jday = np.int(np.floor(this_cell.time))
        hr, mn, sec = int(this_cell.time_str[0:2]), int(this_cell.time_str[3:5]), int(this_cell.time_str[6:8])
        ftime = datetime.datetime(year,1,1) + datetime.timedelta(days = jday-1) + datetime.timedelta(hours = hr) + \
                datetime.timedelta(minutes = mn) + datetime.timedelta(seconds = sec)
        times.append(ftime)

    return cells, np.array(times)

#######################################################################################################################

def hist2d(varx=None, vary=None, binsx=None, binsy=None):
# This will just look at the whole volume

    vars_q = (varx is not None) and (vary is not None)
    bins_q = (binsx is not None) and (binsy is not None)

    if vars_q:

# grab all values and unravel them (if that's what it's called), but only between the valid indices
        if (hasattr(varx, 'mask') and hasattr(vary, 'mask')):
            outmask = np.logical_and(varx.mask, vary.mask)
            varx.mask = outmask
            vary.mask = outmask
            varx_ravel = np.ma.compressed(varx)
            vary_ravel = np.ma.compressed(vary)

        else:
            mskx = (np.isfinite(varx))
            msky = (np.isfinite(vary))
            msk = np.logical_and(mskx,msky)
            varx_ravel = np.ravel(varx[msk])
            vary_ravel = np.ravel(vary[msk])


        if not bins_q:
            binsx = 20
            binsy = 20
#         print vary_ravel
#         print varx_ravel
#         print binsy
#         print binsx

        hist, edges = np.histogramdd((vary_ravel, varx_ravel), bins=(binsy, binsx), normed=True)
        hist *= 100.0/np.sum(hist)
        return hist, edges



    else:
        print 'Need to provide proper variables to the varx and vary keywords'
        return 

#######################################################################################################################
def make_cmap(colors, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The rgb values
    may either be in 8bit [0 to 255] (in which bit must be set to True when
    called) or arithmetic [0 to 1] (default). make_cmap returns a cmap with
    equally spaced colors.
    Arange you tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    '''
    import matplotlib
    bit_rgb = np.linspace(0,1,256)
    position = np.linspace(0,1,len(colors))
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    return matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

if __name__ == '__main__':
    fig = plt.figure()
    ax = fig.add_subplot(211)
### Create a list of RGB tuples
    colors = [(255,0,0), (255,255,0), (255,255,255), (0,157,0), (0,0,255)] # This example uses the 8-bit RGB
### Call the function make_cmap which returns your colormap

    my_cmap = make_cmap(colors, bit=True)
### Use your colormap
    plt.pcolor(np.random.rand(25,50), cmap=my_cmap)
    plt.colorbar()

    ax = fig.add_subplot(212)
    colors = [(1,1,1), (0.5,0,0)] # This example uses the arithmetic RGB
### If you are only going to use your colormap once you can
### take out a step.
    plt.pcolor(np.random.rand(25,50), cmap=make_cmap(colors))
    plt.colorbar()
    plt.show()
