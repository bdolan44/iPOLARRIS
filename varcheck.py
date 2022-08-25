import netCDF4 as nc
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

in1 = nc.Dataset('samples/wrfout_06_d03_2015-12-03_21:55:00.POLARRIS.matsui2018.nc')
in2 = nc.Dataset('samples/KLGX_20151203_215639_V06.nc')

z1 = in1.variables['hgt'][:]
z2 = in2.variables['z0'][:]

Z1 = in1.variables['zhh01'][:]
Z2 = np.squeeze(in2.variables['REF'][:])

Z1ml = Z1[3,:,:]
Z2ml = Z2[3,:,:]

print(Z1ml)
print(Z2ml)
