import numpy as np
import h5py as hf
import pandas as pd
import time as tm
import pytz
import datetime

from glob import glob
from spacepy.pycdf import CDF as cdf

'''
dat = hf.File( '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager1_coho1hr_merged_mag_plasma_19770101_20181201_v04.hf' )

tv1 = 1103068800 # December 15, 2004 ==> Termination shock crossing date for Voyager 1

# For Voyager 1 (before Termination Shock)

ind = np.where( dat['datetime'][:] <= tv1 )

datn = {}

for key in dat.keys():

	datn[key] = dat[key][:][ind[0]]

fn = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager1_coho1hr_merged_mag_plasma_19770101_20041215_v04.hf'

hdf=hf.File(fn,'w')

for i in datn.keys():
	hdf.create_dataset(i,data=np.array(datn[i]))
hdf.close()

# For Voyager 1 (after Termination Shock)

ind = np.where( dat['datetime'][:] >= tv1 )

datn = {}

for key in dat.keys():

	datn[key] = dat[key][:][ind[0]]

fn = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager1_coho1hr_merged_mag_plasma_20041215_20181201_v04.hf'

hdf=hf.File(fn,'w')

for i in datn.keys():
	hdf.create_dataset(i,data=np.array(datn[i]))
hdf.close()

dat = []
datn = []
ind = []

# For Voyager 2 (before Termination Shock)

dat = hf.File( '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager2_coho1hr_merged_mag_plasma_19770101_20181201_v04.hf' )

tv2 = 1188432000 # August 30, 2007 ==> Termination shock crossing date for Voyager 2

ind = np.where( dat['datetime'][:] <= tv2 )

datn = {}

for key in dat.keys():

	datn[key] = dat[key][:][ind[0]]

fn = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager2_coho1hr_merged_mag_plasma_19770101_20070830_v04.hf'

hdf=hf.File(fn,'w')

for i in datn.keys():
	hdf.create_dataset(i,data=np.array(datn[i]))
hdf.close()

# For Voyager 2 (after Termination Shock)

ind = np.where( dat['datetime'][:] <= tv2 )

datn = {}

for key in dat.keys():

	datn[key] = dat[key][:][ind[0]]

fn = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v04/voyager2_coho1hr_merged_mag_plasma_20070830_20181201_v04.hf'

hdf=hf.File(fn,'w')

for i in datn.keys():
	hdf.create_dataset(i,data=np.array(datn[i]))
hdf.close()
'''
