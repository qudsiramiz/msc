from spacepy.pycdf import CDF as cdf
from numpy import *
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf

import time

start = time.time()

mu_0 = 4 * pi * 1.e-7
m_p = 1.67e-27

# Name of the spacecraft

msc = 'cassini'

# Location of unprocessed data

data_dir =  '/scratch/maruca_cdaweb/cassini/mag_1min/cdf/'

# Location where the processed data is to be saved

save_dir = '/home/ahmadr/SpacePlasFest/msc/cassini/data/processed/'

# Search and sort all the file in the data_dir

fnames = sort( glob( data_dir + '*.cdf' ) )

print( 'A total of %s relevant files found for spacecraft %s' %( len( fnames ), msc ) )

df = None

# Define a list where all the DataFrames will be stored

ldf = []

count =0

'''
B_comp: CDF_DOUBLE [244643, 3]
B_comp_labl: CDF_CHAR*2 [3] NRV
B_mag: CDF_DOUBLE [244643]
Epoch: CDF_EPOCH [244643]
FGM_mode: CDF_BYTE [244643]
VHM_mode: CDF_BYTE [244643]
num_points: CDF_BYTE [244643]
sc_lat: CDF_DOUBLE [244643]
sc_long: CDF_DOUBLE [244643]
sc_pos: CDF_DOUBLE [244643, 3]
sc_pos_labl: CDF_CHAR*5 [3] NRV
sc_r: CDF_DOUBLE [244643]
'''

for f in fnames :

	# Read all the files and save relevant parameters to a dictionary

	d = {}
	dat = cdf( f )

	d['Epoch'] = dat['Epoch']

	d['br'] = dat['B_comp'][:,0][:]
	d['bt'] = dat['B_comp'][:,1][:]
	d['bn'] = dat['B_comp'][:,2][:]

	d['bm'] = dat['B_mag'][:]

	d['fgm_mode'] = dat['FGM_mode'][:]
	d['vhm_mode'] = dat['VHM_mode'][:]

	d['heliographicLatitude']  = dat['sc_lat'][:]
	d['heliographicLongitude'] = dat['sc_long'][:]

	d['x_hgi'] = dat['sc_pos'][:,0][:]
	d['y_hgi'] = dat['sc_pos'][:,1][:]
	d['z_hgi'] = dat['sc_pos'][:,2][:]

	d['sc_r'] = dat['sc_r'][:]

	# Change the dictionary to a dataframe
	# This is done because further processing of data like resampling and stuff
	# is easier and faster in a dataframe compared to a simple dictionary and/or
	# arrays

	df = pd.DataFrame( d, index=d['Epoch'] )

	# Replace all the bad datapoints with NaN

	df[(df < -9.99e30) | (df > 9.99e30)] = NaN

	# Resample the entire dataset to 60s period

	dfr = df.resample( '3600S' ).median( )

	# Assign the UTC timezone to the resample dataframe

	dfr.index=dfr.index.tz_localize('UTC')

	dfr.insert(0,'ssepoch',pd.Series(dfr.index,index=dfr.index).apply(datetime.datetime.timestamp))

	# Assign the name of the file in which data is to be saved

	fn = save_dir + fnames[0][-48:-8] + '.hf'

#	dfr.to_pickle( fn, protocol=2 )

	hdf=hf.File(fn,'w')
	hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
	pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
	for i in dfr.columns[1:]: hdf.create_dataset(i,data=np.array(dfr[i]))
	hdf.close()

	# Append the dataframe to the list

	ldf.append( dfr )

	count += 1
	print( 'Data saved for %s' %( f[-48:] ) )
	print( count )

# Save all dataframes to one dataframe and save it in a new file

df_t = pd.concat( ldf )

fn = fnames[0][-48:-8] + '_' + fnames[-1][-16:-8] + '.hf'

#dfr.to_pickle( fn, protocol=2 )

hdf=hf.File(fn,'w')
hdf.create_dataset('datetime',data=(pd.Series(df_t.index)-\
pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
for i in df_t.columns[1:]: hdf.create_dataset(i,data=np.array(df_t[i]))
hdf.close()
