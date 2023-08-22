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

msc = 'voyager1'

# Location of unprocessed data

data_dir =  '/scratch/maruca_cdaweb/voyager/' + msc + '/coho1hr_magplasma'

# Location where the processed data is to be saved

save_dir = '/home/ahmadr/SpacePlasFest/msc/voyager/data/processed/v01/'

# Search and sort all the file in the data_dir

fnames = sort( glob( data_dir + '/**/*.cdf' , recursive=True ) )

print( 'A total of %s relevant files found for spacecraft %s' %( len( fnames ), msc ) )

df = None

# Define a list where all the DataFrames will be stored

ldf = []

count = 0

'''
KeysView(<CDF:
ABS_B: CDF_REAL4 [744]
BN: CDF_REAL4 [744]
BR: CDF_REAL4 [744]
BT: CDF_REAL4 [744]
Epoch: CDF_EPOCH [744]
F: CDF_REAL4 [744]
V: CDF_REAL4 [744]
azimuthAngle: CDF_REAL4 [744]
elevAngle: CDF_REAL4 [744]
heliocentricDistance: CDF_REAL4 [744]
heliographicLatitude: CDF_REAL4 [744]
heliographicLongitude: CDF_REAL4 [744]
protonDensity: CDF_REAL4 [744]
protonFlux10_CRS: CDF_REAL4 [744]
protonFlux11_CRS: CDF_REAL4 [744]
protonFlux12_CRS: CDF_REAL4 [744]
protonFlux13_CRS: CDF_REAL4 [744]
protonFlux14_CRS: CDF_REAL4 [744]
protonFlux15_CRS: CDF_REAL4 [744]
protonFlux16_CRS: CDF_REAL4 [744]
protonFlux17_CRS: CDF_REAL4 [744]
protonFlux18_CRS: CDF_REAL4 [744]
protonFlux1_CRS: CDF_REAL4 [744]
protonFlux1_LECP: CDF_REAL4 [744]
protonFlux2_CRS: CDF_REAL4 [744]
protonFlux2_LECP: CDF_REAL4 [744]
protonFlux3_CRS: CDF_REAL4 [744]
protonFlux3_LECP: CDF_REAL4 [744]
protonFlux4_CRS: CDF_REAL4 [744]
protonFlux5_CRS: CDF_REAL4 [744]
protonFlux6_CRS: CDF_REAL4 [744]
protonFlux7_CRS: CDF_REAL4 [744]
protonFlux8_CRS: CDF_REAL4 [744]
protonFlux9_CRS: CDF_REAL4 [744]
protonTemp: CDF_REAL4 [744]
'''

# Figd the index which corresponds to the the Termination shock crossing. And
# separate the fil1e into two parts, before and after the shock.
ind_list = []

for i,f in enumerate(fnames):
	if( f[-16:-8] > '20041201' ) :
		ind_list.append(i)

min_ind = min(ind_list)
max_ind = max(ind_list)

for f in fnames[min_ind:max_ind] :

	# Read all the files and save relevant parameters to a dictionary

	d = {}
	dat = cdf( f )

	d['Epoch'] = dat['Epoch']

	d['br'] = dat['BR']
	d['bt'] = dat['BT']
	d['bn'] = dat['BN']

	d['bm'] = dat['ABS_B']

	d['azimuthAngle'] = dat['azimuthAngle']
	d['elevAngle']   = dat['elevAngle']

	d['sc_r'] = dat['heliocentricDistance']

	d['heliographicLatitude']  = dat['heliographicLatitude']
	d['heliographicLongitude'] = dat['heliographicLongitude']

	d['vp_m'] = dat['V']
	d['np'] = dat['protonDensity']
	d['Tp'] = dat['protonTemp']

	# Change the dictionary to a dataframe
	# This is done because further processing of data like resampling and stuff
	# is easier and faster in a dataframe compared to a simple dictionary and/or
	# arrays

	df = pd.DataFrame( d, index=d['Epoch'] )

	# Replace all the bad datapoints with NaN

	df[(df < -9.99e30) | (df > 9.99e30)] = NaN

	# Compute the components of velocity

	df['vp_r'] = df.vp_m * cos( deg2rad( df.elevAngle ) ) * cos( deg2rad( df.azimuthAngle ) )
	df['vp_t'] = df.vp_m * cos( deg2rad( df.elevAngle ) ) * sin( deg2rad( df.azimuthAngle ) )
	df['vp_n'] = df.vp_m * sin( deg2rad( df.elevAngle ) )

	# Create a partial datafgrame with just a few parameters (magnetic field
	# and density with velocity) which are important for the calculation of
	# Elsasser variables and cross helicity

	dfd = pd.DataFrame( {'br':df.br, 'bt':df.bt, 'bn':df.bn, 
	                     'vp_r':df.vp_r, 'vp_t':df.vp_t, 'vp_n':df.vp_n} )

	# Compute the deviation from mean for all the parameters from partial
	# daraframe. NOTE: Mean is calculated over '24Hours' in this case. This was
	# done so that we always have multiple correlation time lengths and the
	# deviation makes sense from a statistical point

	db = dfd.groupby( pd.Grouper(freq='24H' ) ).apply( lambda x: x-x.mean() )

	# Compute the Elsasser variables and the normalized cross helicity

	df['zpr'] = 1.e3*db.vp_r + 1e-9*db.br/sqrt( mu_0*1.e6*df.np * m_p )
	df['zpt'] = 1.e3*db.vp_t + 1e-9*db.bt/sqrt( mu_0*1.e6*df.np * m_p )
	df['zpn'] = 1.e3*db.vp_n + 1e-9*db.bn/sqrt( mu_0*1.e6*df.np * m_p )
	df['zpm'] = df.zpr**2 + df.zpt**2 + df.zpn**2

	df['zmr'] = 1.e3*db.vp_r - 1e-9*db.br/sqrt( mu_0*1.e6*df.np * m_p )
	df['zmt'] = 1.e3*db.vp_t - 1e-9*db.bt/sqrt( mu_0*1.e6*df.np * m_p )
	df['zmn'] = 1.e3*db.vp_n - 1e-9*db.bn/sqrt( mu_0*1.e6*df.np * m_p )
	df['zmm'] = df.zmr**2 + df.zmt**2 + df.zmn**2

	df['sig_c'] = ( df.zpm - df.zmm )/( df.zpm + df.zmm )

	# Sector rectification of cross helicity (flip sign if B_r >0 )

	for i in range(len(df.sig_c)):

		if(df.br[i] >0):

			df.sig_c[i] = -df.sig_c[i]

	# Resample the entire dataset to 1Hr period

	dfr = df.resample( '3600S' ).median( )

	# Assign the UTC timezone to the resample dataframe

	dfr.index=dfr.index.tz_localize('UTC')

	dfr.insert(0,'ssepoch',pd.Series(dfr.index,index=dfr.index).apply(datetime.datetime.timestamp))

	# Assign the name of the file in which data is to be saved

	fn = save_dir + f[-51:-8] + '_v01.hf'

#	dfr.to_pickle( fn, protocol=2 )

	hdf=hf.File(fn,'w')
	hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
	pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
	for i in dfr.columns[1:]: hdf.create_dataset(i,data=np.array(dfr[i]))
	hdf.close()

	# Append the dataframe to the list

	ldf.append( dfr )

	count += 1
	print( 'Data saved for %s' %( f[-51:] ) )
	print( count )

# Save all dataframes to one dataframe and save it in a new file

df_t = pd.concat( ldf )

fn = '/home/ahmadr/SpacePlasFest/msc/data/merged_1hr/v01/' + fnames[min_ind][-51:-8] + '_' + fnames[max_ind][-16:-8] + '_v01.hf'

#dfr.to_pickle( fn, protocol=2 )

hdf=hf.File(fn,'w')
hdf.create_dataset('datetime',data=(pd.Series(df_t.index)-\
pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
for i in df_t.columns[1:]: hdf.create_dataset(i,data=np.array(df_t[i]))
hdf.close()

print( 'Data saved to file %s' %( fn ) )

print( 'It took %s seconds to run the code' %( time.time() - start ) )
