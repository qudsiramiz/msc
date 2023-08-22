import datetime
import sys
import time
from glob import glob

import h5py as hf
import numpy as np
import pandas as pd
import pytz
from numpy import *
from spacepy.pycdf import CDF as cdf

start = time.time()

mu_0 = 4 * np.pi * 1.e-7
m_p = 1.67e-27

# Name of the spacecraft

msc = 'new-horizons'

# Location of unprocessed data

data_dir =  '../' + msc + '/data/'

# Location where the processed data is to be saved

save_dir = f'../{msc}/data/processed/'

# Search and sort all the file in the data_dir

fnames = sort(glob(data_dir + '/**/*.cdf' , recursive=True))

if(len(fnames) == 0) :

    print('No relevant files found for spacecraft %s. Exiting execution now. Good luck debugging the code LOL' %(msc))

    sys.exit()

else :

    print('A total of %s relevant files found for spacecraft %s' %(len(fnames), msc))

df = None

# Define a list where all the DataFrames will be stored

ldf = []

count = 0

'''
NOTE: Updated with additional data on date: 2022-06-16 at 1935 UTC
KeysView(<CDF:
Epoch: CDF_TIME_TT2000 [239252]
Epoch_MINUS: CDF_INT8 [] NRV
Epoch_PLUS: CDF_INT8 [239252]
NH_HAE_J2000_D_LAT: CDF_REAL4 [239252]
NH_HAE_J2000_D_LON: CDF_REAL4 [239252]
NH_HAE_J2000_D_X: CDF_REAL4 [239252]
NH_HAE_J2000_D_Y: CDF_REAL4 [239252]
NH_HAE_J2000_D_Z: CDF_REAL4 [239252]
NH_HGI_D_LAT: CDF_REAL4 [239252]
NH_HGI_D_LON: CDF_REAL4 [239252]
NH_HGI_D_R: CDF_REAL4 [239252]
NH_HGI_D_X: CDF_REAL4 [239252]
NH_HGI_D_Y: CDF_REAL4 [239252]
NH_HGI_D_Z: CDF_REAL4 [239252]
NH_HG_D_LAT: CDF_REAL4 [239252]
NH_HG_D_LON: CDF_REAL4 [239252]
NH_HG_D_X: CDF_REAL4 [239252]
NH_HG_D_Y: CDF_REAL4 [239252]
NH_HG_D_Z: CDF_REAL4 [239252]
n: CDF_REAL4 [239252]
pdyn: CDF_REAL4 [239252]
pth: CDF_REAL4 [239252]
sp_met: CDF_INT4 [239252]
st_met: CDF_INT4 [239252]
t: CDF_REAL4 [239252]
v: CDF_REAL4 [239252]
>
'''

for f in fnames[1:] :

    # Read all the files and save relevant parameters to a dictionary

    d = {}
    dat = cdf(f)

    d['Epoch'] = dat['Epoch']

#	d['br'] = dat['BR']
#	d['bt'] = dat['BT']
#	d['bn'] = dat['BN']

#	d['bm'] = dat['ABS_B']

#	d['azimuthAngle'] = dat['azimuthAngle']
#	d['elevAngle']   = dat['elevAngle']

    d['sc_r'] = dat['NH_HGI_D_R']

    d['heliographicLatitude']  = dat['NH_HGI_D_LAT']
    d['heliographicLongitude'] = dat['NH_HGI_D_LON']

    d['vp_m'] = dat['v']
    d['np'] = dat['n']
    d['Tp'] = dat['t']

    # Change the dictionary to a dataframe
    # This is done because further processing of data like resampling and stuff
    # is easier and faster in a dataframe compared to a simple dictionary and/or
    # arrays

    df = pd.DataFrame(d, index=d['Epoch'])

    # Remove the unnecessary columns from the dataframe
    df.drop(['Epoch'], axis=1, inplace=True)

    # Replace all the bad datapoints with NaN
    df[(df < -9.99e30) | (df > 9.99e30)] = np.nan

#	df['vp_r'] = df.vp_m * cos(deg2rad(df.elevAngle)) * cos(deg2rad(df.azimuthAngle))
#	df['vp_t'] = df.vp_m * cos(deg2rad(df.elevAngle)) * sin(deg2rad(df.azimuthAngle))
#	df['vp_n'] = df.vp_m * sin(deg2rad(df.elevAngle))
#
#	# Create a partial datafgrame with just a few parameters (magnetic field
#	# and density with velocity) which are important for the calculation of
#	# Elsasser variables and cross helicity
#
#	dfd = pd.DataFrame({'br':df.br, 'bt':df.bt, 'bn':df.bn, 
#	                     'vp_r':df.vp_r, 'vp_t':df.vp_t, 'vp_n':df.vp_n})
#
#	# Compute the deviation from mean for all the parameters from partial
#	# daraframe. NOTE: Mean is calculated over '24Hours' in this case. This was
#	# done so that we always have multiple correlation time lengths and the
#	# deviation makes sense from a statistical point
#
#	db = dfd.groupby(pd.Grouper(freq='24H')).apply(lambda x: x-x.mean())
#
#	# Compute the Elsasser variables and the normalized cross helicity
#
#	df['zpr'] = 1.e3*db.vp_r + 1e-9*db.br/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zpt'] = 1.e3*db.vp_t + 1e-9*db.bt/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zpn'] = 1.e3*db.vp_n + 1e-9*db.bn/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zpm'] = df.zpr**2 + df.zpt**2 + df.zpn**2
#
#	df['zmr'] = 1.e3*db.vp_r - 1e-9*db.br/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zmt'] = 1.e3*db.vp_t - 1e-9*db.bt/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zmn'] = 1.e3*db.vp_n - 1e-9*db.bn/sqrt(mu_0*1.e6*df.np * m_p)
#	df['zmm'] = df.zmr**2 + df.zmt**2 + df.zmn**2
#
#	df['sig_c'] = (df.zpm - df.zmm)/(df.zpm + df.zmm)

    # Resample the entire dataset to 1Hr period

    dfr = df.resample('3600S').median()

    # Assign the UTC timezone to the resample dataframe

    dfr.index=dfr.index.tz_localize('UTC')

    dfr.insert(0,'ssepoch',pd.Series(dfr.index,index=dfr.index).apply(datetime.datetime.timestamp))

    # Assign the name of the file in which data is to be saved

    fn = save_dir + fnames[0][-52:-11] + '.hf'

    dfr.to_pickle(f"{fn[:-3]}.p", protocol=2)

    hdf=hf.File(fn,'w')
    hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
    pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
    for i in dfr.columns[1:]: hdf.create_dataset(i,data=np.array(dfr[i]))
    hdf.close()

    # Append the dataframe to the list

    ldf.append(dfr)

    count += 1
    print('Data saved for %s' %(f[-45:]))
    print(count)

# Save all dataframes to one dataframe and save it in a new file

df_t = pd.concat(ldf)

#fn = fnames[0][-52:-11] + '.hf'

dfr.to_pickle(fn, protocol=2)

hdf=hf.File(fn,'w')
hdf.create_dataset('datetime',data=(pd.Series(df_t.index)-\
pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
for i in df_t.columns[1:]: hdf.create_dataset(i,data=np.array(df_t[i]))
hdf.close()

print('Data saved to file %s' %(fn))

print('It took %s seconds to run the code' %(time.time() - start))
