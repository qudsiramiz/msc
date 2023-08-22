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
au = 1.49e8 # AU in kms

# Name of the spacecraft

msc = 'psp'

# Location of unprocessed data

data_dir =  '/home/ahmadr/SpacePlasFest/msc/data/native_cadence/v00/'

# Location where the processed data is to be saved

save_dir = '/home/ahmadr/SpacePlasFest/msc/data/native_cadence/v01/'

# Search and sort all the file in the data_dir

fnames = sort( glob( data_dir + '/psp*.p' ) )

print( 'A total of %s relevant files found for spacecraft %s' %( len( fnames ), msc ) )

rr = input( 'run? ==>' )

if( rr ) :

    df = None

    # Define a list where all the DataFrames will be stored

    ldf = []

    count = 0

    '''
    Index([u'ssepoch', u'vpr', u'vpt', u'vpn', u'np_moment', u'wp_moment',
    'sc_pr', u'sc_pt', u'sc_pn', u'carr_latitude', u'carr_longitude',
    'br', u'bt', u'bn', u'smooth_np', u'bmag', u'smooth_bmag', u'vA',
    'Rs', u'vpm', u'uv_b_r', u'tp', u'tp_par', u'tp_per', u'beta',
    'beta_par', u'beta_per', u'pa'],
    ='object')
    '''
    # Read all the files and save relevant parameters to a dictionary  

    d1 = pd.read_pickle( fnames[0] )

    d1 = d1.tz_localize( None )

    #d = d1[ ( d1.index > datetime.datetime( 2018, 10, 30 ) ) & ( d1.index < datetime.datetime( 2018, 11, 11 ) ) ]
    d = d1

    df = pd.DataFrame( )

    df['br'] = d.br
    df['bt'] = d.bt
    df['bn'] = d.bn
    df['bm'] = d.bmag

    df['azimuthAngle'] = empty_like( d.sc_pr )
    df['elevAngle']    = empty_like( d.sc_pr )

    df['azimuthAngle'][:] = nan
    df['elevAngle'][:]    = nan

    df['sc_r'] = sqrt( d.sc_pr**2 + d.sc_pt**2 + d.sc_pn**2 )/au


    df['heliographicLatitude']  = d.carr_latitude
    df['heliographicLongitude'] = d.carr_longitude

#	df['heliographicLatitude'][:]  = nan
#	df['heliographicLongitude'][:] = nan

    df['vp_m'] = d.vpm
    df['vp_r'] = d.vpr
    df['vp_t'] = d.vpt
    df['vp_n'] = d.vpn

    df['np'] = d.np_moment
    df['Tp'] = d.tp

    df['vA'] = 1.e3*d.vA

    # Create a partial datafgrame with just a few parameters (magnetic field
    # and density with velocity) which are important for the calculation of
    # Elsasser variables and cross helicity

    dfd = pd.DataFrame( {'br':df.br, 'bt':df.bt, 'bn':df.bn, 'np':df.np,
                         'vp_r':df.vp_r, 'vp_t':df.vp_t, 'vp_n':df.vp_n} )

# Compute the deviation from mean for all the parameters from partial
# daraframe. NOTE: Mean is calculated over '24Hours' in this case. This was
# done so that we always have multiple correlation time lengths and the
# deviation makes sense from a statistical point
                                                                            
db = dfd.groupby( pd.Grouper(freq='24H' ) ).apply( lambda x: x-x.median() )
                                                                            
#dfn = pd.merge_asof( df, db, left_index=True, right_index=True )

# Compute the Elsasser variables and the normalized cross helicity

df[df.np<=0] = NaN
#db[df.np==0] = NaN

dfn = pd.DataFrame()

df['zpr'] = 1.e3*db.vp_r.values + 1e-9*db.br.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zpt'] = 1.e3*db.vp_t.values + 1e-9*db.bt.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zpn'] = 1.e3*db.vp_n.values + 1e-9*db.bn.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zpm'] = df.zpr**2 + df.zpt**2 + df.zpn**2

df['zmr'] = 1.e3*db.vp_r.values - 1e-9*db.br.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zmt'] = 1.e3*db.vp_t.values - 1e-9*db.bt.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zmn'] = 1.e3*db.vp_n.values - 1e-9*db.bn.values/sqrt( mu_0*1.e6*df.np.values * m_p )
df['zmm'] = df.zmr**2 + df.zmt**2 + df.zmn**2

sig_c = ( df.zpm - df.zmm )/( df.zpm + df.zmm )

# Implement sector correction to cross-helicity (flip sign if Br >0)

df['sig_c'] = -sign(df.br)*sig_c

'''
# Compute the Elsasser variables and the normalized cross helicity without
# taking the mean difference
                                                                            
#df['zpra'] = 1.e3*dfd.vp_r + 1e-9*dfd.br/sqrt( mu_0*1.e6*df.np * m_p )
#df['zpta'] = 1.e3*dfd.vp_t + 1e-9*dfd.bt/sqrt( mu_0*1.e6*df.np * m_p )
#df['zpna'] = 1.e3*dfd.vp_n + 1e-9*dfd.bn/sqrt( mu_0*1.e6*df.np * m_p )
#df['zpma'] = df.zpra**2 + df.zpta**2 + df.zpna**2
#
#df['zmra'] = 1.e3*dfd.vp_r - 1e-9*dfd.br/sqrt( mu_0*1.e6*df.np * m_p )
#df['zmta'] = 1.e3*dfd.vp_t - 1e-9*dfd.bt/sqrt( mu_0*1.e6*df.np * m_p )
#df['zmna'] = 1.e3*dfd.vp_n - 1e-9*dfd.bn/sqrt( mu_0*1.e6*df.np * m_p )
#df['zmma'] = df.zmra**2 + df.zmta**2 + df.zmna**2
#                                                                            
#df['sig_ca'] = ( df.zpma - df.zmma )/( df.zpma + df.zmma )

'''

# Resample the entire dataset to 1Hr period

dfr = df.resample( '3600S' ).median( )

# Assign the UTC timezone to the resample dataframe

dfr.index=dfr.index.tz_localize('UTC')

dfr.insert(0,'ssepoch',pd.Series(dfr.index,index=dfr.index).apply(datetime.datetime.timestamp))

# Assign the name of the file in which data is to be saved

fn = save_dir + fnames[0][-46:-1] + '_v01.hf'

#dfr.to_pickle( fn, protocol=2 )

hdf=hf.File(fn,'w')
hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
for i in dfr.columns[0:]: hdf.create_dataset(i,data=np.array(dfr[i]))
hdf.close()

