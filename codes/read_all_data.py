from spacepy.pycdf import CDF as cdf
import numpy as np
from glob import glob
import pandas as pd
import datetime
import pytz
#import h5py as hf
import time as tm

start = tm.time()


mu_0 = 4 * np.pi * 1.e-7
m_p = 1.67e-27

#

fnames = np.sort( glob('/media/cephadrius/endless/udel_research/msc/data/merged_1hr/v06/*.p'))

fnames = list( fnames )

# Define the length and size of bins

r_min = -1
r_max = 2
r_bin = 76
r_bin = np.logspace( r_min, r_max, r_bin )

# Define a custom list of keys/parameters for which you want the binning to be
# done, since it doesn't make much sense to bin the data for parameters like
# time and elevation angle etc

key_list = ['br', 'bt', 'bn', 'bm', 'sc_r','np', 'Tp', 'vA', 
            'vp_m','vp_r', 'vp_t', 'vp_n', 'heliographicLatitude', 'heliographicLongitude',
            'zpr', 'zpt', 'zpn', 'zpm', 'zmr', 'zmt', 'zmn', 'zmm', 'sig_c', 'alfven_ratio',
            'parker_angle', 'particle_flux', 'proton_beta']

ldf = []
sc_n = []

for f in fnames[-4:-3] :

    #d = {}

    #dat = hf.File( f )
    sc_n.append( f[64:-50])

    #for key in dat.keys() :

    #	d[key] = dat[key]

    #df = pd.DataFrame( d )

    #df.index = pd.to_datetime( df.datetime, unit='s' )

    #df['vA']    = empty_like( df.np )
    #df['vA'][:] = nan

    #try :

    #	df['vA'] = 1e-9*df.bm/sqrt( mu_0*1.e6*df.np * m_p )

    #except :

    #	pass

    df = pd.read_pickle( f )


    print( 'File read and assigned a Dataframe for %s' %( f[36:] ) )

    dfn = None
    dfn= pd.DataFrame()

    for key in df.keys() :

        print( 'Binning date for %s ' %( key ) )

        dfn[key+'_mean'] = np.empty_like( r_bin )
        dfn[key+'_mean'] = np.nan

        dfn[key+'_median'] = np.empty_like( r_bin )
        dfn[key+'_median'] = np.nan

        for ind in range( len( r_bin ) - 1 ) :

            try :

                dfn[key+'_mean'][ind] = df[key][ (df.sc_r > r_bin[ind] ) & ( df.sc_r <= r_bin[ind+1] ) ].mean( )
                dfn[key+'_median'][ind] = df[key][ (df.sc_r > r_bin[ind] ) & ( df.sc_r <= r_bin[ind+1] ) ].median( )

            except :

                    print( key+'_mean values set to NaN' )

    print( 'Data binned and saved to a dataframe' )

    # Assign a new name to the pickle file

    nf = f[:-2] + '%s_binned.p' %( r_bin )

    # Save the binned data to a pickel file using protocol 2. Avoid other
    # protocols, since protocol greater than 2 doesn't work for Python 2.*

    pd.DataFrame.to_pickle( dfn, nf, protocol=2 )

    ldf.append( dfn )

    print( 'Data saved to pickel file for %s' %( nf[36:] ) )

df_t = pd.concat( ldf )

print( 'It took %s seconds to run this code' %( round( tm.time() - start, 2 ) ) )
