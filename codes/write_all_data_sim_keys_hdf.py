from spacepy.pycdf import CDF as cdf
from numpy import *
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf
import sys

import time

start = time.time()

mu_0 = 4 * pi * 1.e-7
m_p = 1.67e-27

# Location of unprocessed data

data_dir =  '../data/merged_1hr/v07/'

# Location where the processed data is to be saved

save_dir = '../data/merged_1hr/v07/'

# Search and sort all the file in the data_dir

fnames = sort( glob( data_dir + '*.p' ) )

if( len( fnames ) == 0 ) :

    print( 'No relevant files found for spacecraft %s. Exiting execution now. Good luck debugging the code LOL' %( msc ) )

    sys.exit()

else :

    print( 'A total of %s relevant files found ' %( len( fnames )) )

df = None

# Define a list where all the DataFrames will be stored

ldf = []

count = 0

for f in fnames[0:] :

    # Read all the files and save relevant parameters to a dataframe

    df = pd.read_pickle( f )

    fn = f[:-1] + 'hf'

    hdf=hf.File(fn,'w')
#	hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
#	pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
    for i in df.columns[0:]: hdf.create_dataset(i,data=np.array(df[i]))
    hdf.close

    print( 'Data saved to file %s' %( fn ) )

print( 'It took %s seconds to run the code' %( time.time() - start ) )
