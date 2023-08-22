from spacepy.pycdf import CDF as cdf
import numpy as np
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf
import time as tm

start = tm.time()

mu_0 = 4 * np.pi * 1.e-7
m_p = 1.67e-27
kb = 1.38e-23

#
data_dir = '../data/merged_1hr/v01/'
save_dir = '../data/merged_1hr/v08/'

fnames = np.sort(glob(data_dir + '*.hf'))

fnames = list(fnames)

# Define a custom list of keys/parameters for which you want the binning to be
# done, since it doesn't make much sense to bin the data for parameters like
# time and elevation angle etc

key_list = [ 'br', 'bt', 'bn', 'bm', 'sc_r','np', 'Tp', 'vA', 
             'vp_m','vp_r', 'vp_t', 'vp_n', 'heliographicLatitude', 
             'heliographicLongitude', 'alfven_ratio', 'parker_angle',
             'particle_flux', 'zpr', 'zpt', 'zpn', 'zpm', 'zmr', 
             'zmt', 'zmn', 'zmm', 'sig_c', 'proton_beta' ]

ldf = []
sc_n = []

for f in fnames[:] :

    print(f)
    d = {}

    dat = hf.File(f)

    sc_n.append(f[23:-51])

    for key in dat.keys() :

        d[key] = dat[key]

    df = pd.DataFrame(d)

    df.index = pd.to_datetime(df.datetime, unit='s')

    try :

        df[df.np<=0] = np.nan

    except :

        pass

    if('vA' in df.keys()) :

        print('Alfven velocity already defined')

    else :

        df['vA']    = np.empty_like(df[df.columns[0]])
        df['vA'][:] = np.nan

        try :

            df['vA'] = 1e-9*df.bm/np.sqrt(mu_0 * 1.e6 * df.np * m_p)

        except :

            pass

    try :

        df['alfven_ratio'] = 1.e3*df.vp_m/df.vA

    except :

        pass

    try :

        df['proton_beta'] = 1.e6*df.np*kb*df.Tp*2*mu_0/(1e-9*df.bm)**2

    except :

        pass

    try :

        # Adding Parker Angle to the dataframe (units: radians)

        bbm = df.br.values**2 + df.bt.values**2 + df.bn.values**2
        sbm = np.array([np.sqrt(xx) for xx in bbm ])

        df['parker_angle'] = np.arccos(abs(df.br)/sbm)

    except :

        pass

    try:

        # Adding particle flux to the dataframe (units: 1/second)

        df['particle_flux'] = (df.np*1e6) * (df.vp_m*1e3) * (
                                             (df.sc_r*1.49e11)**2)

    except :

        pass

    # Set the value of rest of the keys to np.nan

    for key in key_list :

        if(key in df.keys()) :

            print(key + ' already present in the file')

        else :

            df[key] = np.empty_like(df[df.columns[0]])
            df[key] = np.nan

            print(key + ' value set to np.nan for file %s' %(f))

    print('Saving data to a pickle file for all keys')

    # Assign a new name to the pickle file

    # For omni: [56:-7]
    # For ulysses: [59:-7]
    # For all other spacecrafts: [51:-7]
    nf = save_dir + f[23:-7] + '_v08.p'

    # Save the binned data to a pickel file using protocol 2. Avoid other
    # protocols, since protocol greater than 2 doesn't work for Python 2.*

    pd.DataFrame.to_pickle(df, nf, protocol=2)

    print('Data saved to pickle file for %s' %(nf[51:]))

print('It took %s seconds to run this code' %(round(tm.time() - start, 2)))