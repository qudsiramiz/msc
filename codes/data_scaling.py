import numpy as np
from glob import glob
import pandas as pd
import datetime
import time

from scipy.stats import linregress as lrg
#import lmfit

#from scipy.stats import iqr

## The following function computes the value of each parameter and scales it 
## against the value of the same parameter at 1 AU using an intrapolation 
## technique and using Omni data and returns a Pandas DataFrame

def trace_func(df=None, key_list=None) :

    '''The following function computes the value of each parameter and scales it 
    against the value of the same parameter at 1 AU using an intrapolation 
    technique and using Omni datai and returns a Pandas DataFrame
    '''
    import pandas as pd

    from scipy.interpolate import interp1d as spl

    au = 1.496e8 # 1 AU in kms

    df_omni_o = pd.read_pickle('/content/gdrive/My Drive/Studies/Research/Active_Research/MSC/data/omni/omni_coho1hr_merged_mag_plasma_19630101_20191201_updated.p')

    df_omni_o['particle_flux'] = (df_omni_o.np * 1e6) * (
                      df_omni_o.vp_m * 1e3) * ((1.49e11)**2)

    #df_omni = df_omni_o.resample('4W').median()
    df_omni = df_omni_o.rolling(window='28D', min_periods=100).median()


    t_omni = (pd.Series(df_omni.index) - pd.datetime(1970,1,1,0,0,0,0)).dt.total_seconds()

    try :

        t_sc_unix = ((df.index - pd.datetime(1970,1,1,0,0,0,0)).total_seconds())

        print('Computing time done')

    except :

        pass

    # Compute the time at Omni corresponding to time at the spacecraft

    if(np.isnan(df.vp_r).all()) :

        # If all the velocity values in the data is NaN, then we fall back to
        # using the velocity from OMNI data. We compute the minimum and maximum
        # time values for the spacecraft observation, and using that, we find
        # out the average solar wind speed during that time at 1 AU. We use that
        # solar wind speed to find the approximate time. A better method would
        # be to find the approximate time at 1 AU instead of just using the
        # time of spacecraft. Will have to implement that part in a later
        # version

        # TODO: Change the way time is being used to compute the average speed

        t_max = df.index.max()
        t_min = df.index.min()
        dfv = df_omni.vp_r[ (df_omni.index > t_min) & (df_omni.index < t_max) ]
        vg = dfv.groupby(pd.Grouper(freq='4W')).apply(lambda x: x.median())

        new_v = np.zeros_like(df.Tp)
        new_v[:] = np.nan
        ind1 = 0

        for i in range(len(vg)) :

            ind2 = max(np.where(df.index < vg.index[i])[0])

            new_v[ind1:ind2] = vg[i]
            ind1 = ind2

        t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1)*au/new_v)

    else :

        t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1)*au/df.vp_r)

    df_intrp = pd.DataFrame(index=df.index)

    df_intrp['sc_r'] = df.sc_r

    for key in key_list :

        intrp_func = spl(t_omni, df_omni[key], fill_value='extrapolate')

        data_new = intrp_func(t_omni_sc)

        df_intrp[key] = df[key]/data_new

        print('Scaling done for %s' %(key))

    return df_intrp