from cProfile import label
import pandas as pd
import numpy as np
import importlib
import glob
import matplotlib.pyplot as plt
#from data_binning import trace_func as tfnc
#importlib.reload(tfnc)


def trace_func(df=None, key_list=None) :
    '''The following function computes the value of each parameter and scales it against the value
    of the same parameter at 1 AU using an interpolation  technique and using Omni data and returns
    a Pandas DataFrame
    '''
    import pandas as pd
    import pytz
    from scipy.interpolate import interp1d as spl

    # Set the timezone to UTC for df
    try:
        df.index = df.index.tz_localize(pytz.utc)
    except:
        pass

    au = 1.496e8 # 1 AU in kms

    df_omni_o = pd.read_pickle('/mnt/cephadrius/udel_research/msc/omni/data/processed/v03/omni_coho1hr_merged_mag_plasma_19630101_20211201_v03.p')

    # For any key in df_omni_o, if the value is smaller/larger than -/+1e-20, then set it to NaN
    for key in df_omni_o.keys() :
            df_omni_o[key][(df_omni_o[key] < -1e20) | (df_omni_o[key] > 1e20)] = np.nan

    df_omni_o['particle_flux'] = (df_omni_o.np*1e6) * (df_omni_o.vp_m*1e3) * ((1.49e11)**2)

    #df_omni = df_omni_o.resample('4W').median()
    df_omni = df_omni_o.rolling(window='28D', min_periods=100).median()

    t_omni = (pd.Series(df_omni.index) - pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds()

    try :
        t_sc_unix = ((df.index - pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).total_seconds())
        print('Computing time done')
    except :
        pass

    # Compute the time at Omni corresponding to time at the spacecraft
    #if(np.isnan(df.vp_r).all()) :
    for somethingh in range(1):
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
        vg = dfv.groupby(pd.Grouper(freq='4W')).apply(lambda x: np.nanmedian(x))

        new_v = np.zeros_like(df.Tp)
        new_v[:] = np.nan
        ind1 = 0

        for i in range(len(vg)) :
            ind2 = max(np.where(df.index < vg.index[i])[0])
            new_v[ind1:ind2] = vg[i]
            ind1 = ind2
        t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1)*au/new_v)
    #else :
    #    t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1)*au/df.vp_r)

    #print(f"Length of t_omni_sc: {len(t_omni_sc)}\n")
    #print(f"Length of df.index: {len(df.index)}\n")
    df_intrp = pd.DataFrame(index=df.index)
    df_intrp['sc_r'] = df.sc_r

    for key in key_list[0:2]:
        intrp_func = spl(t_omni, df_omni[key], fill_value='extrapolate')
        data_new = intrp_func(t_omni_sc)
        df_intrp[key] = df[key]/data_new
        print('Scaling done for %s' %(key))

    return df_intrp

"""
df_omni_o = pd.read_pickle('/mnt/cephadrius/udel_research/msc/omni/data/processed/v03/omni_coho1hr_merged_mag_plasma_19630101_20211201_v03.p')

df_omni_n = df_omni_o.copy()
for key in df_omni_o.keys() :
    df_omni_n[key][(df_omni_n[key] < -1e20) | (df_omni_n[key] > 1e20)] = np.nan
# Check all the places where values are finite for each key in df_omni_o
data_len_dict = {}
for key in df_omni_o.keys() :
    #print(key)
    #print(len(df_omni_o[key][np.isfinite(df_omni_o[key])]))
    data_len_dict[key] = len(df_omni_o[key][np.isfinite(df_omni_o[key])])
    #print('\n')

data_len_dict2 = {}
for key in df_omni_o.keys() :
    #print(key)
    #print(len(df_omni_o[key][np.isfinite(df_omni_o[key])]))
    data_len_dict2[key] = len(df_omni_n[key][(np.isfinite(df_omni_n[key])) & (~np.isnan(df_omni_n[key]))])
    #print('\n')

#for key in data_len_dict.keys() :
#    print(f"{key} : {data_len_dict[key], data_len_dict2[key]}")
"""

data_read = '1'
if data_read:
    fnames = np.sort(glob.glob("/mnt/cephadrius/udel_research/msc/data/merged_1hr/v06/*.p"))

    f_mariner = fnames[4]

    df_mariner = pd.read_pickle(f_mariner)
    n_key_list = ['Tp', 'bm', 'br', 'bt', 'bn', 'np', 'sig_c', 'vp_m', 'vp_r', 'vp_t', 'vp_n',
                          'zmm', 'zmr', 'zmt', 'zmn', 'zpm', 'zpr', 'zpt', 'zpn', 'vA', 'particle_flux']

    # Find the number of unique values for key sc_r
    #print(df_mariner.sc_r.unique())
    # Select df where time is between '2018-10-01' and '2018-11-30'
    dfr = df_mariner.copy()
    dfr = dfr[(dfr.index >= '2018-10-01') & (dfr.index <= '2018-11-30')]

    df_intrp  = trace_func(dfr, key_list=n_key_list)
    # Check all the places where values are finite for each key in df_mariner

    data_len_dict_mariner = {}
    data_len_dict_intrp = {}
    for key in df_intrp.keys() :
        #print(key)
        #print(len(df_mariner[key][np.isfinite(df_mariner[key])]))
        data_len_dict_mariner[key] = len(df_mariner[key][(np.isfinite(df_mariner[key])) & (~np.isnan(df_mariner[key]))])
        data_len_dict_intrp[key] = len(df_intrp[key][(np.isfinite(df_intrp[key]) & (~np.isnan(df_intrp[key])))] )
        #print('\n')

    for key in data_len_dict_mariner.keys() :
        print(f"{key} : {data_len_dict_mariner[key], data_len_dict_intrp[key]}")


plt.close("all")
plt.figure()
plt.plot(df_mariner.index, df_mariner.bm, 'rd', ms=2, alpha=1, label='B_mag_mariner')
plt.yscale('log')
plt.legend(loc=2)
# Rotate x-axis labels by 45 degrees
plt.xticks(rotation=45)

# Set the x limits in the plot between two dates
plt.xlim(pd.Timestamp('2018-10-01'), pd.Timestamp('2018-12-30'))
axs = plt.twinx()
axs.plot(df_intrp.index, df_intrp.bm, 'b.', ms=1, alpha=0.8, label='B_mag_intrp')
axs.plot(df_omni_n.index, df_omni_n.bm, 'g.', ms=1, alpha=0.8, label='B_mag_omni')
axs.set_yscale('log')
axs.set_xlim(pd.Timestamp('2018-10-01'), pd.Timestamp('2018-11-30'))

plt.legend(loc=1)
#plt.show()
plt.savefig('../figures/B_mag_mariner_vs_intrp_v2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)
#df_new_intrp = trace_func(dfr, key_list=n_key_list)


plt.figure()
plt.plot(df_mariner.index, df_mariner.Tp, 'rd', ms=2, alpha=1, label='Tp_mag_mariner')
plt.yscale('log')
plt.legend(loc=2)

# Rotate x-axis labels by 45 degrees
plt.xticks(rotation=45)

# Set the x limits in the plot between two dates
plt.xlim(pd.Timestamp('2018-10-01'), pd.Timestamp('2018-11-30'))
axs = plt.twinx()
axs.plot(df_intrp.index, df_intrp.Tp, 'b.', ms=1, alpha=0.8, label='Tp_mag_intrp')
axs.plot(df_omni_n.index, df_omni_n.Tp, 'g.', ms=1, alpha=0.8, label='Tp_mag_omni')
axs.set_yscale('log')
plt.legend(loc=1)
axs.set_xlim(pd.Timestamp('2018-10-01'), pd.Timestamp('2018-11-30'))
#plt.show()
plt.savefig('../figures/Tp_mag_mariner_vs_intrp_v2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.1)