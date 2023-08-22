#from spacepy.pycdf import CDF as cdf
from glob import glob
import pandas as pd
import h5py as hf
import time as tm
import numpy as np
import os

# Suppress warnings
import warnings
warnings.filterwarnings('ignore')

#from fitting_func import fit_func, trace_func

start = tm.time()


def trace_func(df=None, key_list=None) :
    '''The following function computes the value of each parameter and scales it against the value
    of the same parameter at 1 AU using an interpolation  technique and using Omni data and returns
    a Pandas DataFrame
    '''
    import pandas as pd
    import pytz
    from scipy.interpolate import interp1d as spl

    # Set the timezone to UTC for df
    df.index = df.index.tz_localize(pytz.utc)

    au = 1.496e8 # 1 AU in kms
    mu_0 = 4 * np.pi * 1.e-7
    kb = 1.38e-23

    df_omni_o = pd.read_pickle('/mnt/cephadrius/udel_research/msc/omni/data/processed/v03/omni_coho1hr_merged_mag_plasma_19630101_20211201_v03.p')

    # For any key in df_omni_o, if the value is smaller/larger than -/+1e-20, then set it to NaN
    for key in df_omni_o.keys() :
            df_omni_o[key][(df_omni_o[key] < -1e20) | (df_omni_o[key] > 1e20)] = np.nan

    df_omni_o['particle_flux'] = (df_omni_o.np*1e6) * (df_omni_o.vp_m*1e3) * ((1.49e11)**2)

    # Compute the proton beta
    df_omni_o["proton_beta"] = 1.e6 * df_omni_o.np * kb * df_omni_o.Tp * 2 * mu_0 / (
                               1e-9 * df_omni_o.bm)**2

    # Compute the alfven mach numver
    df_omni_o['alfven_ratio'] = 1.e3*df_omni_o.vp_m/df_omni_o.vA
 
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
    ##for somethingh in range(1):
    #    # If all the velocity values in the data is NaN, then we fall back to
    #    # using the velocity from OMNI data. We compute the minimum and maximum
    #    # time values for the spacecraft observation, and using that, we find
    #    # out the average solar wind speed during that time at 1 AU. We use that
    #    # solar wind speed to find the approximate time. A better method would
    #    # be to find the approximate time at 1 AU instead of just using the
    #    # time of spacecraft. Will have to implement that part in a later
    #    # version
    #    # TODO: Change the way time is being used to compute the average speed
    #    t_max = df.index.max()
    #    t_min = df.index.min()
    #    dfv = df_omni.vp_r[ (df_omni.index > t_min) & (df_omni.index < t_max) ]
    #    vg = dfv.groupby(pd.Grouper(freq='4W')).apply(lambda x: np.nanmedian(x))
#
    #    new_v = np.zeros_like(df.Tp)
    #    new_v[:] = np.nan
    #    ind1 = 0
#
    #    for i in range(len(vg)) :
    #        ind2 = max(np.where(df.index < vg.index[i])[0])
    #        new_v[ind1:ind2] = vg[i]
    #        ind1 = ind2
    #    t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1)*au/new_v)
    ##else :
    t_omni_sc = np.array(t_sc_unix) - np.array((df.sc_r - 1) * au / df.vp_r)

    df_intrp = pd.DataFrame(index=df.index)
    df_intrp['sc_r'] = df.sc_r

    for key in key_list :
        intrp_func = spl(t_omni, df_omni[key], fill_value='extrapolate')
        data_new = intrp_func(t_omni_sc)
        df_intrp[key] = df[key]/data_new
        print('Scaling done for %s' %(key))

    return df_intrp


def bin_data(df=None, filename=None, filetype='pickle', n_bin=100, file_version='v20',
             key_list=None) :

    '''Define the function which takes a Pandas DataFrame and file name with number
    of bins to bin the data and save the files

    df = dataframe
    filename = name of the file to which binned data is to be saved
    filetype = [pickle or hdf], default is pickle
    n_bin = number of bins between the smallest and largest distance, 
    default is 100
    file_version = file version number, default is 100
    '''

    dfn= pd.DataFrame()

    for key in key_list :

        #print('Binning done for %s ' %(key))

        #dfn[key+'_mean'] = empty_like(r_bin)
        #dfn[key+'_mean'] = nan

        dfn[key+'_median'] = np.full(n_bin, np.nan)
        dfn[key+'_iqr_10'] = np.full(n_bin, np.nan)
        dfn[key+'_iqr_25'] = np.full(n_bin, np.nan)
        dfn[key+'_iqr_50'] = np.full(n_bin, np.nan)
        dfn[key+'_iqr_75'] = np.full(n_bin, np.nan)
        dfn[key+'_iqr_90'] = np.full(n_bin, np.nan)
        dfn[key+'_num'] = np.full(n_bin, np.nan)

        for ind in range(len(r_bin) - 1) :
            try :
                dat = df[key][(df.sc_r > r_bin[ind]) & 
                              (df.sc_r <= r_bin[ind+1]) &
                              (~np.isnan(df[key]))
                              ]

                #dfn[key+'_mean'][ind] = dat.mean()
                dfn[key+'_median'][ind] = dat.median()
                dfn[key+'_iqr_10'][ind] = dat.quantile(0.1)
                dfn[key+'_iqr_25'][ind] = dat.quantile(0.25)
                dfn[key+'_iqr_50'][ind] = dat.quantile(0.50)
                dfn[key+'_iqr_75'][ind] = dat.quantile(0.75)
                dfn[key+'_iqr_90'][ind] = dat.quantile(0.90)
                dfn[key+'_num'][ind] = len(dat)
            except :
                #print(key+'_mean values set to NaN')
                pass

    #print('Data binned and saved to a dataframe')

    # Assign a new name to the pickle file

    #nf = filename[:-2] + '_%s_binned_scaled.p' %(n_bin)

    if('pickle' in filetype) :
        # Save the binned data to a pickle file using protocol 2. Avoid other
        # protocols, since protocol greater than 2 doesn't work for Python 2.*

        fn = f"{filename}_{n_bin}_binned_{file_version}.p"
        pd.DataFrame.to_pickle(dfn, fn, protocol=2)
        print('Data written for Pickle file')

    if('hdf' in filetype) :

        fn = f"{filename}_{n_bin}_binned_{file_version}.hf"
        hdf=hf.File(fn,'w')
        for i in dfn.columns[0:]: hdf.create_dataset(i,data=np.array(dfn[i]))
        hdf.close()
        print('Data written for HDF file')

    #else :

    
    #raise TypeError('File type must be either hdf or pickle')

    #print('Data saved to pickle file for %s' %(nf[54:]))

    return(dfn)

r_min = -1.2
r_max = 2
n_bin = 80
r_bin = np.logspace(r_min, r_max, n_bin)

scaled = True
binning = True
fnames = np.sort(glob("/mnt/cephadrius/udel_research/msc/data/merged_1hr/v07/*.p"))

file_version = 'v20'
df_all_l = []

for f in fnames[10:11]:
    df = pd.read_pickle(f)

    if scaled:
        n_key_list = ['Tp', 'bm', 'br', 'bt', 'bn', 'np', 'sig_c', 'vp_m', 'vp_r', 'vp_t', 'vp_n',
                      'zmm', 'zmr', 'zmt', 'zmn', 'zpm', 'zpr', 'zpt', 'zpn', 'vA', 'particle_flux',
                      'proton_beta', 'alfven_ratio']

        dfr = df.copy()

        df_intrp  = trace_func(dfr, key_list=n_key_list)

        df_all_l.append(df_intrp)

        if(binning) :
            save_dir1 = "/mnt/cephadrius/udel_research/msc/data/binned/scaled/v20/"
            key_list = list(df.keys())
            df_temp = bin_data(df_intrp, filename=save_dir1 + f[54:-2]+'_%s_binned_scaled' %(n_bin),
                            filetype=['hdf', 'pickle'], n_bin=n_bin, file_version=file_version,
                            key_list=key_list)
    else:
        save_dir1 = "/mnt/cephadrius/udel_research/msc/data/binned/unscaled/v20/"
        key_list = list(df.keys())
        df_all_l.append(df)
        if binning:
            df_temp = bin_data(df, save_dir1 + f[54:-2], filetype=['hdf', 'pickle'], n_bin=n_bin,
                               file_version=file_version, key_list=key_list)

    print(f"Code finished running for file {f}")

df_all = pd.concat(df_all_l, sort=False)

save_dir = "/mnt/cephadrius/udel_research/msc/data/all_data/v20/"
version = 'v20'
key_list = list(df_all.keys())
if(scaled):
    hdf=hf.File(save_dir+'all_spcaecraft_data_scaled_%s.hf' %(version),'w')
    for i in df_all.columns[0:]: hdf.create_dataset(i,data=np.array(df_all[i]))
    hdf.close()
    pd.DataFrame.to_pickle(df_all, save_dir+'all_spcaecraft_data_scaled_%s.p' %(version),
                           protocol=2)
    dfn = bin_data(df_all, filename=save_dir+'all_spacecraft_data_binned_scaled',
                                           filetype=['hdf', 'pickle'], n_bin=n_bin,file_version=file_version, key_list=key_list)
else :
    hdf=hf.File(save_dir+'all_spcaecraft_data_%s.hf' %(version),'w')
    for i in df_all.columns[0:]: hdf.create_dataset(i,data=np.array(df_all[i]))
    hdf.close()
    pd.DataFrame.to_pickle(df_all, save_dir+'all_spcaecraft_data_%s.p' %(version), protocol=2)
    dfn = bin_data(df_all, filename=save_dir+'all_spacecraft_data_binned',
                                           filetype=['hdf', 'pickle'], n_bin=n_bin, file_version=file_version, key_list=key_list)
