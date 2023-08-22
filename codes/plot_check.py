import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_check(scaled=False, key_1='sc_r_median', key_2='bm_median', fn1=None, fn2=None):

    if fn1 is None:
        if scaled is True:
            fn1 = '/media/cephadrius/endless/udel_research/msc/data/binned/scaled/v14/psp_coho1hr_merged_mag_plasma_20180101_20210701_v06_76_binned_scaled_76_binned_v14.p'
        else:
            fn1 = '/media/cephadrius/endless/udel_research/msc/data/binned/unscaled/v14/psp_coho1hr_merged_mag_plasma_20180101_20210701_v06_76_binned_v14.p'
    else:
        fn1 = fn1

    if fn2 is None:
        if scaled is True:
            fn2 = '/media/cephadrius/endless/udel_research/msc/data/binned/scaled/v13/parker_coho1hr_merged_mag_plasma_20180101_20200101_v04_76_binned_scaled_v13.p'
        else:
            fn2 = '/media/cephadrius/endless/udel_research/msc/data/binned/unscaled/v13/parker_coho1hr_merged_mag_plasma_20180101_20200101_v04_76_binned_v13.p'
        
    else:
        fn2 = fn2

    df1 = pd.read_pickle(fn1)
    df2 = pd.read_pickle(fn2)

    #plt.close('all')
    plt.clf()
    #plt.figure(figsize=(8,4))
    # add a subpplot
    ax1 = plt.subplot(121)
    ax1.plot(df1[key_1], df1[key_2], 'r.', ms=3, label='v14')
    ax1.plot(df2[key_1], df2[key_2], 'b.', ms=3, label='v13')

    ax1.legend()
    ax1.set_xlabel('r (AU)', fontsize=14)
    ax1.set_ylabel(f'{key_2}', fontsize=14)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.05, 2)

    axs = plt.subplot(122)
    axs.plot(df1[key_1], df1[key_2[:-6]+'num'], 'g.', ms=3, label='v14_num')
    axs.plot(df2[key_1], df2[key_2[:-6]+'num'], 'm.', ms=3, label='v13_num')
    axs.legend()
    axs.set_xlabel('r (AU)', fontsize=14)
    # Set y label on right side
    axs.set_ylabel(f'{key_2[:-6]}num', fontsize=14)
    axs.yaxis.tick_right()
    axs.yaxis.set_label_position("right")
    axs.set_xscale('log')
    axs.set_yscale('log')
    axs.set_xlim(0.05, 2)
    plt.tight_layout()
    plt.show()

code_inputs = {
    'scaled': True,
    'key_1': 'sc_r_median',
#    'key_2': 'Tp_median',
    'fn1': None,
    'fn2': None
}

key_list = ['bm_median', 'np_median', 'Tp_median', 'vA_median', 'particle_flux_median']

plt.close('all')
for key in key_list:
    #plt.figure()
    plt.figure(figsize=(8,4))
    plot_check(**code_inputs, key_2=key)