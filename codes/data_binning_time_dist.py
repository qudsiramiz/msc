from glob import glob
import pandas as pd
import time as tm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

# # Suppress warnings
# import warnings
# warnings.filterwarnings('ignore')

# Set the font style to Times New Roman
font = {'family': 'serif', 'weight': 'normal', 'size': 10}
plt.rc('font', **font)
plt.rc('text', usetex=True)

data_dir =  '../data/merged_1hr/v07/'
fnames = np.sort(glob(data_dir + '*.p'))

df_list = []
for f in fnames:
    df = pd.read_pickle(f)
    df_list.append(df)

# Concatenate all the dataframes
df = pd.concat(df_list)

rebin_data = '1'
if rebin_data:
    r_bins = np.logspace(np.log10(df.sc_r.min()), np.log10(df.sc_r.max()), 80)
    t_bins = np.linspace(df.index.min().value, df.index.max().value, 80)
    t_bins = pd.to_datetime(t_bins)

    Tp = np.full((len(r_bins), len(t_bins)), np.nan)
    Tp_num = np.full((len(r_bins), len(t_bins)), np.nan)
    for i in range(len(r_bins)-1):
        for j in range(len(t_bins)-1):
            tk = np.where((df.sc_r >= r_bins[i]) & (df.sc_r < r_bins[i+1]) & (
                           df.index >= t_bins[j]) & (
                           df.index < t_bins[j+1]))[0]
            if len(tk) > 0:
                Tp[i,j] = np.nanmedian(df.Tp[tk])
                Tp_num[i,j] = len(tk)
        print(i)

    Tp_num[Tp_num == 0] = np.nan

norm_plot_list = [True, False]

for norm_plot in norm_plot_list:
    plt.close('all')

    labelsize = 15
    tick_labelsize = 10
    tick_length = 6
    tick_width = 2
    m_tick_length = 4
    m_tick_width = 1
    m_tick_label_size = 8
    c_tick_labelsize = 10
    pad = 5


    fig = plt.figure(figsize=(6, 6), dpi=100, facecolor='w', edgecolor='gray')
    fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    axs = fig.add_subplot(1, 1, 1)

    # For each column of Tp, divide each bin by the maximum in that columns
    Tp_max = np.nanmax(Tp, axis=0)
    Tp_norm = np.full(Tp.shape, np.nan)
    for i in range(Tp.shape[1]):
        Tp_norm[:,i] = Tp[:,i] / Tp_max[i]


    if norm_plot:
        plot_data = Tp_norm
        cbar_norm = mpl.colors.LogNorm(vmin=1e-3, vmax=1)
    else:
        plot_data = Tp
        cbar_norm = mpl.colors.LogNorm(vmin=1e3, vmax=2e5)
    axs.pcolormesh(r_bins, t_bins, np.transpose(plot_data), cmap='viridis',
                   norm=cbar_norm, alpha=0.5, shading="nearest", rasterized=True)

    axs.set_xlabel(r'$r$ (AU)', fontsize=labelsize)
    axs.set_ylabel(r'Time (UTC)', fontsize=labelsize)
    axs.set_title(r'$T_p$')
    axs.set_xlim(r_bins[0], r_bins[-1])
    axs.set_ylim(t_bins[0], t_bins[-1])
    axs.set_xscale('log')
    axs.set_yscale('linear')
    axs.set_aspect('auto')
    axs.grid(True)
    #axs.set_facecolor('w')
    #axs.set_edgecolor('w')
    axs.tick_params(axis='both', which='major', labelsize=tick_labelsize)
    axs.tick_params(axis='both', which='minor', labelsize=tick_labelsize)
    axs.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    axs.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    #axs.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    #axs.yaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    axs.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=10))
    axs.yaxis.set_major_locator(mpl.ticker.LinearLocator(numticks=10))
    #axs.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=10, subs=np.arange(2,10)))
    #axs.yaxis.set_minor_locator(mpl.ticker.LinearLocator(numticks=10))
    axs.xaxis.set_tick_params(which='both', direction='in', labelsize=tick_labelsize, pad=pad)
    axs.yaxis.set_tick_params(which='both', direction='in', labelsize=tick_labelsize, pad=pad)
    axs.xaxis.set_tick_params(which='minor', direction='in', pad=pad)
    axs.yaxis.set_tick_params(which='minor', direction='in', pad=pad)
    axs.xaxis.set_ticks_position('both')
    axs.yaxis.set_ticks_position('both')
    
    # Add a colorbar to axes axs
    cbar = fig.colorbar(axs.collections[0], ax=axs, shrink=0.95, aspect=20, extend='both')
    cbar.set_label(r'$T_p$ (K)')
    cbar.ax.tick_params(labelsize=c_tick_labelsize)
    cbar.ax.yaxis.set_ticks_position('both')
    cbar.ax.yaxis.set_tick_params(which='both', direction='in', pad=10)
    cbar.ax.yaxis.set_tick_params(which='minor', direction='in', pad=10)
    cbar.ax.yaxis.set_ticks_position('both')

    # plt.savefig(f"../figures/Tp_map_median_{f.split('/')[-1].split('.')[0]}.pdf")
    plt.savefig(f"/home/cephadrius/Dropbox/Tp_map_median_all_spc_time{int(norm_plot)}.pdf",
                bbox_inches='tight', pad_inches=0.1)

