import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

fnames = np.sort(glob.glob('../data/all_data/v19/*data_scaled_v19.p'))
#fnames = np.array(['../data/all_data/v19/all_spcaecraft_data_v19.p',
#                   '../data/all_data/v19/all_spcaecraft_data_scaled_v19.p'])

data_read = ""
key1 = "datetime"
key2 = "alfven_ratio"
if data_read:
    df1 = pd.read_pickle(fnames[0])
    df2 = pd.read_pickle(fnames[1])
    df1 = df1.sort_index()
    df2 = df2.sort_index()
    dfr1 =df1.rolling(window='365D', min_periods=100).median()
    dfr2 =df2.rolling(window='365D', min_periods=100).median()


plt.close("all")
plt.figure()
legends = ["unscaled", "scaled"]
colors = ['r', 'b']
plt.plot(df1.index, df1[key2], ".", ms=1, alpha=0.01, c=colors[0])
plt.plot(dfr1.index, dfr1[key2], "-", lw=2, label="unscaled", alpha=1, c=colors[0])

plt.plot(df2.index, df2[key2], ".", ms=1, alpha=0.01, c=colors[1])
plt.plot(dfr2.index, dfr2[key2], "-", lw=2, label="scaled", alpha=1, c=colors[1])
plt.legend()
"""
for i, f in enumerate(fnames[:]):
    df = pd.read_pickle(f)
    df = df.sort_index()
    # Check i all values of the key are nan
    if np.all(np.isnan(df[key2])):
        continue
    else:
        dfr = df.rolling(window='365D', min_periods=100).median()
        plt.plot(dfr.index, dfr[key2], ".", ms=6, label=f[23:-50], alpha=1, c=colors[i])
        plt.plot(df.index, df[key2], ".", ms=6, label=f[23:-50], alpha=0.2, c=colors[i])
        #plt.plot(df.sc_r, df.vA/1e3, 'b.', ms=3, label='vA')
        plt.legend(legends[i])
        #axs = plt.twinx()
        #axs.plot(df.sc_r, df.alfven_ratio, 'g.', ms=4, label='ratio', alpha=0.8)
        #axs.set_yscale("log")
        #axs.set_xscale("log")
        #axs.legend(loc=2)
"""
#plt.xscale("log")
plt.yscale("log")
plt.ylim(0.10, 50)
plt.xlabel(key1)
plt.ylabel(key2)
plt.savefig(f'/home/cephadrius/Dropbox/all_spc_scaled_{key1}_{key2}.png')
