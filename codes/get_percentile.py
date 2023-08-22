import numpy as np
import pandas as pd
import pickle

fname = '/mnt/cephadrius/udel_research/msc/data/all_data/v19/all_spacecraft_data_binned_scaled_80_binned_v19.p'


sel_key_list = ['sc_r', 'bm', 'np', 'vp_m', 'Tp', 'proton_beta', 'vA', 'alfven_ratio']
all_keys = dict()
for keys in sel_key_list:
    all_keys[keys] = [f"{keys}_iqr_10", f"{keys}_iqr_25", f"{keys}_iqr_50", f"{keys}_iqr_75", f"{keys}_iqr_90"]

r_min = -1.2
r_max = 2
n_bin = 80
r_bin = np.logspace(r_min, r_max, n_bin)

df = pd.read_pickle(fname)

dfn = pd.DataFrame()
for keys in sel_key_list:
    relevant_keys = all_keys[keys]
    for key in relevant_keys:
        dfn[key] = df[key]

# Set dfn index to r_bin
dfn.index = r_bin