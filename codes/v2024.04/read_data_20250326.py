# Suppress warnings
import warnings
from glob import glob
from pathlib import Path

import h5py as hf
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

data_folder = "/mnt/cephadrius/udel_research/msc/data/v2025.03/all_data/"

scaled_file = "all_spcaecraft_data_scaled_v2025.03.p"
unscaled_file = "all_spcaecraft_data_v2025.03.p"
unscaled_binned_file = "all_spacecraft_data_85_binned_v2025.03.p"
scaled_binned_file = "all_spacecraft_data_scaled_v2025.03.p"

scaled_file_path = Path(data_folder) / scaled_file
unscaled_file_path = Path(data_folder) / unscaled_file
unscaled_binned_file_path = Path(data_folder) / unscaled_binned_file
scaled_binned_file_path = Path(data_folder) / scaled_binned_file

# Load the scaled data
df_scaled = pd.read_pickle(scaled_file_path)
df_scaled.name = "scaled_unbinned"
# Load the unscaled data
df_unscaled = pd.read_pickle(unscaled_file_path)
df_unscaled.name = "unscaled_unbinned"

# Load the unscaled binned data
df_unscaled_binned = pd.read_pickle(unscaled_binned_file_path)
df_unscaled_binned.name = "unscaled_binned"
# Load the scaled binned data
df_scaled_binned = pd.read_pickle(scaled_binned_file_path)
df_scaled_binned.name = "scaled_binned"


# Find the minimum and maximum value of sc_r, where vp_r is not NaN for all the dataframes
def find_min_max(df, key):
    min_val = df["sc_r"][~df[key].isna()].min()
    max_val = df["sc_r"][~df[key].isna()].max()
    return min_val, max_val


def print_min_max(df, key):
    min_val, max_val = find_min_max(df, key)
    print(f"{key} min: {min_val:0.3f}, max: {max_val:0.3f} for {df.name}")


print_min_max(df_scaled, "vp_r")
print_min_max(df_scaled, "vp_m")

print_min_max(df_unscaled, "vp_r")
print_min_max(df_unscaled, "vp_m")

# print_min_max(df_scaled_binned, "vp_r")
# print_min_max(df_scaled_binned, "vp_m")
#
# print_min_max(df_unscaled_binned, "vp_r")
# print_min_max(df_unscaled_binned, "vp_m")


psp_file = "/mnt/cephadrius/udel_research/msc/data/merged_1hr/v2024.05/psp_coho1hr_merged_mag_plasma_20180101_20231001_v2024.05.p"
df_psp = pd.read_pickle(psp_file)
df_psp.name = "psp_merged"
print_min_max(df_psp, "vp_r")
print_min_max(df_psp, "vp_m")

r_min = -1.4
r_max = 2

n_bin = 85
r_bin = np.logspace(r_min, r_max, n_bin)

# Select all the location where sc_r is between r_bin[2] and r_bin[3]
df_psp_selected = df_psp[(df_psp["sc_r"] > r_bin[2]) & (df_psp["sc_r"] < r_bin[3])].copy()
df_psp_selected.name = "psp_merged_selected"

selected_indices = df_psp_selected.index

# Convert these indices to a datetime object list
selected_dates = [pd.to_datetime(i) for i in selected_indices]

# Find the corresponding indices in the scaled and unscaled dataframes
scaled_indices = df_scaled.index[df_scaled.index.isin(selected_dates)]
unscaled_indices = df_unscaled.index[df_unscaled.index.isin(selected_dates)]

# Find the values of vp_r and vp_m for the selected indices in the scaled and unscaled dataframes
scaled_vp_r = df_scaled.loc[scaled_indices, "vp_r"]
scaled_vp_m = df_scaled.loc[scaled_indices, "vp_m"]
unscaled_vp_r = df_unscaled.loc[unscaled_indices, "vp_r"]
unscaled_vp_m = df_unscaled.loc[unscaled_indices, "vp_m"]
