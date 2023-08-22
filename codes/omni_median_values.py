import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle

df = pd.read_pickle(
    "../omni/data/processed/v03/omni_coho1hr_merged_mag_plasma_19630101_20211201_v03.p"
    )

# Remove all values larger than 10**30 or smaller that -10**30
for key in df.keys():
    df[key][(df[key] < -1e20) | (df[key] > 1e20)] = np.nan

# Change alfven velocity to units of km/s
df["vA"] = df["vA"] * 1e-3

# Add alfven ratio, particle flux and proton beta
mu_0 = 4 * np.pi * 1.e-7
kb = 1.38e-23

df["particle_flux"] = (df.np * 1e6) * (df.vp_m * 1e3) * ((1.49e11) ** 2)
df["proton_beta"] = (df.np * 1e6 * kb * df.Tp * 2 * mu_0) / (1e-9 * df.bm) ** 2
df["alfven_ratio"] = df.vp_m / df.vA

# Dates for Solar cycles:
c24 = ["2008-01-01", "2019-12-31"]
c23 = ["1996-05-01", "2008-01-31"]
c22 = ["1986-09-01", "1996-04-30"]
c21 = ["1976-03-01", "1986-08-31"]
c20 = ["1964-10-01", "1976-02-29"]
call = ["1964-10-01", "2019-12-31"]

# Select data from January 2008 to December 2019, corresponding to 24th Solar Cycle
df_24 = df.loc[c24[0]:c24[1]]
df_23 = df.loc[c23[0]:c23[1]]
df_22 = df.loc[c22[0]:c22[1]]
df_21 = df.loc[c21[0]:c21[1]]
df_20 = df.loc[c20[0]:c20[1]]
df_all = df.loc[call[0]:call[1]]

keys = ['br', 'bt', 'bn', 'bm', 'vp_m', 'np', 'Tp', 'vp_r', 'vp_t', 'vp_n', 'vA', 'particle_flux',
        'proton_beta', 'alfven_ratio']

# For each solar cyucle, calculate the median values of the selected keys
df_24_med = df_24[keys].median()
df_23_med = df_23[keys].median()
df_22_med = df_22[keys].median()
df_21_med = df_21[keys].median()
df_20_med = df_20[keys].median()
df_all_med = df_all[keys].median()

# Create a dataframe with the median values of each solar cycle with keys as index
df_med = pd.DataFrame([df_24_med, df_23_med, df_22_med, df_21_med, df_20_med, df_all_med], index=[24, 23, 22, 21, 20, 'all'])

# Set the index name to 'Solar Cycle'
df_med.index.name = 'Solar Cycle'

# For each key, calculate the mean values of the selected keys
df_24_mean = df_24[keys].mean()
df_23_mean = df_23[keys].mean()
df_22_mean = df_22[keys].mean()
df_21_mean = df_21[keys].mean()
df_20_mean = df_20[keys].mean()
df_all_mean = df_all[keys].mean()

# Create a dataframe with the mean values of each solar cycle with keys as index
df_mean = pd.DataFrame([df_24_mean, df_23_mean, df_22_mean, df_21_mean, df_20_mean, df_all_mean], index=[24, 23, 22, 21, 20, 'all'])

# Set the index name to 'Solar Cycle'
df_mean.index.name = 'Solar Cycle'

# Compute the 10, 25, 50, 75, 90 percentiles of the selected keys
df_24_per = df_24[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
df_23_per = df_23[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
df_22_per = df_22[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
df_21_per = df_21[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
df_20_per = df_20[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])
df_all_per = df_all[keys].quantile([0.1, 0.25, 0.5, 0.75, 0.9])


# Create a dataframe with the percentiles of each solar cycle with keys as index
#df_per = pd.DataFrame([df_24_per, df_23_per, df_22_per, df_21_per, df_20_per], index=[24, 23, 22, 21, 20])

# For each solar cycle, for each key create a dictionary with mean, median and percentiles
# and save it to a dataframe
df_24_dict = {}
df_23_dict = {}
df_22_dict = {}
df_21_dict = {}
df_20_dict = {}
df_all_dict = {}
dict_list = [df_24_dict, df_23_dict, df_22_dict, df_21_dict, df_20_dict, df_all_dict]
solar_cycle = [24, 23, 22, 21, 20, 'all']

for key in keys:
    print(key)
    print(np.round(df[key].mean(), 3))
    for i, df in enumerate([df_24, df_23, df_22, df_21, df_20, df_all]):
        dict_list[i][key] = {
            'solar_cycle': solar_cycle[i],
            "mean": np.round(df[key].mean(), 3),
            "median": np.round(df[key].median(), 3),
            "10%": np.round(df[key].quantile(0.1), 3),
            "25%": np.round(df[key].quantile(0.25), 3),
            "50%": np.round(df[key].quantile(0.5), 3),
            "75%": np.round(df[key].quantile(0.75), 3),
            "90%": np.round(df[key].quantile(0.9), 3),
        }

# Save the dictionaries to a pickle file

with open("../omni/data/processed/v03/omni_solar_cycle_average_values_v02.p", "wb") as f:
    pickle.dump(dict_list, f)

# Load the pickle file
with open("../omni/data/processed/v03/omni_solar_cycle_average_values_v02.p", "rb") as f:
    omni_av_data = pickle.load(f)