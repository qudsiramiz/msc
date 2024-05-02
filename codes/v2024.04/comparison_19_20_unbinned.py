import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fn_19 = "../../data/all_data/merged_1hr/v07/psp_coho1hr_merged_mag_plasma_20180101_20211001_v07.p"
fn_20 = "../../data/all_data/merged_1hr/v2024.05/psp_coho1hr_merged_mag_plasma_20180101_20231001_v2024.05.p"

df_19 = pd.read_pickle(fn_19)
df_20 = pd.read_pickle(fn_20)

key = "bm"
plt.figure(figsize=(6, 6))

# If sc_r is smaller than 1e-10, set it to NaN.
# indx = np.where(df_20["sc_r"][:] < 0)[0]
# df_20["sc_r"][indx] = np.nan

for key in df_19.keys():
    if key not in df_20.keys():
        print(f"{key} not in df_20")
        continue
    plt.figure(figsize=(6, 6))
    plt.plot(df_19["datetime"][:], df_19[f"{key}"][:], "o", label="v19", ms=7)
    plt.plot(
        df_20["datetime"][:],
        df_20[f"{key}"][:],
        "o",
        label="v2024.05",
        alpha=0.5,
        ms=5,
    )
    # plt.yscale("log")
    # On the twinx axis, plot the sc_r data.
    # plt.twinx()
    # plt.plot(df_19["sc_r"][:], df_19[f"{key}"][:], "d", color="c", ms=10)
    # plt.plot(df_20["sc_r"][:], df_20[f"{key}"][:], "d", alpha=0.5, color="k", ms=5)
    # plt.xscale("log")
    # plt.yscale("log")
    plt.legend(loc="upper left", fontsize=12)
    # plt.text(0.9, 0.75, f"{key}", transform=plt.gca().transAxes, fontsize=15)
    plt.xlabel("Datetime", fontsize=15)
    plt.ylabel(f"{key}", fontsize=15)
    plt.savefig(f"../../figures/comparison_19_20_{key}_time_series.png")
    # plt.show()
    plt.close("all")
