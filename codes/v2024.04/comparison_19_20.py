import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py as hf

p1_19 = "/home/vetinari/Desktop/git/msc/data/binned_unscaled_v19/v19/"
p2_19 = "_coho1hr_merged_mag_plasma_"
p3_19 = "_v07_80_binned_v19.p"

sc_fl_19 = [
    p1_19 + "psp" + p2_19 + "20180101_20211001" + p3_19,
    p1_19 + "helios1" + p2_19 + "19740101_19811201" + p3_19,
    p1_19 + "helios2" + p2_19 + "19760101_19801201" + p3_19,
    p1_19 + "uy" + p2_19 + "19900101_19920201" + p3_19,
    p1_19 + "mariner2" + p2_19 + "19620830_19621116" + p3_19,
    p1_19 + "mariner10" + p2_19 + "19731103_19740918" + p3_19,
    p1_19 + "cassini" + p2_19 + "20000101_20040101" + p3_19,
    p1_19 + "pioneer10" + p2_19 + "19720101_19950901" + p3_19,
    p1_19 + "pioneer11" + p2_19 + "19730101_19941201" + p3_19,
    p1_19 + "new_horizons" + p2_19 + "20081010_20200127" + p3_19,
    p1_19 + "voyager1" + p2_19 + "19770101_20181201" + p3_19,
    p1_19 + "voyager2" + p2_19 + "19770101_20181201" + p3_19,
]


p1_20 = "../../data/v2024.04/individual_spc/binned_unscaled/"
p2_20 = "_coho1hr_merged_mag_plasma_"
p3_20 = "_v2024.04.p"

sc_fl_20 = [
    p1_20 + "psp" + p2_20 + "20180101_20231001" + p3_20,
    p1_20 + "helios1" + p2_20 + "19740101_19811201" + p3_20,
    p1_20 + "helios2" + p2_20 + "19760101_19801201" + p3_20,
    p1_20 + "uy" + p2_20 + "19900101_19920201" + p3_20,
    p1_20 + "mariner2" + p2_20 + "19620830_19621116" + p3_20,
    p1_20 + "mariner10" + p2_20 + "19731103_19740918" + p3_20,
    p1_20 + "cassini" + p2_20 + "20000101_20040101" + p3_20,
    p1_20 + "pioneer10" + p2_20 + "19720101_19950901" + p3_20,
    p1_20 + "pioneer11" + p2_20 + "19730101_19941201" + p3_20,
    p1_20 + "newhorizons" + p2_20 + "20081010_20230731" + p3_20,
    p1_20 + "voyager1" + p2_20 + "19770101_20181201" + p3_20,
    p1_20 + "voyager2" + p2_20 + "19770101_20181201" + p3_20,
    # p1_20 + "solo" + p2_20 + "20200101_20231201" + p3_20,
]


# fn_19 = sc_fl_19[9]
# fn_20 = sc_fl_20[9]

fn_19 = "/home/vetinari/Desktop/git/msc/data/all_data/merged_1hr/psp_coho1hr_merged_mag_plasma_20180101_20211001_v01.hf"
fn_20 = "/home/vetinari/Desktop/git/msc/data/all_data/merged_1hr/psp_coho1hr_merged_mag_plasma_20180101_20231001_v01.hf"

df_19 = hf.File(fn_19, "r+")

df_20 = hf.File(fn_20, "r+")

# df_19 = pd.read_pickle(fn_19)

# df_20 = pd.read_pickle(fn_20)

arr_key = [
    "sc_r_iqr_50",
    "b_iqr_50",
    "np_iqr_50",
    "vp_iqr_50",
    "vpr_iqr_50",
    "tp_iqr_50",
    "va_iqr_50",
    "na_iqr_50",
    "loss_iqr_50",
    # "ang_iqr_50",
    "proton_beta_iqr_50",
    "alfven_ratio_iqr_50",
]
key = "bm"
plt.figure(figsize=(10, 10))

# If sc_r is smaller than 1e-10, set it to NaN.
indx = np.where(df_20["sc_r"][:] < 0)[0]
df_20["sc_r"][indx] = np.nan

"""
plt.plot(df_19["sc_r_iqr_50"], df_19[f"{key}_iqr_50"], "o", label="v19")
plt.plot(df_20["sc_r_iqr_50"], df_20[f"{key}_iqr_50"], "d", label="v2024.04", alpha=0.5)
plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper right", fontsize=15)
plt.text(0.9, 0.75, f"{key}", transform=plt.gca().transAxes, fontsize=15)
plt.savefig("../../figures/comparison_19_20.png")
# plt.close()
"""
# plt.plot(df_19["datetime"][:], df_19[f"{key}"][:], "o", label="v19", color="red", ms=10)
# plt.plot(
#     df_20["datetime"][:],
#     df_20[f"{key}"][:],
#     "o",
#     label="v2024.04",
#     alpha=0.5,
#     color="blue",
#     ms=5,
# )
# plt.yscale("log")
# # On the twinx axis, plot the sc_r data.
# plt.twinx()
plt.plot(df_19["sc_r"][:], df_19["bm"][:], "d", color="c", ms=10)
plt.plot(df_20["sc_r"][:], df_20["bm"][:], "d", alpha=0.5, color="k", ms=5)

plt.xscale("log")
plt.yscale("log")
plt.legend(loc="upper right", fontsize=15)
plt.text(0.9, 0.75, f"{key}", transform=plt.gca().transAxes, fontsize=15)
plt.savefig("../../figures/comparison_19_20.png")
plt.show()
