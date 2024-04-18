from spacepy.pycdf import CDF as cdf
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf
import sys
import numpy as np
from pathlib import Path

import time

start = time.time()

mu_0 = 4 * np.pi * 1.0e-7
m_p = 1.67e-27

# Name of the spacecraft

msc = "omni"

# Location of unprocessed data

data_dir = "/mnt/cephadrius/udel_research/msc/omni/coho1hr_magplasma"

# Location where the processed data is to be saved

save_dir = "/mnt/cephadrius/udel_research/msc/omni/data/processed/v2024.1"

# Search and sort all the file in the data_dir

fnames = np.sort(glob(data_dir + "/**/*.cdf", recursive=True))

if len(fnames) == 0:
    print(
        "No relevant files found for spacecraft %s. Exiting execution now. Good luck debugging the code LOL"
        % (msc)
    )
    sys.exit()
else:
    print("A total of %s relevant files found for spacecraft %s" % (len(fnames), msc))

df = None

# Define a list where all the DataFrames will be stored
ldf = []
count = 0

"""
KeysView(<CDF:
ABS_B: CDF_REAL4 [744]
BN: CDF_REAL4 [744]
BR: CDF_REAL4 [744]
BT: CDF_REAL4 [744]
Epoch: CDF_EPOCH [744]
N: CDF_REAL4 [744]
T: CDF_REAL4 [744]
V: CDF_REAL4 [744]
azimuthAngle: CDF_REAL4 [744]
elevAngle: CDF_REAL4 [744]
heliographicLatitude: CDF_REAL4 [744]
heliographicLongitude: CDF_REAL4 [744]
>)
"""

for f in fnames[0:]:
    # Read all the files and save relevant parameters to a dictionary
    d = {}
    dat = cdf(f)

    d["Epoch"] = dat["Epoch"]

    d["br"] = dat["BR"]
    d["bt"] = dat["BT"]
    d["bn"] = dat["BN"]

    d["bm"] = dat["ABS_B"]

    d["azimuthAngle"] = dat["azimuthAngle"]
    d["elevAngle"] = dat["elevAngle"]

    d["heliographicLatitude"] = dat["heliographicLatitude"]
    d["heliographicLongitude"] = dat["heliographicLongitude"]

    d["vp_m"] = dat["V"]
    d["np"] = dat["N"]
    d["Tp"] = dat["T"]

    # Change the dictionary to a dataframe
    # This is done because further processing of data like resampling and stuff
    # is easier and faster in a dataframe compared to a simple dictionary and/or
    # arrays

    df = pd.DataFrame(d, index=d["Epoch"])

    # Replace all the bad datapoints with NaN

    # df[(df < -9.99e30) | (df > 9.99e30)] = NaN
    df.replace(-1.0e31, np.nan, inplace=True)
    df["vp_r"] = (
        df.vp_m * np.cos(np.deg2rad(df.elevAngle)) * np.cos(np.deg2rad(df.azimuthAngle))
    )
    df["vp_t"] = (
        df.vp_m * np.cos(np.deg2rad(df.elevAngle)) * np.sin(np.deg2rad(df.azimuthAngle))
    )
    df["vp_n"] = df.vp_m * np.sin(np.deg2rad(df.elevAngle))

    # Create a partial datafgrame with just a few parameters (magnetic field
    # and density with velocity) which are important for the calculation of
    # Elsasser variables and cross helicity
    dfd = pd.DataFrame(
        {
            "br": df.br,
            "bt": df.bt,
            "bn": df.bn,
            "vp_r": df.vp_r,
            "vp_t": df.vp_t,
            "vp_n": df.vp_n,
        }
    )

    # Set the index of the partial dataframe to the same as the main dataframe
    dfd.index = df.index

    # Compute the deviation from mean for all the parameters from partial
    # daraframe. NOTE: Mean is calculated over '24Hours' in this case. This was
    # done so that we always have multiple correlation time lengths and the
    # deviation makes sense from a statistical point

    db = dfd.groupby(pd.Grouper(freq="24h")).apply(lambda x: x - x.mean())

    # Set the index of the deviation dataframe to the same as the main dataframe
    db.index = df.index

    # Compute the Elsasser variables and the normalized cross helicity

    df["zpr"] = 1.0e3 * db.vp_r + 1e-9 * db.br / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zpt"] = 1.0e3 * db.vp_t + 1e-9 * db.bt / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zpn"] = 1.0e3 * db.vp_n + 1e-9 * db.bn / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zpm"] = df.zpr**2 + df.zpt**2 + df.zpn**2

    df["zmr"] = 1.0e3 * db.vp_r - 1e-9 * db.br / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zmt"] = 1.0e3 * db.vp_t - 1e-9 * db.bt / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zmn"] = 1.0e3 * db.vp_n - 1e-9 * db.bn / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
    df["zmm"] = df.zmr**2 + df.zmt**2 + df.zmn**2

    df["sig_c"] = (df.zpm - df.zmm) / (df.zpm + df.zmm)

    # Sector rectification of cross helicity (flip sign if B_r >0)

    # for i in range(len(df.sig_c)):
    # If df.br is greater than 0, then flip the sign of df.sig_c
    df.loc[df.br > 0, "sig_c"] = -df.loc[df.br > 0, "sig_c"]

    # Resample the entire dataset to 1Hr period

    dfr = df.resample("3600s").median()

    # Assign the UTC timezone to the resample dataframe

    dfr.index = dfr.index.tz_localize("UTC")

    dfr.insert(
        0,
        "ssepoch",
        pd.Series(dfr.index, index=dfr.index).apply(datetime.datetime.timestamp),
    )

    # Add Alfven velocity to Omni data

    dfr["vA"] = 1e-9 * dfr.bm / np.sqrt(mu_0 * 1.0e6 * dfr.np * m_p)

    # Assign the name of the file in which data is to be saved

    fn = save_dir + "/indiv_files/" + Path(f).name[:-8] + "_v2024.1.p"
    fn_hf = save_dir + "/indiv_files/" + Path(f).name[:-8] + "_v2024.1.hf"

    dfr.to_pickle(fn)

    hdf = hf.File(fn_hf, "w")
    hdf.create_dataset(
        "datetime",
        data=(
            pd.Series(dfr.index) - pd.to_datetime("1970-01-01 00:00:00", utc=True)
        ).dt.total_seconds(),
    )

    time_keys = ["ssepoch", "datetime", "Epoch"]

    for i in dfr.columns[:]:
        if i not in time_keys:
            hdf.create_dataset(i, data=np.array(dfr[i]))
    hdf.close()

    # hdf=hf.File(fn,'w')
    # hdf.create_dataset('datetime',data=(pd.Series(dfr.index)-\
    # pd.datetime(1970,1,1,0,0,0,0, pytz.UTC)).dt.total_seconds())
    # for i in dfr.columns[1:]: hdf.create_dataset(i,data=np.array(dfr[i]))
    # hdf.close()

    # Append the dataframe to the list

    ldf.append(dfr)

    count += 1
    # print('Data saved for %s' %(f[-45:]))
    print(count)

# Save all dataframes to one dataframe and save it in a new file

df_t = pd.concat(ldf)

fn = (
    save_dir
    + "/"
    + Path(fnames[0]).name[:-8]
    + "_"
    + Path(fnames[-1]).name[-16:-8]
    + "_v2024.1.p"
)
fn_hf = (
    save_dir
    + "/"
    + Path(fnames[0]).name[:-8]
    + "_"
    + Path(fnames[-1]).name[-16:-8]
    + "_v2024.1.hf"
)

# If the files fn and fn_hf already exist, then delete them
if Path(fn).exists():
    Path(fn).unlink()
if Path(fn_hf).exists():
    Path(fn_hf).unlink()


df_t.to_pickle(fn)

hdf = hf.File(fn_hf, "w")
hdf.create_dataset(
    "datetime",
    data=(
        pd.Series(df_t.index) - pd.to_datetime("1970-01-01 00:00:00", utc=True)
    ).dt.total_seconds(),
)

time_keys = ["ssepoch", "datetime", "Epoch"]

for i in df_t.columns[:]:
    if i not in time_keys:
        hdf.create_dataset(i, data=np.array(df_t[i]))
hdf.close()

print("Data saved to file %s" % (fn))

print("It took %s seconds to run the code" % (time.time() - start))
