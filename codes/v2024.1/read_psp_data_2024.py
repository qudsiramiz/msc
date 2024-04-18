from spacepy.pycdf import CDF as cdf
from pathlib import Path
import numpy as np
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf
import sys

import time

start = time.time()

mu_0 = 4 * np.pi * 1.0e-7
m_p = 1.67e-27

# Name of the spacecraft

msc = "psp"

# Location of unprocessed data
data_dir = f"/mnt/cephadrius/udel_research/msc/{msc}/coho1hr/"

# Location where the processed data is to be saved
save_dir = f"/mnt/cephadrius/udel_research/msc/{msc}/processed/"

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
<CDF:
B: CDF_REAL4 [744]
BN: CDF_REAL4 [744]
BR: CDF_REAL4 [744]
BT: CDF_REAL4 [744]
Epoch: CDF_EPOCH [744]
ProtonSpeed: CDF_REAL4 [744]
VN: CDF_REAL4 [744]
VR: CDF_REAL4 [744]
VT: CDF_REAL4 [744]
flow_lon: CDF_REAL4 [744]
flow_theta: CDF_REAL4 [744]
heliographicLatitude: CDF_REAL4 [744]
heliographicLongitude: CDF_REAL4 [744]
protonDensity: CDF_REAL4 [744]
protonTemp: CDF_REAL4 [744]
radialDistance: CDF_REAL4 [744]
"""

# Loop over all the files and save the relevant data in a DataFrame
for fname in fnames[:]:
    d = {}
    dat = cdf(fname)

    d["Epoch"] = dat["Epoch"][:]
    d["bm"] = dat["B"][:]
    d["br"] = dat["BR"][:]
    d["bt"] = dat["BT"][:]
    d["bn"] = dat["BN"][:]

    d["np"] = dat["protonDensity"][:]
    d["vp_m"] = dat["ProtonSpeed"][:]
    d["Tp"] = dat["protonTemp"][:]

    d["sc_r"] = dat["radialDistance"][:]

    d["heliographicLatitude"] = dat["heliographicLatitude"][:]
    d["heliographicLongitude"] = dat["heliographicLongitude"][:]

    # TODO: Check if the following is correct
    d["azimuthAngle"] = dat["flow_theta"][:]
    d["elevAngle"] = dat["flow_lon"][:]

    d["vp_r"] = dat["VR"][:]
    d["vp_t"] = dat["VT"][:]
    d["vp_n"] = dat["VN"][:]

    # Change the dictionary to a DataFrame and set the index to the Epoch
    # This is done because further processing of data like resampling and stuff is easier and faster
    # in a dataframe compared to a simple dictionary and/or arrays
    df = pd.DataFrame(d, index=d["Epoch"])

    # Replace all the bad data points with NaN
    df.replace(-1.0e31, np.nan, inplace=True)
    # df[(df < -9.99e30) | (df > 9.99e30)] = np.nan

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
    db = dfd.groupby(pd.Grouper(freq="24H")).apply(lambda x: x - x.mean())

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

    # Sector rectification of cross helicity (flip sign if B_r >0 )
    for i in range(len(df.sig_c)):
        if df.br[i] > 0:
            df.sig_c[i] = -df.sig_c[i]

    # Resample the entire dataset to 1Hr period
    dfr = df.resample("3600S").median()

    # Assign the UTC timezone to the resample dataframe
    dfr.index = dfr.index.tz_localize("UTC")

    dfr.insert(
        0,
        "ssepoch",
        pd.Series(dfr.index, index=dfr.index).apply(datetime.datetime.timestamp),
    )

    # Get the name of the file from "fname"
    file_name = Path(fname).name

    # Add "save_dir" to the file name
    fn = save_dir + file_name[:-8] + "_v01.hf"
    hdf = hf.File(fn, "w")
    hdf.create_dataset(
        "datetime",
        data=(
            pd.Series(dfr.index) - pd.to_datetime("1970-01-01 00:00:00", utc=True)
        ).dt.total_seconds(),
    )
    for i in dfr.columns[2:]:
        hdf.create_dataset(i, data=np.array(dfr[i]))
    hdf.close()

    # Append the dataframe to the list
    ldf.append(dfr)

    count += 1
    print("Data saved for %s" % (fname[-50:]))
    print(count)

# Save all dataframes to one dataframe and save it in a new file
df_t = pd.concat(ldf)


fn = f"{Path(fnames[0]).name[:-8]}_{fnames[-1][-16:-8]}_v01.hf"

hdf = hf.File(fn, "w")
hdf.create_dataset(
    "datetime",
    data=(
        pd.Series(df_t.index) - pd.to_datetime("1970-01-01 00:00:00", utc=True)
    ).dt.total_seconds(),
)
for i in df_t.columns[2:]:
    hdf.create_dataset(i, data=np.array(df_t[i]))
hdf.close()

print(f"Data saved to file {fn}")

print(f"It took {np.round(time.time() - start, 2)} seconds to run the code")
