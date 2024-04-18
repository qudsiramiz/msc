from spacepy.pycdf import CDF as cdf
import numpy as np
from glob import glob
import pandas as pd
import datetime
import pytz
import h5py as hf
import time as tm
from pathlib import Path

start = tm.time()

mu_0 = 4 * np.pi * 1.0e-7
m_p = 1.67e-27
kb = 1.38e-23

#
data_dir = "/mnt/cephadrius/udel_research/msc/data/merged_1hr/v01/2024_version/"
save_dir = "/mnt/cephadrius/udel_research/msc/data/merged_1hr/v2024.1"

fnames = np.sort(glob(data_dir + "*.hf"))

fnames = list(fnames)

# Define a custom list of keys/parameters for which you want the binning to be
# done, since it doesn't make much sense to bin the data for parameters like
# time and elevation angle etc

key_list = [
    "br",
    "bt",
    "bn",
    "bm",
    "sc_r",
    "np",
    "Tp",
    "vA",
    "vp_m",
    "vp_r",
    "vp_t",
    "vp_n",
    "heliographicLatitude",
    "heliographicLongitude",
    "alfven_ratio",
    "parker_angle",
    "particle_flux",
    "zpr",
    "zpt",
    "zpn",
    "zpm",
    "zmr",
    "zmt",
    "zmn",
    "zmm",
    "sig_c",
    "proton_beta",
]

ldf = []
sc_n = []

for f in fnames[:]:
    print(f)
    d = {}
    dat = hf.File(f)

    # Get the spacecraft name from the file name
    sc_n = (Path(f).name).split("_")[0]

    for key in dat.keys():
        d[key] = dat[key]

    df = pd.DataFrame(d)
    df.index = pd.to_datetime(df.datetime, unit="s")

    try:
        df[df.np <= 0] = np.nan
    except Exception:
        pass

    if "vA" in df.keys():
        print("Alfven velocity already defined")
    else:
        df["vA"] = np.empty_like(df[df.columns[0]])
        df["vA"][:] = np.nan
        try:
            df["vA"] = 1e-9 * df.bm / np.sqrt(mu_0 * 1.0e6 * df.np * m_p)
        except Exception:
            pass

    try:
        df["alfven_ratio"] = 1.0e3 * df.vp_m / df.vA
    except Exception:
        pass

    try:
        df["proton_beta"] = 1.0e6 * df.np * kb * df.Tp * 2 * mu_0 / (1e-9 * df.bm) ** 2
    except Exception:
        pass

    try:
        # Adding Parker Angle to the dataframe (units: radians)
        bbm = df.br.values**2 + df.bt.values**2 + df.bn.values**2
        sbm = np.array([np.sqrt(xx) for xx in bbm])
        df["parker_angle"] = np.arccos(abs(df.br) / sbm)
    except Exception:
        pass

    try:
        # Adding particle flux to the dataframe (units: 1/second)
        df["particle_flux"] = (
            (df.np * 1e6) * (df.vp_m * 1e3) * ((df.sc_r * 1.49e11) ** 2)
        )
    except Exception:
        pass

    # Set the value of rest of the keys to np.nan
    for key in key_list:
        if key in df.keys():
            print(key + " already present in the file")
        else:
            df[key] = np.empty_like(df[df.columns[0]])
            df[key] = np.nan
            print(key + " value set to np.nan for file %s" % (f))

    print("Saving data to a pickle file for all keys")

    # Assign a new name to the pickle file
    # For omni: [56:-7]
    # For ulysses: [59:-7]
    # For all other spacecrafts: [51:-7]
    nf = save_dir + "/" + Path(f).name[:-7] + "_v2024.1.p"
    nf_hf = save_dir + "/" + Path(f).name[:-7] + "_v2024.1.hf"

    # Save the binned data to a pickel file using protocol 2. Avoid other
    # protocols, since protocol greater than 2 doesn't work for Python 2.*
    pd.DataFrame.to_pickle(df, nf, protocol=2)

    hdf = hf.File(nf_hf, "w")
    hdf.create_dataset(
        "datetime",
        data=(
            pd.Series(df.index) - pd.to_datetime("1970-01-01 00:00:00")
        ).dt.total_seconds(),
    )
    for i in df.columns[:]:
        if i != "datetime":
            hdf.create_dataset(i, data=np.array(df[i]))
    hdf.close()

    print("Data saved to pickle file for %s" % (nf))

print("It took %s seconds to run this code" % (round(tm.time() - start, 2)))
