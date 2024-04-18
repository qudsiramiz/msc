# Load necessary module
import os.path
from glob import glob
from ftplib import FTP_TLS
import numpy as np
import warnings
from glob import glob

import matplotlib.pyplot as plt


def load_date(
    localpath,
    mission="mariner10",
    instrument="fields",
    instrument_1=None,
    level="l2",
    years=[2018],
    months=None,
    dates=None,
):
    """Define the function to download data from the NASA CDAWeb  for PSP
    mission. It takes following parameters as input:
    localpath:
    mission:
    instrument:
    instrument_1:
    level:
    years:
    months:
    dates:
    """

    for yy in years:
        for mm in months:
            for dd in dates:

                if mm < 10:
                    mm = str(mm).zfill(2)
                if dd < 10:
                    dd = str(dd).zfill(2)

                directory = "https://cdaweb.gsfc.nasa.gov/pub/data/psp/coho1hr_magplasma/cdf/"
                    fn = mission + "_fld_" + level + "_mag_rtn_4_sa_per_cyc_"

                fl0 = fn + str(yy) + str(mm) + str(dd) + "_v??.cdf"
                fl0_path = os.path.join(localpath, fl0)
                gb = glob(fl0_path)
                try:
                    ftp = FTP_TLS("spdf.gsfc.nasa.gov")
                    ftp.login()
                    ftp.cwd(directory)
                    ls = ftp.nlst(fl0)
                    fl = ls[-1]
                    savepath = os.path.join(localpath, fl)
                    if os.path.isfile(savepath):
                        print(
                            f"file already exists for {yy}/{mm}/{dd}, dowloading the next file in the list"
                        )
                    else:
                        ftp.retrbinary("RETR " + fl, open(savepath, "wb").write)
                        print(f"Data downloaded for {yy}/{mm}/{dd}")
                except:
                    print(
                        f"Data not downloaded for {yy}/{mm}/{dd}. Something went wrong, obviously!"
                    )
                    pass
    print(f"All files downloaded")


localpath = "/content/drive/Shared drives/Shadowfax/GIT/SpacePlasFest/qudsi/psp/data"
years = [2018]
months = np.arange(10, 12)
dates = np.arange(1, 5)
load_date(
    localpath=localpath,
    mission="psp",
    instrument="sweap",
    instrument_1="spc",
    level="l2",
    years=years,
    months=months,
    dates=dates,
)
