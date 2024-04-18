import h5py as hf
import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import os

from glob import glob

fn = "/mnt/cephadrius/udel_research/msc/voyager/data/processed/v00/voyager1_coho1hr_merged_mag_plasma_19770301.hf"

# dat = hf.File(fn, "r+")

# Print the keys of the file
# print(dat.keys())
# dat_diff = dat["datetime"][1:] - dat["datetime"][:-1]

# Plot the difference in time with respect to the index

# plt.plot(dat_diff, ".")
# plt.xlabel("Index")
# plt.ylabel("Time difference")
# plt.title("Time difference between consecutive data points")
# plt.savefig("time_diff.png")


def download_files(
    data_link, file_name_style, start_year, end_year, start_month, end_month, local_path
):

    for year in range(start_year, end_year + 1):
        for month in range(start_month, end_month + 1):

            year = str(year).zfill(4)
            month = str(month).zfill(2)

            # Define the file name
            file_name = f"{file_name_style}_{year}{month}01_v01.cdf"

            final_file_name = f"{data_link}/{year}/{file_name}"
            # Check if file_name exists in the local_path directory
            # Download the file
            # Check if the file already exists in the local path
            if os.path.isfile(local_path + "/" + file_name):
                print(
                    f"File already exists for {year}/{month}, dowloading the next file in the list"
                )
            else:
                try:
                    os.system(f"wget -P {local_path} {final_file_name}")
                    print("File downloaded to " + local_path + file_name)
                except Exception as e:
                    print(e)
            # os.system("wget -P " + local_path + " " + data_link + file_name)


## For PSP
# Link to downaload the data from
data_link = "https://cdaweb.gsfc.nasa.gov/pub/data/psp/coho1hr_magplasma/cdf"
file_name_style = "psp_coho1hr_merged_mag_plasma"
start_year = 2022
end_year = 2024
start_month = 1
end_month = 12
local_path = "/mnt/cephadrius/udel_research/msc/psp/coho1hr"
# local_path = "/home/cephadrius/Desktop/git/msc/codes/v2024.1/"

# download_files(
#     data_link, file_name_style, start_year, end_year, start_month, end_month, local_path
# )

## For Solar Orbiter
# Link to downaload the data from
data_link = "https://cdaweb.gsfc.nasa.gov/pub/data/solar-orbiter/coho1hr_magplasma"
file_name_style = "solo_coho1hr_merged_mag_plasma"
start_year = 2020
end_year = 2024
start_month = 1
end_month = 12
local_path = "/mnt/cephadrius/udel_research/msc/solo"

# download_files(
#     data_link, file_name_style, start_year, end_year, start_month, end_month, local_path
# )

## For OMNI
# Link to downaload the data from
data_link = "https://cdaweb.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/coho1hr_magplasma/"
file_name_style = "omni_coho1hr_merged_mag_plasma"
start_year = 2021
end_year = 2024
start_month = 1
end_month = 12
local_path = "/mnt/cephadrius/udel_research/msc/omni/coho1hr_magplasma"

download_files(
    data_link, file_name_style, start_year, end_year, start_month, end_month, local_path
)
