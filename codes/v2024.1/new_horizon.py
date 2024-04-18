import pandas as pd
import spacepy.pycdf as cdf
import numpy as np

new_horizon_folder = "/mnt/cephadrius/udel_research/msc/new-horizons/data/validsum/"

file_0 = "new_horizons_swap_validsum_20081010210700_v1.0.3.cdf"
file_1 = "new_horizons_swap_validsum_20081010210700_v1.0.5.cdf"
file_2 = "new_horizons_swap_validsum_20081010210700_v1.0.7.cdf"

# Read the files
dat0 = cdf.CDF(new_horizon_folder + file_0)
dat1 = cdf.CDF(new_horizon_folder + file_1)
dat2 = cdf.CDF(new_horizon_folder + file_2)

# Print the keys of the file
print(dat0.keys())
print(dat1.keys())
print(dat2.keys())
