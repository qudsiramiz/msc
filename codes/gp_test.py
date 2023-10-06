import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py as hf
import matplotlib as mpl

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern
from sklearn.model_selection import train_test_split

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats

# Activate the latex text rendering.
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

"""
# data = hf.File("../data/uy_coho1hr_merged_mag_plasma_19900101_20091201.hf", "r")
df = pd.read_pickle("../data/all_data/v19/all_spcaecraft_data_v19.p")

# Select data in the range where heliographic latitude is between -30 and 30 and sc_r is between 0.1
# and 30.
lat_lim = 5
r_lim = 10
df = df[
    (df["heliographicLatitude"] >= -lat_lim) & (df["heliographicLatitude"] <= lat_lim)
]
df = df[(df["sc_r"] >= 0.1) & (df["sc_r"] <= r_lim)]

# Convert data to a pandas dataframe.
# df = pd.DataFrame()
# for key in data.keys():
#     df[key] = data[key][:]

# Get rid of all the rows where Tp is NaN.
df = df[np.isfinite(df["Tp"][:])]

# Get the x and y data.
x_i = abs(df["sc_r"][:] * np.cos(np.deg2rad(df["heliographicLatitude"][:])))
y_i = abs(df["sc_r"][:] * np.sin(np.deg2rad(df["heliographicLatitude"][:])))

# Exclude all datapoints in the dataframe where x_val is less than 0.01 and y_val is less than 0.01.
df = df[(x_i >= 0.01) & (y_i >= 0.01)]

x_val = abs(df["sc_r"][:] * np.cos(np.deg2rad(df["heliographicLatitude"][:])))
y_val = abs(df["sc_r"][:] * np.sin(np.deg2rad(df["heliographicLatitude"][:])))

x_val = np.log10(x_val.values)
y_val = np.log10(y_val.values)

x_min = x_val.min() * 0.9
x_max = x_val.max() * 1.1
y_min = y_val.min() * 0.9
y_max = y_val.max() * 1.1

nv = 20
xv = np.linspace(x_min, x_max, nv)
yv = np.linspace(y_min, y_max, nv)
x_grid, y_grid = np.meshgrid(xv, yv)

# Get the Tp data.
Tp_val = df["Tp"][:].values

# Get the median value of Tp for each grid cell.
Tp_grid = np.zeros((nv, nv))
Tp_grid[:] = np.nan
for i in range(nv - 1):
    for j in range(nv - 1):
        # Get the indices of the data that lie within the grid cell.
        ind = (
            (x_val >= xv[i])
            & (x_val < xv[i + 1])
            & (y_val >= yv[j])
            & (y_val < yv[j + 1])
        )
        # If there is at least one data point in the grid cell, get the median value of Tp.
        if np.sum(ind) > 0:
            Tp_grid[i, j] = np.median(Tp_val[ind])

# Get the x_clean and y_clean only where Tp_grid is not NaN.
x = x_grid[~np.isnan(Tp_grid)]
y = y_grid[~np.isnan(Tp_grid)]
z = Tp_grid[~np.isnan(Tp_grid)]

X = np.array([x, y]).T
Z = z

# Split the data into training and testing sets.
guess_l = (1.0, 1.0)
bounds_l = ((1e-3, 1e0), (1e-3, 1e0))
guess_n = 1
bounds_n = (1e-10, 1e1)
kernel = ConstantKernel(guess_n, bounds_n) * Matern(guess_l, bounds_l)

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Z, test_size=0.02, random_state=0
)

# Fit the GP model.
gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100)
gp.fit(X_train, Y_train)

# Get the predicted values.
Zfit, Zsigma = gp.predict(X, return_std=True)

Xall = np.column_stack((x_grid.reshape(-1), y_grid.reshape(-1)))
Zall, ZsigmaA_all = gp.predict(Xall, return_std=True)
Zall = Zall.reshape(x_grid.shape)
Zsigma_all = ZsigmaA_all.reshape(x_grid.shape)

# Ignnore values lower than 1e3.
Zall[Zall < 1e3] = np.nan
"""
# Plot the results.
fig, axs = plt.subplots(1, 2, figsize=(8, 4.5), sharex=True)
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

# Plot the Zall
im1 = axs[0].pcolormesh(
    10 ** (x_grid),
    10 ** (y_grid),
    Zall.T,
    cmap="viridis_r",
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)
axs[0].set_xlabel(r"$x$ (AU)", fontsize=10)
axs[0].set_ylabel(r"$y$ (AU)", fontsize=10)
axs[0].set_title(r"$T_p$ (K)", fontsize=10)
# Set both axes to log scale.
axs[0].set_xscale("log")
axs[0].set_yscale("log")
divider1 = make_axes_locatable(axs[0])
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(im1, cax=cax1)
cbar1.ax.tick_params(labelsize=10)

# Plot the Zsigma_all
im2 = axs[1].pcolormesh(
    10 ** (x_grid),
    10 ** (y_grid),
    Tp_grid.T,
    cmap="viridis_r",
    norm=mpl.colors.LogNorm(),
)
axs[1].set_xlabel(r"$x$ (AU)", fontsize=10)
axs[1].set_ylabel(r"$y$ (AU)", fontsize=10)
axs[1].set_title(r"$\sigma$", fontsize=10)
# Set both axes to log scale.
axs[1].set_xscale("log")
axs[1].set_yscale("log")
divider2 = make_axes_locatable(axs[1])
cax2 = divider2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(im2, cax=cax2)
cbar2.ax.tick_params(labelsize=10)

# Save the figure.
plt.savefig("../figures/gp_ext_all.png", dpi=300, bbox_inches="tight")
