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


def gp_ext(x, y, z, x_grid, y_grid, kernel, alpha=0.0):
    """
    Extend the data using a gaussian process.
    """
    # Initialize the gaussian process.
    gp = GaussianProcessRegressor(kernel=kernel, alpha=alpha)

    # Fit the gaussian process.
    gp.fit(np.array([x, y]).T, z)

    # Predict the values of the grid.
    z_grid, sigma = gp.predict(
        np.array([x_grid.flatten(), y_grid.flatten()]).T, return_std=True
    )

    # Reshape the grid.
    z_grid = z_grid.reshape(x_grid.shape)
    sigma = sigma.reshape(x_grid.shape)

    return z_grid, sigma


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
x = abs(df["sc_r"][:] * np.cos(np.deg2rad(df["heliographicLatitude"][:])))
y = abs(df["sc_r"][:] * np.sin(np.deg2rad(df["heliographicLatitude"][:])))

x = x.values
y = y.values

X = np.array([x, y]).T
z = df["Tp"][:]
z = z.values

# Get a meshgrid of x and y in linear space.
x_min = x.min() * 0.9
x_max = x.max() * 1.1
y_min = y.min() * 0.9
y_max = y.max() * 1.1
# Get the grid in logspace.
nbins = 40
x_lin = np.linspace(x_min, x_max, nbins)
y_lin = np.linspace(y_min, y_max, nbins)
x_grid, y_grid = np.meshgrid(x_lin, y_lin)

Tp_grid = np.full(x_grid.shape, np.nan)  # Initialize the grid.

# Compute the median Tp for each bin.
for i in range(len(x_grid) - 1):
    for j in range(len(y_grid) - 1):
        tk = np.where(
            (x >= x_grid[i, j])
            & (x < x_grid[i + 1, j])
            & (y >= y_grid[i, j])
            & (y < y_grid[i, j + 1])
        )[0]
        if len(tk) > 0:
            Tp_grid[i, j] = np.nanmedian(df["Tp"].values[tk])
        else:
            Tp_grid[i, j] = np.nan

X = np.column_stack((x_grid.reshape(-1), y_grid.reshape(-1)))
Z = Tp_grid.reshape(-1, 1)
guess_l = (2.0, 2)
bounds_l = ((1e-3, 1e3), (1e-3, 1e3))
guess_n = 1
bounds_n = (1e-3, 1e3)
kernel = ConstantKernel(guess_n, bounds_n) * RBF(guess_l, bounds_l) + Matern(
    guess_l, bounds_l
)

X_train, X_test, Y_train, Y_test = train_test_split(
    X, Z, test_size=0.95, random_state=0
)
gpr = GaussianProcessRegressor(
    kernel=kernel, alpha=0.01, n_restarts_optimizer=10, normalize_y=True
)
gpr.fit(X_train, Y_train)
Zfit, Zsigma = gpr.predict(X, return_std=True)
Zfit = Zfit.reshape(x_grid.shape)
Zsigma = Zsigma.reshape(x_grid.shape)

# Set up the figure
fig, axs = plt.subplots(figsize=(4, 6), sharex=True)
axs.set_xlabel(r"$x$ (AU)", fontsize=10)
axs.set_ylabel(r"$y$ (AU)", fontsize=10)
axs.set_xlim([x_min, x_max])
axs.set_ylim([y_min, y_max])
axs.set_xscale("log")
axs.set_yscale("log")

# Do the plotting
lev = np.logspace(0, 250, 10)
axs.contour(x, y, lev, colors="k", linewidths=0.5)
axs.plot(*X_train.T, "k.", markersize=1)
axs.contour(x, y, Zfit, lev, colors="r", linewidths=0.5, ls="--")
plt.tight_layout()
plt.savefig("../figures/all_Tp_gp.png", dpi=300, bbox_inches="tight", pad_inches=0.1)

"""
# Flatten x_grid, y_grid, and Tp_grid.
x_flat = x_grid.flatten()
y_flat = y_grid.flatten()
Tp_flat = Tp_grid.flatten()

# Get rid of all the NaNs.
x_flat = x_flat[~np.isnan(Tp_flat)]
y_flat = y_flat[~np.isnan(Tp_flat)]
Tp_flat = Tp_flat[~np.isnan(Tp_flat)]

# Define a combination of kernels.
kernel = ConstantKernel(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-3, 1e3))
# Carry out the gaussian process.
z_grid, sigma = gp_ext(
    x_flat,
    y_flat,
    Tp_flat,
    x_grid,
    y_grid,
    kernel=kernel,
    alpha=0.01,
)

# Plot the z-grid and the sigma on a 2 by 1 subplot.
fig, axs = plt.subplots(2, 1, figsize=(4, 6), sharex=True)
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
im1 = axs[0].pcolormesh(
    x_grid,
    y_grid,
    z_grid,
    cmap="viridis_r",
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)
im2 = axs[1].pcolormesh(
    x_grid,
    y_grid,
    sigma,
    cmap="viridis_r",
    norm=mpl.colors.Normalize(vmin=0, vmax=1),
)

# Add a colorbar to the plot.
divider1 = make_axes_locatable(axs[0])
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(im1, cax=cax1)
divider2 = make_axes_locatable(axs[1])
cax2 = divider2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(im2, cax=cax2)

# Set the x and y scales to log.
axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[1].set_xscale("log")
axs[1].set_yscale("log")

# Set the x and y limits.
axs[0].set_xlim([x_min, x_max])
axs[0].set_ylim([y_min, y_max])
axs[1].set_xlim([x_min, x_max])
axs[1].set_ylim([y_min, y_max])

# Set the x and y labels.
axs[0].set_ylabel(r"$y$ (AU)", fontsize=10)
axs[1].set_xlabel(r"$x$ (AU)", fontsize=10)

# Set the titles.
axs[0].set_title(r"$T_p$ (K)", fontsize=10)
axs[1].set_title(r"$\sigma$", fontsize=10)

# Save the figure.
plt.savefig("../figures/uy_Tp_gp.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
"""
