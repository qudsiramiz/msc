import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py as hf
from scipy import stats
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

# Imprt gaussian processed from sklearn.
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

# Activate the latex text rendering.
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

data = hf.File("../data/uy_coho1hr_merged_mag_plasma_19900101_20091201.hf", "r")

# Convert data to a pandas dataframe.
df = pd.DataFrame()
for key in data.keys():
    df[key] = data[key][:]

# Get rid of all the rows where Tp is NaN.
df = df[np.isfinite(df["Tp"][:])]

# Get the x and y data.
x = abs(df["sc_r"][:] * np.cos(np.deg2rad(df["heliographicLatitude"][:])))
y = abs(df["sc_r"][:] * np.sin(np.deg2rad(df["heliographicLatitude"][:])))

x = x.values
y = y.values

# Get a meshgrid of x and y in logspace.
x_min = 0.01
x_max = np.log10(x.max())
y_min = 0.01
y_max = np.log10(y.max())
x_grid, y_grid = np.mgrid[x_min:x_max:40j, y_min:y_max:40j]
x_grid = 10**x_grid
y_grid = 10**y_grid
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

# Create a 2d interpolation of the data using scipy.
x_grid, y_grid = np.mgrid[0:6:40j, 0:3:40j]
Tp_interp = griddata(
    np.array([x, y]).T,
    df["Tp"][:],
    (x_grid, y_grid),
    method="cubic",
)

# Plot the data.
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.pcolormesh(x_grid, y_grid, Tp_grid, cmap="RdBu_r", norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5))

# Set x and y scales to log.
ax.set_xscale("log")
ax.set_yscale("log")

# Make a histogram of Tp.
# x_grid, y_grid = np.mgrid[0:6:10j, 0:3:10j]
# Tp_grid = stats.binned_statistic_2d(
#     x, y, df["Tp"][:], statistic="median", bins=[x_grid[:, 0], y_grid[0, :]]
# )[0]

# Using the Gaussian Processes interpolate the data to every point in the grid.
# Define the kernel.
# kernel = C(1.0, (1e-3, 1e3)) * RBF(10, (1e-2, 1e2))
# Define the Gaussian Process.
# gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=9)
# Fit the data.
# gp.fit(np.array([x, y]).T, df["Tp"][:])
# Predict the data.
#  Tp_grid = gp.predict(np.array([x_grid.flatten(), y_grid.flatten()]).T).reshape(
#      x_grid.shape
# )

# Plot the data.
# fig, ax = plt.subplots(figsize=(8, 8))
# im = ax.pcolormesh(x_grid, y_grid, Tp_grid, cmap="RdBu_r", norm=mpl.colors.LogNorm())

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(r"$T_p$ (K)")

# Set the axis labels.
ax.set_aspect("equal")
ax.set_xlabel(r"$x$ (AU)")
ax.set_ylabel(r"$y$ (AU)")
ax.set_title(r"$T_p$ (K)")

#
plt.savefig("../figures/uy_Tp_log.png", dpi=300, bbox_inches="tight", pad_inches=0.1)

# Plot the interpolated data.
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.pcolormesh(x_grid, y_grid, Tp_interp, cmap="RdBu_r", norm=mpl.colors.LogNorm())

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, cax=cax)
cbar.ax.tick_params(labelsize=10)
cbar.ax.set_ylabel(r"$T_p$ (K)")

# Set the axis labels.
ax.set_aspect("equal")
ax.set_xlabel(r"$x$ (AU)")
ax.set_ylabel(r"$y$ (AU)")
ax.set_title(r"$T_p$ (K)")

# Save the figure.
plt.savefig("../figures/uy_Tp_intrp.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
# Plot the interpolated data.
# fig, ax = plt.subplots(figsize=(8, 8))
# im = ax.pcolormesh(x_grid, y_grid, Tp_grid, cmap="RdBu_r", norm=mpl.colors.LogNorm())
#
# divider = make_axes_locatable(ax)
# cax = divider.append_axes("right", size="5%", pad=0.05)
# cbar = fig.colorbar(im, cax=cax)
# cbar.ax.tick_params(labelsize=10)
# cbar.ax.set_ylabel(r"$T_p$ (K)")
#
# # Set the axis labels.
# ax.set_aspect("equal")
# ax.set_xlabel(r"$x$ (AU)")
# ax.set_ylabel(r"$y$ (AU)")
# ax.set_title(r"$T_p$ (K)")
#
# # Save the figure.
# plt.savefig("../figures/uy_Tp_gp.png", dpi=300, bbox_inches="tight", pad_inches=0.1)
# # Plot the heliographic latitude with respect to the radial distance.
# fig, ax = plt.subplots(figsize=(8, 8))
# ax.plot(df["sc_r"][:], abs(df["heliographicLatitude"][:]), "k.")
# ax.set_xlabel(r"$r$ (AU)")
# ax.set_ylabel(r"$\theta$ (deg)")
# ax.set_title(r"$\theta$ vs $r$")
# plt.savefig(
#     "../figures/uy_theta_vs_r.png", dpi=300, bbox_inches="tight", pad_inches=0.1
# )
