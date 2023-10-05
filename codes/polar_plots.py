import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import h5py as hf
from scipy import stats
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

from histogram_plots import brzl_hist

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


# Make a histogram of Tp.
x_grid, y_grid = np.mgrid[0:6:40j, 0:3:40j]
Tp_grid = stats.binned_statistic_2d(
    x, y, df["Tp"][:], statistic="median", bins=[x_grid[:, 0], y_grid[0, :]]
)[0]


x_grid, y_grid = np.mgrid[0:6:40j, 0:3:40j]
Tp_interp_cubic = griddata(
    np.array([x, y]).T,
    df["Tp"][:],
    (x_grid, y_grid),
    method="cubic",
)

Tp_interp_nearest = griddata(
    np.array([x, y]).T,
    df["Tp"][:],
    (x_grid, y_grid),
    method="nearest",
)

Tp_interp_linear = griddata(
    np.array([x, y]).T,
    df["Tp"][:],
    (x_grid, y_grid),
    method="linear",
)

cmap = "viridis_r"
# In a 2 by 2 subplot, plot the histogram of Tp and the histogram of the interpolated Tp.
fig, axs = plt.subplots(2, 2, figsize=(8, 4.5), sharex=True)
fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
im1 = axs[0, 0].pcolormesh(
    x_grid,
    y_grid,
    Tp_grid,
    cmap=cmap,
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)
im2 = axs[0, 1].pcolormesh(
    x_grid,
    y_grid,
    Tp_interp_cubic,
    cmap=cmap,
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)
im3 = axs[1, 0].pcolormesh(
    x_grid,
    y_grid,
    Tp_interp_nearest,
    cmap=cmap,
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)
im4 = axs[1, 1].pcolormesh(
    x_grid,
    y_grid,
    Tp_interp_linear,
    cmap=cmap,
    norm=mpl.colors.LogNorm(vmin=1e4, vmax=3e5),
)

divider1 = make_axes_locatable(axs[0, 0])
cax1 = divider1.append_axes("right", size="5%", pad=0.05)
cbar1 = plt.colorbar(im1, cax=cax1)
cbar1.ax.tick_params(labelsize=10)
# cbar1.set_label(r"$T_p$ (K)", fontsize=10)

divider2 = make_axes_locatable(axs[0, 1])
cax2 = divider2.append_axes("right", size="5%", pad=0.05)
cbar2 = plt.colorbar(im2, cax=cax2)
cbar2.ax.tick_params(labelsize=10)
cbar2.set_label(r"$T_p$ (K)", fontsize=10)

divider3 = make_axes_locatable(axs[1, 0])
cax3 = divider3.append_axes("right", size="5%", pad=0.05)
cbar3 = plt.colorbar(im3, cax=cax3)
cbar3.ax.tick_params(labelsize=10)
# cbar3.set_label(r"$T_p$ (K)", fontsize=10)

divider4 = make_axes_locatable(axs[1, 1])
cax4 = divider4.append_axes("right", size="5%", pad=0.05)
cbar4 = plt.colorbar(im4, cax=cax4)
cbar4.ax.tick_params(labelsize=10)
cbar4.set_label(r"$T_p$ (K)", fontsize=10)

# Set the axis labels.
axs[0, 0].set_aspect("equal")
axs[0, 0].set_xlabel(r"$X$ [AU]", fontsize=10)
axs[0, 0].set_ylabel(r"$Z$ [AU]", fontsize=10)
axs[0, 0].set_title("Original", fontsize=10)

axs[0, 1].set_aspect("equal")
axs[0, 1].set_xlabel(r"$X$ [AU]", fontsize=10)
# axs[0, 1].set_ylabel(r"$Z$ [AU]", fontsize=10)
axs[0, 1].set_title("Cubic", fontsize=10)

axs[1, 0].set_aspect("equal")
axs[1, 0].set_xlabel(r"$X$ [AU]", fontsize=10)
axs[1, 0].set_ylabel(r"$Z$ [AU]", fontsize=10)
axs[1, 0].set_title("Nearest", fontsize=10)

axs[1, 1].set_aspect("equal")
axs[1, 1].set_xlabel(r"$X$ [AU]", fontsize=10)
# axs[1, 1].set_ylabel(r"$Z$ [AU]", fontsize=10)
axs[1, 1].set_title("Linear", fontsize=10)

# Save the figure.
plt.savefig("../figures/uy_Tp_interp.pdf", dpi=300, bbox_inches="tight", pad_inches=0.1)
