import numpy as np
import matplotlib.pyplot as plt

"""
# Constants for r_bin1
r_min1 = -1.2
r_max1 = 2
n_bin1 = 80
step_bin1 = (10**r_max1 / 10**r_min1) ** (1 / (n_bin1 - 1))

# Constants for r_bin2
r_min2 = -1.4
r_max2 = 2

# Determine the best `x`
min_diff = float("inf")
best_x = None

# Iterate over possible `x` values to find the one closest to `step_bin1`
for x in range(40, 121):  # iterate from 40 to 120 bins
    step_bin2 = (10**r_max2 / 10**r_min2) ** (1 / (x - 1))
    diff = abs(step_bin2 - step_bin1)

    if diff < min_diff:
        min_diff = diff
        best_x = x
"""
print("Best x value:", best_x)

r_min = -1.2
r_max = 2
n_bin1 = 80
r_bin1 = np.logspace(r_min, r_max, n_bin1)

# Append 5 values to the start of r_bin1 as NaN
r_bin1 = np.append(np.full(5, np.nan), r_bin1)

r_min = -1.4
r_max = 2
n_bin2 = 85
r_bin2 = np.logspace(r_min, r_max, n_bin2)

# Make a comparison plot
plt.figure(figsize=(6, 6))
plt.plot(r_bin2 - r_bin1, label="r_bin2 - r_bin1", marker="o")
# plt.plot(r_bin2, label="r_bin2")
# plt.plot(r_bin2, np.ones_like(r_bin2), label="r_bin2")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("r_bin")
plt.legend()
plt.show()
