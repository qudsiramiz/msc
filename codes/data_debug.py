import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


df = pd.read_pickle("../data/all_data/all_spacecraft_data.p")
df_binned = pd.read_pickle("../data/all_data/v16/all_spacecraft_data_binned_76_v16.p")


#plt.close("all")

plt.figure()
key = "np"
plt.plot(df.sc_r, df[f"{key}"], 'r.', ms=2, alpha=1, label=f"{key}")
plt.plot(df_binned.sc_r_median, df_binned[f"{key}_median"], 'd', ms=3, alpha=1, label="binned_median")
#plt.plot(df_binned.sc_r_median, df_binned.Tp_median, 'd', ms=3, alpha=1, label="binned_median")
plt.xscale('log')
plt.yscale('log')

plt.xlabel("Distance from the Sun")
plt.ylabel("Temperature")
plt.legend(loc="best")
plt.xlim(0.06, 0.7)
plt.show()
