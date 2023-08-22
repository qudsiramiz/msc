import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#df_omni_v03 = pd.read_pickle('/mnt/cephadrius/udel_research/msc/omni/data/processed/v03/omni_coho1hr_merged_mag_plasma_19630101_20211201_v03.p')
#
#df_omni_v02 = pd.read_pickle('/mnt/cephadrius/udel_research/msc/omni/data/processed/v02/omni_coho1hr_merged_mag_plasma_19630101_20211201_v02.p')
#
## Get rid of data where the values are less than 1.e-29
#df_omni_v02 = df_omni_v02[df_omni_v02.np > 1.e-29]
#df_omni_v03 = df_omni_v03[df_omni_v03.np > 1.e-29]

#plt.figure()
plt.clf()
plt.plot(df_omni.index, df_omni.np, 'rd', ms=1, label='v02')
plt.plot(df_omni_v03.index, df_omni_v03.np, 'k', ms=1, label='v03')

# Set the time limit of x-axis
plt.xlim(pd.Timestamp('2020-12-01'), pd.Timestamp('2021-01-01'))
plt.legend()
plt.show()