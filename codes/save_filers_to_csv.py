"""
The command should be something like the following.
 
<dataframe>.to_csv(<file name>,
sep=”,”,
na_rep=”NaN”,
float_format=”%.8e”, # 8 decimal place float in exponential notation,
header=True,
index=True,
compression=None)
 
 
You’ll want to manually inspect the output before you finalize it, so I recommend something like the following to select the first 5 entries to keep your screen clean.
 
<dataframe>.iloc[:5].to_csv(<file name>,
sep=”,”,
na_rep=”NaN”,
float_format=”%.8e”,
header=True,
index=True,
compression=None)
"""
import glob as glob
import numpy as np
import pandas as pd

data_dir_read = "/mnt/cephadrius/udel_research/msc/data/all_data/v19/"
data_dir_save = "/mnt/cephadrius/udel_research/msc/data/all_data/v19/csv_files/"

fnames = np.sort(glob.glob(data_dir_read + "*.p"))
key_type_list = ['_iqr_10', '_iqr_25', '_iqr_50', '_iqr_75', '_iqr_90', '_num']
key_list = ['bm', 'np', 'vp_m', 'Tp', 'vA', 'alfven_ratio', 'proton_beta']

ind1 = 1
ind2 = ind1 + 1

for fname in fnames[ind1:ind2]:
    df = pd.read_pickle(fname)
    for key in key_list:
        dfn = pd.DataFrame()
        for key_type in key_type_list:
            dfn[key + key_type] = df[key + key_type]
            # if key + '_num' is zero, then set it to NaN
        dfn.loc[dfn[key + '_num'] == 0, key + key_type] = np.nan

        # Add sc_r_iqr_50 to dfn in scientific notation and set it to index
        dfn['sc_r_iqr_50'] = df['sc_r_iqr_50']
        dfn.index = dfn['sc_r_iqr_50'].apply(lambda x: '%.2e' % x)
        dfn = dfn.drop(columns=['sc_r_iqr_50'])
        # Remove rows which have all columns as NaN
        dfn = dfn.dropna(how='all')
        # Save dfn to csv file with name of the key in the file name
        dfn.to_csv(data_dir_save + key + "_" + fname.split("/")[-1].split(".")[0] + ".csv",
                   sep=",", na_rep="", float_format="%.2e", header=True, index=True,
                   compression=None)

    # dfn.to_csv(data_dir_save + fname.split("/")[-1].split(".")[0] + ".csv", sep=",",
    #            na_rep="NaN", float_format="%.3e", header=True, index=True, compression=None)
    
    print(fname.split("/")[-1].split(".")[0] + ".csv")
    # print(dfn.iloc[:5].to_csv(sep=",", na_rep="NaN", float_format="%.2e", header=True, index=True,
    #       compression=None))