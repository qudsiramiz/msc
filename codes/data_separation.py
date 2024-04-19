import pandas as pd

fn = "../data/all_data/v19/all_spacecraft_data_binned_scaled_80_binned_v19.p"

df = pd.read_pickle(fn)

key_list = df.keys()

desired_keys = [
    "sc_r_median",
    "Tp_median",
    "bm_median",
    "np_median",
    "vp_m_median",
    "sc_r_iqr_10",
    "Tp_iqr_10",
    "bm_iqr_10",
    "np_iqr_10",
    "vp_m_iqr_10",
    "sc_r_iqr_90",
    "Tp_iqr_90",
    "bm_iqr_90",
    "np_iqr_90",
    "vp_m_iqr_90",
]

# Create a new dataframe with only the desired keys.
df_new = pd.DataFrame()
for key in desired_keys:
    df_new[key] = df[key][:]
# Save the new dataframe to a pickle file.
fn_new = fn.replace(".p", "_selected_scalar.p")
df_new.to_pickle(fn_new)
print(f"Dataframe saved to {fn_new}")
