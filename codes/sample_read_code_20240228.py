import pandas as pd
import matplotlib.pyplot as plt

"""
List of keys in the dataframe:
'datetime': Date and time of the observation.
'br': Radial component of the magnetic field in nanoTesla.
'bt': Tangential component of the magnetic field in nanoTesla.
'bn': Normal component of the magnetic field in nanoTesla.
'bm': Magnetic field magnitude in nanoTesla.
'vp_r': Proton radial velocity in kilometers per second.
'vp_t': Proton tangential velocity in kilometers per second.
'vp_n': Proton normal velocity in kilometers per second.
'vp_m': Protom velocity magnitude in kilometers per second.
'heliographicLatitude': HelioGraphic Latitude in degrees.
'heliographicLongitude': Heliographic Longitude in degrees.
'sc_r': Scapecraft radial distance from the sun in astronomical units.
'vA': Alfven velocity in kilometers per second.
'np': Proton number density in numbers per cubic centimeter.
'Tp': Proton temperature in Kelvin.
"""

fn = "../data/all_data/v19/all_spcaecraft_data_v19.p"
fns = "../data/all_data/v19/all_spcaecraft_data_scaled_v19.p"

df = pd.read_pickle(fn)
dfs = pd.read_pickle(fns)

# Get the data corresponding to the list of keys.
df_new = df[
    [
        "datetime",
        "br",
        "bt",
        "bn",
        "bm",
        "vp_r",
        "vp_t",
        "vp_n",
        "vp_m",
        "heliographicLatitude",
        "heliographicLongitude",
        "sc_r",
        "vA",
        "np",
        "Tp",
    ]
]
dfs_new = dfs[
    [
        "br",
        "bt",
        "bn",
        "bm",
        "vp_r",
        "vp_t",
        "vp_n",
        "vp_m",
        "sc_r",
        "vA",
        "np",
        "Tp",
    ]
]


# Save the new dataframe to a pickle file.
fn_new = fn.replace(".p", "_selected_scalar.p")
fns_new = fns.replace(".p", "_selected_scalar.p")
df_new.to_pickle(fn_new)
dfs_new.to_pickle(fns_new)
print(f"Dataframe saved to {fn_new}")
