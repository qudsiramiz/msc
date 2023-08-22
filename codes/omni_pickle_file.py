import numpy as np
import pandas as pd
import datetime
import pytz

fn = '/media/cephadrius/endless/udel_research/msc/omni/data/processed/v02/test.p'

# Read the pickle file
df = pd.read_pickle(fn)

# set the index to the time column
df.set_index(pd.to_datetime(df['datetime']), inplace=True)
