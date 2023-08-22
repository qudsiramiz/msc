import pickle

# Load the pickle file
with open("../omni/data/processed/v03/omni_solar_cycle_average_values_v02.p", "rb") as f:
    omni_av_data = pickle.load(f)

# Get the average values for the 24th solar cycle
avg_val_24 = omni_av_data[0]

# Print the average values for the 24th solar cycle
avg_type = "mean"  # Different types: mean, median, 10%, 25%, 50%, 75%, 90%
print(f"Printing {avg_type} values")
for key in avg_val_24.keys():
    print(f"{key}: {avg_val_24[key][avg_type]}")
