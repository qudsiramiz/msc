import os
import re

# Path to the root folder
root_folder = "/mnt/cephadrius/udel_research/msc/data/v2024.05/"

# Define patterns for old filenames and their replacements
# patterns = {
#     r"(.*)v2024\.5\.p$": r"\1v2024.05.p",
#     r"(.*)v2024\.5\.hf$": r"\1v2024.05.hf",
# }
# Replace all __ with _
patterns = {
    r"(.*)__(.*)": r"\1_\2",
}
# Walk through all files in the folder and subfolders
for dirpath, dirnames, filenames in os.walk(root_folder):
    for filename in filenames:
        for pattern, replacement in patterns.items():
            # If the filename matches the pattern, rename it
            if re.match(pattern, filename):
                old_file_path = os.path.join(dirpath, filename)
                new_filename = re.sub(pattern, replacement, filename)
                new_file_path = os.path.join(dirpath, new_filename)
                os.rename(old_file_path, new_file_path)
                print(f"Renamed: {old_file_path} -> {new_file_path}")
