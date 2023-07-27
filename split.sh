#!/bin/bash

set -x

# Load Modules 
module purge
module use $HOME/MyModules
module load miniconda3/latest
which python3
python --version 
which pip
conda list

# Set the input file paths
bacteria_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/BACTERIA/bacteria_selected.csv"
phage_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/PHAGE/phage_selected.csv"

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

# Set the output file paths for train, test, and dev splits
output_path="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT"
train_file="$output_path/train.csv"
test_file="$output_path/test.csv"
dev_file="$output_path/dev.csv"

# Set the train-test-dev split ratios
train_ratio=0.8
test_ratio=0.1
dev_ratio=0.1

# Run Python script to perform train-test-dev split
python_script=$(cat << END
import pandas as pd
from sklearn.model_selection import train_test_split

# Load the input CSV files
bacteria_data = pd.read_csv("$bacteria_file")
phage_data = pd.read_csv("$phage_file")

# Perform train-test-dev split
bacteria_train_data, bacteria_temp_data = train_test_split(bacteria_data, train_size=$train_ratio, test_size=$test_ratio+ $dev_ratio)
phage_train_data, phage_temp_data = train_test_split(phage_data, train_size=$train_ratio, test_size=$test_ratio+ $dev_ratio)


# Perform test-dev split on the remaining data
bacteria_test_data, bacteria_dev_data = train_test_split(bacteria_temp_data, train_size=$test_ratio/($test_ratio + $dev_ratio), test_size=$dev_ratio/($test_ratio + $dev_ratio))
phage_test_data, phage_dev_data = train_test_split(phage_temp_data, train_size=$test_ratio/($test_ratio + $dev_ratio), test_size=$dev_ratio/($test_ratio + $dev_ratio))

# Save train, test, and dev data to separate CSV files
bacteria_train_data.to_csv("$output_path/bacteria_train.csv", index=False)
bacteria_test_data.to_csv("$output_path/bacteria_test.csv", index=False)
bacteria_dev_data.to_csv("$output_path/bacteria_dev.csv", index=False)

phage_train_data.to_csv("$output_path/phage_train.csv", index=False)
phage_test_data.to_csv("$output_path/phage_test.csv", index=False)
phage_dev_data.to_csv("$output_path/phage_dev.csv", index=False)

# # Concat and shuffle separate train, test, and dev files
# train1 = pd.read_csv("$output_path/bacteria_train.csv")
# train2 = pd.read_csv("$output_path/phage_train.csv")
# train_data = pd.concat([train1, train2], ignore_index=True)
# train_data = train_data.sample(frac=1, random_state=42)

# test1 = pd.read_csv("$output_path/bacteria_test.csv")
# test2 = pd.read_csv("$output_path/phage_test.csv")
# test_data = pd.concat([test1, test2], ignore_index=True)
# test_data = test_data.sample(frac=1, random_state=42)

# dev1 = pd.read_csv("$output_path/bacteria_dev.csv")
# dev2 = pd.read_csv("$output_path/phage_dev.csv")
# dev_data = pd.concat([dev1, dev2], ignore_index=True)
# dev_data = dev_data.sample(frac=1, random_state=42)

# # Save train, test, and dev data to separate CSV files
# train_data.to_csv("$train_file", index=False)
# test_data.to_csv("$test_file", index=False)
# dev_data.to_csv("$dev_file", index=False)

END
)

# Execute the Python script
python3 -c "$python_script"

# Concat and shuffle train, test, and dev files
cd $output_path
cat bacteria_train.csv phage_train.csv | shuf > "$train_file"
cat bacteria_test.csv phage_test.csv | shuf > "$test_file"
cat bacteria_dev.csv phage_dev.csv | shuf > "$dev_file"

# Set the headers for train, test, and dev files
temp_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/temp.csv"

echo "sequence,label" | cat - "$train_file" > "$temp_file"
mv "$temp_file" "$train_file"

echo "sequence,label" | cat - "$test_file" > "$temp_file"
mv "$temp_file" "$test_file"

echo "sequence,label" | cat - "$dev_file" > "$temp_file"
mv "$temp_file" "$dev_file"

echo "TIME: End: = `date +"%Y-%m-%d %T"`" 