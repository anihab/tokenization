#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=soc-gpu-np
#SBATCH --partition=soc-gpu-np
#SBATCH --job-name=tokenization
#SBATCH --time=1:00:00
#SBATCH -o /uufs/chpc.utah.edu/common/home/u1323098/Split%j.outerror

set -x

# module load anaconda3
# echo "starting dna env on conda"
# conda activate dna
module purge
module use $HOME/MyModules
module load miniconda3/latest
which python3
python --version
which pip
conda list

# Set the input file path
input_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/formatted.csv"

echo "TIME: Start: = `date +"%Y-%m-%d %T"`"

# Set the output file paths for train, test, and dev splits
train_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/train.csv"
test_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/test.csv"
dev_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/dev.csv"

# Set the train-test-dev split ratios
train_ratio=0.8
test_ratio=0.1
dev_ratio=0.1

# Run Python script to perform train-test-dev split
python_script=$(cat << END
import pandas as pd
from sklearn.model_selection import train_test_split

# Load the input CSV file
data = pd.read_csv("$input_file")

# Perform train-test-dev split
train_data, temp_data = train_test_split(data, train_size=$train_ratio, test_size=$test_ratio+ $dev_ratio)

# Perform test-dev split on the remaining data
test_data, dev_data = train_test_split(temp_data, train_size=$test_ratio/($test_ratio + $dev_ratio), test_size=$dev_ratio/($test_ratio + $dev_ratio))

# Save train, test, and dev data to separate CSV files
train_data.to_csv("$train_file", index=False)
test_data.to_csv("$test_file", index=False)
dev_data.to_csv("$dev_file", index=False)

END
)

# Execute the Python script
python3 -c "$python_script"

# Set the headers for train, test, and dev files
temp_file="/uufs/chpc.utah.edu/common/home/sundar-group2/ANISA/RAW_DATA/SPLIT/temp.csv"

echo "sequence,label" | cat - "$train_file" > "$temp_file"
mv "$temp_file" "$train_file"

echo "sequence,label" | cat - "$test_file" > "$temp_file"
mv "$temp_file" "$test_file"

echo "sequence,label" | cat - "$dev_file" > "$temp_file"
mv "$temp_file" "$dev_file"

echo "TIME: End: = `date +"%Y-%m-%d %T"`" 