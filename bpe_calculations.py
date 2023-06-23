import sys
import csv

## For analyzing bpe tokenization results

def get_tokens(file_path):
    tokens = []

    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row:
                tokens.append(row[0])

    return tokens

def calculate_average_word_length(values):
    total_length = 0
    word_count = 0

    for value in values:
        words = value.split()
        for word in words:
            total_length += len(word)
            word_count += 1

    if word_count == 0:
        return 0

    average_length = total_length / word_count
    return average_length

csv.field_size_limit(sys.maxsize)

# Example usage
file_path = '/uufs/chpc.utah.edu/common/home/u1323098/anisa/TOKENIZED_DATA/BPE/KR011061_tokenized.csv'  # Replace with the path to your CSV file

# Get values from the first column
first_column_values = get_tokens(file_path)

# Calculate average word length
average_length = calculate_average_word_length(first_column_values)
print(f"The average length of all tokens is: {average_length}")
