# -*- coding: utf-8 -*-

# Libraries
import argparse
import gzip
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from tokenizers import Tokenizer, models, trainers, normalizers
from transformers import PreTrainedTokenizerFast

# Globals
MAX_TOKENS = 510
VOCAB_SIZE = 5000
MAX_VOCAB_TOKEN_LENGTH = 2000

## Tokenize sequences given directory input

"""\
Given a bacteria directory and a phage directory, tokenize all files
according to method of choice.

Arguments:
  bacteria_dir -- str, path to directory of bacteria fasta files
  phage_dir -- str, path to directory of phage fasta files
  method -- str, tokenization method of choice
  b_out -- str, path to directory for all or bacteria output
  p_out -- str, path to directory for phage ouput
  k -- int, length of k if using kmer tokenization
  vocab -- str, path to directory of fasta files to build model vocab
"""
def read_files(bacteria_dir, phage_dir, method, **kwargs):
  b_out = kwargs.get('b_out', None)
  p_out = kwargs.get('p_out', b_out)
  k = kwargs.get('k', None)
  vocab = kwargs.get('vocab', None)

  # Ensure that input is correct
  if method == 'kmer' and k == None:
    print("Missing argument, kmer tokenization method requires parameter \'k\'.")
  if method == 'bpe' and vocab == None:
    print("Missing argument, bpe tokenization method requires parameter \'vocab\'.")

  # Build model vocabulary if using bpe
  if method == 'bpe':
    build_vocab(vocab)

  # Tokenize all files in bacteria_dir 
  for filename in os.listdir(bacteria_dir):
    f = os.path.join(bacteria_dir, filename)
    if os.path.isfile(f) and not is_file_read(b_out, filename):
      tokenize(f, 0, method, b_out, k=k)
      
  # Tokenize all files in phage_dir
    for filename in os.listdir(phage_dir):
      f = os.path.join(phage_dir, filename)
      if os.path.isfile(f) and not is_file_read(p_out, filename):
        tokenize(f, 1, method, p_out, k=k)

"""\
Runs fasta files through tokenizer and adds the label of 1 for phage and
0 for bacteria. Then shuffles the rows in the dataframe and saves to CSV 

Arguments:
  filepath -- str, path to fasta file
  label -- int, 0 for bacteria or 1 for phage
  method -- str, tokenization method of choice
  output_dir -- str, path to directory for output
  k -- int, length of k if using kmer tokenization
"""
def tokenize(filepath, label, method, output_dir, **kwargs):
  sequences = []
  tokens = []

  k = kwargs.get('k', None)
  filename = os.path.basename(filepath)
  filename = filename.split('.')[0]

  # Calculate appropriate max segment length
  if method == 'codon':
    max_length = MAX_TOKENS * 3
  elif method == 'kmer':
    max_length = MAX_TOKENS - (k - 1)
  elif method == 'bpe':
    max_length = None

  # Process data to get sequences of appropriate length 
  df = preprocess_data(filepath, max_length)
  sequences = df['sequence'].values.tolist()

  # Tokenize according to chosen method
  for seq in range(len(sequences)):
    if method == 'codon':
      tokens.append(seq2codon(sequences[seq]))
    elif method == 'kmer':
      tokens.append(seq2kmer(sequences[seq], k))
    elif method == 'bpe':
      tokens.append(seq2bpe(sequences[seq]))
  df['tokenized'] = tokens
  df['label'] = [label] * len(tokens)
  
  # Shuffle and save to csv
  df = df.sample(frac=1).reset_index(drop=True)
  write_csv(filename, df, output_dir)
  return df

"""\
Read fasta file and truncate sequences to appropriate length, returns dataframe

Arguments:
  filepath -- str, path to fasta file
  max_length -- int, maximum sequence length
Returns:
  df -- dataframe, includes the file title, > input line, start position, and sequence
""" 
def preprocess_data(filepath, max_length): 
  records = []

  f = filepath
  if filepath.endswith('.gz'):
    f = gzip.open(filepath, 'rt', encoding='utf-8')

  try:
    for record in SeqIO.parse(f, 'fasta'):
      filename = os.path.basename(filepath)
      name = filename.split('.')[0]
      segment = str(record.name)
      seq = str(record.seq).upper()
      pos = 0 
      # Truncate sequences if longer than max_length
      if max_length != None:
        while len(seq) > max_length:
          records.append(                  # add subsequence up to max_length
            {
              'name': name,
              'segment': segment,
              'start': pos,
              'sequence': seq[:max_length]
            }
          )
          seq = seq[max_length:]           # sequence continuing from max_length
          pos += max_length
      records.append(
          {
            'name': name,
            'segment': segment,
            'start': pos,
            'sequence': seq
          }
      )
  finally:
    df = pd.DataFrame(data=records)
    return df

"""\
Read in sequences and tokens to attach labels and return dataframe

Arguments:
  sequences -- list, original sequences
  tokens -- list, tokenized sequences
  label -- int, 1 for phage or 0 for bacteria
Returns:
  df -- dataframe
""" 
def attach_labels(sequences, tokens, label):
  d = []
  for i in range(len(tokens)):
    d.append(
        {
          'sequence': sequences[i],
          'tokenized': tokens[i],
          'label': label
        }
    )
  df = pd.DataFrame(data=d)
  return df

"""\
Save the given dataframe to two separate csv files:
1. _full.csv includes the name, start position, sequence, tokenized
   sequence, and label.
2. _tokenized.csv includes the tokenized sequence and the label.

Arguments:
  filename -- str, name of file being tokenized
  df -- dataframe, full dataframe of tokenized sequences
  output_dir -- str, path to directory for output
""" 
def write_csv(filename, df, output_dir):
  # df.to_csv(directory + "/" + filename + '_full.csv', encoding='utf-8', index=False)
  tokenized = df[['tokenized', 'label']]
  tokenized.to_csv(output_dir + "/" + filename + '_tokenized.csv', encoding='utf-8', index=False, header=False, sep='\t')
    

## Different tokenization methods

"""\
Convert a sequence to codons

Arguments:
  seq -- str, original sequence
Returns:
  codons -- str, codons separated by space
"""
def seq2codon(seq):
  codon = [seq[i:i+3] for i in range(0,len(seq),3)]
  codons = " ".join(codon)
  return codons

"""\
Convert a sequence to kmers

Arguments:
  seq -- str, original sequence
  k -- int, kmer of length k
Returns:
  kmers -- str, kmers separated by space
"""
def seq2kmer(seq, k):
  kmer = [seq[i:i+k] for i in range(len(seq)+1-k)]
  kmers = " ".join(kmer)
  return kmers

"""\
Convert a sequence to byte pair encodings

Arguments:
  seq -- str, original sequence
Returns:
  output -- str, decoded tokens separated by a space
"""
def seq2bpe(sequence):                         
  tokenizer = PreTrainedTokenizerFast(tokenizer_file="dna_tokenizer.json")
  encoded_input = tokenizer(sequence, return_tensors="pt")
  token_ids = encoded_input.input_ids
  output = " ".join(tokenizer.batch_decode(token_ids)) 
  return output

"""\
Build the vocabulary for the BPE model

Arguments:
  vocab -- str, directory or list of files to build model vocabulary on
"""
def build_vocab(vocab):
  sequences = parse_vocab(vocab)
  tokenizer = Tokenizer(models.BPE())
  # Customize the tokenizer to handle DNA sequences
  tokenizer.normalizer = normalizers.Sequence([normalizers.NFKC()])
  # Train the tokenizer on your DNA sequences
  trainer = trainers.BpeTrainer(VOCAB_SIZE)
  tokenizer.train_from_iterator(sequences, trainer=trainer)
  tokenizer.save("dna_tokenizer.json")

"""\
Parse vocabulary for the BPE model to work with

Arguments:
  vocab -- str, directory or list of files to build model vocabulary on
Returns:
  sequences -- list, sequences to train model on
"""
def parse_vocab(vocab):
  sequences = []
  # If the input vocabulary is a directory
  if os.path.isdir(vocab):
    for filename in os.listdir(vocab):
      f = os.path.join(vocab, filename)
      if os.path.isfile(f):
        if f.endswith('.gz'):
          f = gzip.open(f, 'rt', encoding='utf-8')
        for record in SeqIO.parse(f, 'fasta'):
          # Truncate sequences if longer than max_length
          seq = str(record.seq).upper()
          while len(seq) > MAX_VOCAB_TOKEN_LENGTH:
            sequences.append(seq[:MAX_VOCAB_TOKEN_LENGTH])
            seq = seq[MAX_VOCAB_TOKEN_LENGTH:]  
          sequences.append(seq)
  # If the input vocabulary is a list of files
  elif os.path.isfile(vocab):
    with open(vocab, 'r') as list:
      for f in list:
        f = f.strip()
        filename = os.path.basename(f)
        if os.path.isfile(f):
          if f.endswith('.gz'):
            f = gzip.open(f, 'rt', encoding='utf-8')
          for record in SeqIO.parse(f, 'fasta'):
            seq = str(record.seq).upper()
            # Truncate sequences if longer than max_length
            while len(seq) > MAX_VOCAB_TOKEN_LENGTH:
              sequences.append(seq[:MAX_VOCAB_TOKEN_LENGTH])
              seq = seq[MAX_VOCAB_TOKEN_LENGTH:]  
            sequences.append(seq)
  return sequences

## Keep track of files

'''\
Determines whether or not a file has already been processed by checking
if the output filename exists in the output directory and has a size
greater than 0.

Arguments:
  filename -- str, the full name of the file to check
Returns:
  isfile -- bool, whether the file has already been tokenized
'''
def is_file_read(directory, filename):
  file_path = os.path.join(directory, filename.split('.')[0]  + '_tokenized.csv')
  if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
    return True
  else:
    return False

## MAIN

def main():
  parser = argparse.ArgumentParser()
  # Parameters
  parser.add_argument(
        "--b", default=None, type=str, required=True, help="The input bacteria directory."
    )
  parser.add_argument(
        "--p", default=None, type=str, required=True, help="The input phage directory."
    )
  parser.add_argument(
        "--b_out", default=None, type=str, required=False, help="The first output directory, for bacteria if using both."
    )
  parser.add_argument(
        "--p_out", default=None, type=str, required=False, help="The second output directory, for phage if using both."
    )
  parser.add_argument(
        "--method", default=None, type=str, required=True, help="The tokenization method of choice: kmer, codon, or bpe."
    )
  parser.add_argument(
        "--k", default=None, type=int, required=False, help="Length k for kmer tokenization."
    )
  parser.add_argument(
        "--vocab", default=None, type=str, required=False, help="The directory or list of files to build the model vocabulary, if using bpe."
    )
  args = parser.parse_args()

  # Tokenize files
  read_files(bacteria_dir=args.b, phage_dir=args.p, method=args.method, b_out=args.b_out, p_out=args.p_out, k=args.k, vocab=args.vocab)

if __name__ == "__main__":
    main()