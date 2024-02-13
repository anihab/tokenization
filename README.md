This repository includes Colab notebooks for building and running a Bacterial DNA tokenizer, as well as analyzing and comparing results.

---

### Generate Random Samples

**randomSamples.ipynb** generates a set of random samples and processes files.

### Build Vocabulary

**vocabulary.ipynb** trains a BPE tokenizer on the full bacterial training directory to create the vocabulary json file.

Must define vocabulary size, input data path, and output directory path. The `INPUT_PATH` should be a path to a directory of randomly sampled and processed bacterial fasta files.

### Tokenize

**tokenize.ipynb** tokenizes all files given:

- a directory or list of selected bacteria sequences in csv files
- a directory or list of phage sequences in csv files
- a tokenizer (vocabulary json file)

---

### Statistics

**statistics.ipynb** gives a collection of functions to analyze tokenized output and produce figures.

---

### Old 

A subdirectory of my old tokenization scripts. Implements 3 tokenization methods:

1. Kmer Tokenization ~1*500~500 nucleotides per 500 tokens​
2. Codon Tokenization ~3*500~1500 nucleotides per 500 tokens​
3. Byte Pair Encoding Tokenization ~8*500~4000 nucleotides per 500 tokens
