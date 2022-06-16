# ntHits

ntHits is a method for identifying repeats in high-throughput DNA sequencing data. 

# Installation

Requirements:

- [cmake](https://cmake.org/) version 3.22 or higher.
- C++ compiler with c++17 support.

To build ntHits, download the latest release and run the following command in the project's root directory to create a buildsystem in the `release` folder:

```shell
cmake -S . -B release
```

Then, compile ntHits and its depedencies ([btllib](https://github.com/bcgsc/btllib) and [argparse](https://github.com/p-ranav/argparse)) using:

```shell
cmake --build release --target all
```

This will generate an executable binary `ntHits` in the `release` folder.

# Usage

```
Usage: ntHits [options] files 

Reports the most frequent k-mers in input files.

Positional arguments:
files           Input files [required]

Optional arguments:
-t --threads    Number of parallel threads [default: 16]
-k --kmer       k-mer length [default: 64]
-h --hashes     Number of hashes to generate per k-mer/spaced seed [default: 4]
-c --cutoff     k-mer cutoff threshold [required]
-p --prefix     Output files' prefix [default: "repeat"]
-s --seeds      If specified, use given spaced seeds separated by commas (e.g. 10101,11011)
--outbloom      Output the most frequent k-mers in a Bloom filter [default: false]
--solid         Output the solid k-mers (non-erroneous k-mers) [default: false]
--long-mode     Optimize data reader for long sequences (>5kbp) [default: false]
-b --bit        [default: 16]
-F              
-f              
-r              
Copyright 2019 Canada's Michael Smith Genome Science Centre
```
