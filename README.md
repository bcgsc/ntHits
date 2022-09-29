# ntHits

ntHits is a method for identifying repeats in high-throughput DNA sequencing data.

# Dependencies

- C++ compiler with c++17 and OpenMP support
- [Meson](https://mesonbuild.com/)
- [btllib](https://github.com/bcgsc/btllib)

To install the dependencies using conda, run `conda install -c bioconda -c conda-forge --file requirements.txt`.

ntHits uses [argparse](https://github.com/p-ranav/argparse) for command-line argument parsing which is built-in as a submodule (no further installation required).

# Installation

Download the latest release and run the following command in the project's root directory to create a buildsystem in the `build` folder:

```shell
meson build
```

Then, `cd` into the `build` folder and compile ntHits using:

```shell
ninja
```

This will generate an executable binary `nthits` in the `build` folder.

# Usage

```
Usage: nthits [options] files 

Reports the most frequent k-mers in input files.

Positional arguments:
files           Input files [required]

Optional arguments:
-t --threads    Number of parallel threads [default: 16]
-k --kmer       k-mer length [default: 64]
-h --hashes     Number of hashes to generate per k-mer/spaced seed [default: 4]
-c --cutoff     k-mer cutoff threshold [required]
-o --out     	  Output file's name [required]
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
