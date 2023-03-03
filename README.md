# ntHits

ntHits is a tool for efficiently counting and filtering k-mers based on their frequencies.

# Dependencies

- C++ compiler with c++17 and OpenMP support
- [Meson](https://mesonbuild.com/)
- [btllib](https://github.com/bcgsc/btllib) (>=1.5.0, please install from source and update your `CPPFLAGS` and `LDFLAGS`)
- [Catch2](https://github.com/catchorg/Catch2), only for running tests

ntHits uses [argparse](https://github.com/p-ranav/argparse) for command-line argument parsing which is built-in as a submodule (no further installation required).

# Installation

Download the latest release and run the following command in the project's root directory to create a buildsystem in the `build` folder:

```shell
meson setup build
```

Then, `cd` into the `build` folder and compile ntHits using:

```shell
ninja
```

This will generate two binary files in the `build` folder: `nthits` for generating the desired data structure containing the k-mers and if possible, their counts; and `nthits-bfq` for querying the output if it's a (counting) Bloom filter.

# Usage

```
       _          _ _              
 _ __ | |_  /\  /(_) |_ ___       
| '_ \| __|/ /_/ / | __/ __|      
| | | | |_/ __  /| | |_\__ \     
|_| |_|\__\/ /_/ |_|\__|___/  

Usage: nthits --frequencies VAR [--min-count VAR] [--max-count VAR] [--kmer-length VAR] [-h] [--error-rate VAR] [--seeds VAR] [--threads VAR] [--solid] [--long-mode] --out-file VAR out_type files

Filters k-mers based on counts (cmin <= count <= cmax) in input files

Positional arguments:
  out_type              Output format: Bloom filter 'bf', counting Bloom filter ('cbf'), or table ('table') [required]
  files                 Input files [nargs: 0 or more] [required]

Optional arguments:
  -f, --frequencies     Frequency histogram file (e.g. from ntCard) [required]
  -cmin, --min-count    Minimum k-mer count (>=1), ignored if using --solid [default: 1]
  -cmax, --max-count    Maximum k-mer count (<=254) [default: 254]
  -k, --kmer-length     k-mer length, ignored if using spaced seeds (-s) [default: 64]
  -h, --num-hashes      Number of hashes to generate per k-mer/spaced seed [default: 3]
  -p, --error-rate      Target Bloom filter error rate [default: 0.0001]
  -s, --seeds           If specified, use spaced seeds (separate with commas, e.g. 10101,11011) 
  -t, --threads         Number of parallel threads [default: 4]
  --solid               Automatically tune 'cmin' to filter out erroneous k-mers 
  --long-mode           Optimize data reader for long sequences (>5kbp) 
  -v                    Level of details printed to stdout (-v: normal, -vv detailed) 
  -o, --out-file        Output file's name [required]

Copyright 2022 Canada's Michael Smith Genome Science Centre
```

If the output data structure is a Bloom filter (or CBF), it can be queried by either using the `nthits-bfq` tool, or using btllib's API:

## ntHits Bloom Filter Query Tool

```none
Usage: nthits-bfq [-h] [--cbf] [--seeds VAR] [--silent] bf_path

Query tool for ntHits' output Bloom filter

Positional arguments:
  bf_path       Input Bloom filter file [required]

Optional arguments:
  -h, --help    shows help message and exits 
  -v, --version prints version information and exits 
  --cbf         Treat input file as a counting Bloom filter and output k-mer counts 
  -s, --seeds   Spaced seed patterns separated with commas (e.g. 10101,11011) 
  --silent      Don't print logs to stdout 

Copyright 2022 Canada's Michael Smith Genome Science Centre
```

## btllib's Bloom Filter API

C++ example:

```c++
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <string>

int main() {
  btllib::KmerBloomFilter bf(path_to_bloom_filter);
  // or btllib::KmerCountingBloomFilter8 
  std::string kmer = "AGCTATCAGTCGA";
  std::cout << bf.contains(kmer) << std::endl;
  return 0;
}

```

Python example:

```python
import btllib

bf = btllib.KmerBloomFilter(path_to_bloom_filter)
# or btllib.KmerCountingBloomFilter8
kmer = "AGCTATCAGTCGA"
print(bf.contains(kmer))
```

If using spaced seeds, btllib's `BloomFilter` and `CountingBloomFilter` classes should be used instead. In this case, refer to btllib's docs and examples to query the Bloom filters using hashes generated from a `SeedNtHash` object.
