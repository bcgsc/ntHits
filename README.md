ntHits 
=
ntHits is a method for identifying reapeats in high-throughput DNA sequencing data. 

Compiling ntHits from GitHub
===========================

When installing ntHits from GitHub source the following tools are
required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

To generate the configure script and make files:

	./autogen.sh
 
Compiling ntHits from source
===========================
To compile and install ntHits in /usr/local:

```
$ ./configure
$ make 
$ sudo make install 
```

To install ntHits in a specified directory:

```
$ ./configure --prefix=/opt/ntHits
$ make 
$ make install 
```

ntHits uses OpenMP for parallelization, which requires a modern compiler such as GCC 4.2 or greater. If you have an older compiler, it is best to upgrade your compiler if possible. If you have multiple versions of GCC installed, you can specify a different compiler:

```
$ ./configure CC=gcc-xx CXX=g++-xx 
```

For the best performance of ntHits, pass `-O3` flag:  

```
$ ./configure CFLAGS='-g -O3' CXXFLAGS='-g -O3' 
```


To run ntHits, its executables should be found in your PATH. If you installed ntHits in /opt/ntHits, add /opt/ntHits/bin to your PATH:

```
$ PATH=/opt/ntHits/bin:$PATH
```

Run ntHits
==========
```
nthits [OPTIONS] ... [FILE]
```
Parameters:
  * `-k`,  `--kmer=SIZE`: the length of *k*-mer `[64]`
  * `-j`,  `--threads=N`: use N parallel threads `[16]`
  * `-c`,  `--cutoff=N`: the maximum coverage of *k*-mer in output `[40]`
  * `-p`,  `--pref=STRING`: the prefix for output file name `[repeat]`
  * `FILE`: input file or set of files seperated by space, in fasta, fastq, sam, and bam formats. The files can also be in compressed (`.gz`, `.bz2`, `.xz`) formats . A list of files containing file names in each row can be passed with `@` prefix.
  
For example to run nthits on a test file `reads.fastq` with `k=50`:
```
$ nthits -k50 reads.fastq 
```
As another example, to run nthits on `5` input files file_1.fq.gz, file_2.fa, file_3.sam, file_4.bam, file_5.fq with `k=64` and 32 threads and repeat regions with frequency `>100``:
```
$ nthits -k64 -c100 -j64 file_1.fq.gz file_2.fa file_3.sam file_4.bam file_5.fq
```

If we have a list of input files `lib.in` with input file names in each row and want to run ntHits with `k=144` and 12 threads:
```
$ nthits -k144 -j12 @lib.in 
```
