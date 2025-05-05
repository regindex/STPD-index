# STPD-index

## Download and Install

~~~~
git clone https://github.com/regindex/STPD-index
cd STPD-index
mkdir build
cd build
cmake ..
make
~~~~

### Requirements

The \texttt{STPD-index} tool requires
* A Linux or MacOS 64-bit operating system.
* A modern Python 3 release version 3.7 or higher.
* A modern C++11\14 compiler such as `g++` version 4.9 or higher.
* The zlib library installed on your system.

### Usage

```console
// Construct the ST colex- index
./sources/path-decomp-src/stpd_small -i ../data/toy.txt -o ../data/toy.txt.colex_m -c
./sources/stpd-index-src/build_store_stpd_index -i ../data/toy.txt -o ../data/toy_small_index -v small

// Construct the ST colex+- index
./sources/path-decomp-src/stpd_small -i ../data/toy.txt -o ../data/toy.txt.colex_pm -C
./sources/stpd-index-src/build_store_stpd_index -i ../data/toy.txt -o ../data/toy_large_index -v large

// Run locate all ST colex index
./sources/stpd-index-src/locate -i ../data/toy_small_index -p ../data/toy.reads.fasta
```
