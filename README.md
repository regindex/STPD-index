# STPD-index

This repository contains an implementation of the **Suffix Tree Path Decomposition** index (`STPD-index`).

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

The `STPD-index` tool requires
* A Linux or MacOS 64-bit operating system.
* A modern C++11\14 compiler such as `g++` version 4.9 or higher.

### Usage

You can construct the STPD-index by using the `build_store_stpd_index` executable:
```
build_store_stpd_index [options]
Options:
-h          Print usage info.
-i <arg>    Input text file path. (REQUIRED)
-l <arg>    RLZ reference sequence length (if known). (Def. None)
-o <arg>    Output index file path. (REQUIRED)
```
The current implementation is **optimized for the DNA alphabet**; therefore, the input text must contain only DNA characters (A, C, G, T) and should be provided in ASCII format. <br>
Note that the current path decomposition algorithm computes the explicit suffix tree; therefore, the software **has been tested on small input files** up to a few gigabytes in size.

You can query the STPD-index by using the `locate` executable:
```
locate [options]
Options:
-h          Print usage info.
-i <arg>    Input index filepath. (REQUIRED)
-p <arg>    Patterns FASTA file.  (REQUIRED)
-t <arg>    Maximum number of occurrences to report per pattern. (Def. none)
```
This executable runs **locate all occurrences** queries for all patterns in the file specified with the `-p` option. The pattern file must be provided in FASTA format. The `-t` flag allows you to set the maximum number of occurrences to report for each pattern.
The **output is written to a file named after the pattern file**, with the `.occs` extension.

### Run on Example Data

```console
// Construct the STPD-index
./build/sources/stpd-index-src/build_store_stpd_index -i toy_data/yeast.txt -o toy_data/yeast.ci

// Run locate all occurrence queries using the STPD-index
./build/sources/stpd-index-src/locate -i toy_data/yeast.ci -p toy_data/yeast_patt_100.fasta
```

### External resources

Below is a list of **external software resources** used in this software.

* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
* [sux](https://github.com/vigna/sux)

## Reference and citation 

[1] Ruben Becker, Davide Cenzato, Travis Gagie, Sung-Hwan Kim, Ragnar Groot Koerkamp, Giovanni Manzini, Nicola Prezza: Compressing Suffix Trees by Path Decompositions , ArXiv 2025. ([go to the paper](https://arxiv.org/abs/2506.14734))

## Contacts

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davidecenzato Unive email.

## Funding

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe research and innovation programme, project REGINDEX, grant agreement No. 101039208.