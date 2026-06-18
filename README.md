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

If you want to compile the competing software used in the benchmarks, please use
~~~~
cmake -DSTPD_BUILD_BENCH_CLI=ON ..
make
~~~~

### Requirements

The `STPD-index` tool requires
* A Linux or MacOS 64-bit operating system.
* A modern C++20 compiler such as `g++` version 10 or higher.

### Usage

You can construct the STPD-index by using the `CLI/stpd-build` executable:
```
stpd-build [options]
Options:
-i <arg>    Input text base path.                   (REQUIRED)
-o <arg>    Output index file path.                 (REQUIRED)
-p <arg>    Phi-function machinery: (r-index|move). (Def. r-index)
-t <arg>    Text oracle variant: (RLZ|bitpacked).   (Def. RLZ)
-b <arg>    Tabulation sequences length.            (Def. 15)
-l <arg>    RLZ reference prefix length.            (Def. None)
-k          Keep temporary files.                   (Def. False)
-h          Print usage info.
```
The current implementation is **optimized for the DNA alphabet**; therefore, the input text must contain only DNA characters (A, C, G, T) and should be provided in ASCII format. <br>
Note that the current path decomposition algorithm computes the explicit suffix tree; therefore, the software **has been tested on small input files** up to a few gigabytes in size.

You can query the STPD-index by using the `CLI/stpd-locate` executable:
```
stpd-locate [options]
Options:
-h          Print usage info.
-i <arg>    Input index filepath. (REQUIRED)
-p <arg>    Patterns FASTA file.  (REQUIRED)
-t <arg>    Maximum number of occurrences to report per pattern. (Def. none)
-c          Check occurrences correctness. (Def. False)
-b          Run queries in benchmark mode. (Def. False)
```
This executable runs **locate all occurrences** queries for all patterns in the file specified with the `-p` option. The pattern file must be provided in FASTA format. The `-t` flag allows you to set the maximum number of occurrences to report for each pattern.
The **output is written to a file named after the pattern file**, with the `.occs` extension.

### Run on Example Data 

```console
// Construct the STPD-index
./build/CLI/stpd-build -i data/yeast.txt -o data/yeast.stpdi

// Run locate all occurrence queries using the STPD-index
./build/CLI/stpd-locate -i data/yeast.stpdi -p data/yeast_patt_100.fasta
```

### External resources

Below is a list of **external software resources** used in this software.

* [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
* [sux](https://github.com/vigna/sux)

## Reference and citation 

[1] Ruben Becker, Davide Cenzato, Travis Gagie, Sung-Hwan Kim, Ragnar Groot Koerkamp, Giovanni Manzini, Nicola Prezza: Compressing Suffix Trees by Path Decompositions, ArXiv 2025. ([go to the paper](https://arxiv.org/abs/2506.14734))

## Contacts

If you notice any bugs, please feel free to report them by opening a Git issue or by contacting us at davide[dot]cenzato[at]unive[dot]it.

## Funding

This project has received funding from the European Research Council (ERC) under the European Union’s Horizon Europe research and innovation programme, project REGINDEX, grant agreement No. 101039208.