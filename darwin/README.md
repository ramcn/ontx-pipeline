## Intro
Darwin is a genome sequence alignment tool which implements a novel dynamic programming technique called GACT
and a novel hashing scheme called D-SOFT. More details refer the paper: https://stanford.edu/~yatisht/pubs/darwin.pdf
This repo is a clone of https://github.com/yatisht/darwin and experiments of adding few heuristics on it.
## Changes done over vanilla darwin
1. Implement edlib heuristic to optimize GACT
2. Modify the command line arguments to pass configuration parametes like bin_size, threshold, first tile size and first tile threshold
3. Add a new hash table implementation called L1 hash table as an alternative for current seed position table method.
4. L2 hash table code changes 
5. Convert the align kernel to OpenCL
6. Initial changes to support Intel FPGA
## Setup
1. Intel or GCC compiler
2. Intel OpenCL setup
## Getting started
Run the following commands to perform a reference-guided assembly of long reads (sampled_reads.fa) aligned to a reference (sample_ref.fa) using darwin. The sample reference is chrI of the yeast genome (sacCer3) and sample reads are generated using PBSIM (https://code.google.com/archive/p/pbsim/) for 20X coverage of the reference. Output alignments are in Multiple Alignment Format (MAF).
```
    $ make // build the host program darwin-xl
    $ aoc -march=emulator -v align.cl  -o align.aocx // to build for emulator
    $ aoc -board=pac_a10 -v align.cl  -o align.aocx // to build for board
## Host program arguments
    $  darwin-xl reference.fasta reads.fasta mode l1l2
    mode is 1 for cpu implementation of darwin-xl, 2 for fpga implementation of darwin-xl.
    l1l2 is to enable or disable l1l2 hash table. 
## Sample runs
    $ env CL_CONTEXT_EMULATOR_DEVICE_INTELFPGA=1   ./darwin-xl data/escherichia_coli.fa data/reads_ecoli.fa 2 0 > out.maf // run in emulator mode
    $ ./darwin-xl data/escherichia_coli.fa data/reads_ecoli.fa 2 0 > out.maf // run on fpga
```
## Dataset

The directory data/ has 2 datasets. A simulated dataset called sample_ref and a real ecoli bacteria dataset.

## Testing precision and recall
    $ ./darwin-xl data/escherichia_coli.fa data/reads_ecoli_long.fa 3 0 > outv1.maf // original implementation
    $ ./darwin-xl data/escherichia_coli.fa data/reads_ecoli_long.fa 1 0 > outv1.maf // modified implementation
    $ ./maf-convert.py sam outv1.maf > outv1.sam
    $ cd tool-eval
    $ cp ../outv1.sam evaluation/reads-simulated/OxfordNanopore-pbsim-observed_last-2d-1k/escherichia_coli/DARWIN-v1.sam 
    $ ./run-evaluation-ecoli.py v1 // displays the below output of precision and recall

		OxfordNanopore-pbsim-observed_last-2d-1k	escherichia_coli	36.0 / 34.5




