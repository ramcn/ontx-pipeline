# ONTX pipeline

## Overview

Software pipeline to execute the ONTX hybrid pipeline for acceleration

## Features

* Dataset and work flow to execute a sample hybrid pipeline  

# Getting Started


```bash
git clone https://github.com/ramcn/ontx-pipeline 
cd ontx-pipeline
export PATH=$PATH:ontx-pipeline/x86-bin
./hybrid-workflow.sh
```

### Compilation From Source
Flappie has the following dependences
* [Cmake](https://cmake.org/) for building
* [CUnit](http://cunit.sourceforge.net/) library for unit testing
* [HDF5](https://www.hdfgroup.org/) library
* [OpenBLAS](https://www.openblas.net/) library for linear algebra



```bash
#  ! It is highly recommended that OpenBLAS is run in single threaded mode
export OPENBLAS_NUM_THREADS=1
#  List available models
flappie --model help
#  Basecall reads directory
flappie reads/ > basecalls.fq
#  Basecall using a different model
flappie --model r941_5mC reads/ > basecalls.fq
#  Output to SAM (not compatible with modification calls)
flappie --format sam reads/ > basecalls.sam
#  Output to BAM (not compatible with modification calls)
flappie --format sam reads | samtools view -Sb - > basecalls.bam
#  Dump trace data
flappie --trace trace.hdf5 reads > basecalls.fq
#  Basecall in parallel
find reads -name \*.fast5 | parallel -P $(nproc) -X flappie > basecalls.fq
#  Dump trace in parallel.  One trace per parallel process.
find reads -name \*.fast5 | parallel -P $(nproc) -X flappie --trace trace_{%}.hdf5 {} > basecalls.fq
```


