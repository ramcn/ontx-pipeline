# ONTX pipeline

## Overview

Software pipeline to execute the ONTX hybrid pipeline for acceleration

## Features

* Dataset and work flow to execute a sample hybrid pipeline  


### Compilation From Source
Refer each subdirectories for the dependency. Flappie for example has the below requirement.
* [Cmake](https://cmake.org/) for building
* [CUnit](http://cunit.sourceforge.net/) library for unit testing
* [HDF5](https://www.hdfgroup.org/) library
* [OpenBLAS](https://www.openblas.net/) library for linear algebra

Similarly install the dependency mentioned in the README of spades, darwin and nanopolish. And then proceed forward to compilation as below.


```bash
# Compile flappie basecaller 
cd flappie; make 
# Compile spades assembler 
cd SPAdes-3.31.1; ./spades_compile.sh 
# Compile darwin aligner 
cd darwin; make 
# Compile Nanopolish 
git clone https://github.com/jts/nanopolish; cd nanopolish; make 
```

# Running the pipeline


```bash
git clone https://github.com/ramcn/ontx-pipeline 
cd ontx-pipeline/dataset
export PATH=$PATH:ontx-pipeline/x86-bin
./hybrid-workflow.sh
```
