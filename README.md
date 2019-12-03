# ONTX pipeline

## Overview

Software pipeline to execute the ONTX hybrid pipeline for acceleration

## Features

* Dataset and work flow to execute a sample hybrid pipeline  

# Getting Started


```bash
git clone https://github.com/ramcn/ontx-pipeline 
cd ontx-pipeline/dataset
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
# Compile flappie basecaller 
cd flappie; make 
# Compile spades assembler 
cd SPAdes-3.31.1; ./spades_compile.sh 
# Compile darwin aligner 
cd darwin; make 
# Compile Nanopolish 
git clone https://github.com/jts/nanopolish; cd nanopolish; make 
```


