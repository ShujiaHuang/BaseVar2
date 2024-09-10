BaseVar2
========


## Installation

Build the source codes step-by-step.


### How to install htslib

1. Download BaseVar from github:

```bash
$ git clone --recursive git@github.com:ShujiaHuang/BaseVar2.git
```

> WARNING: Please try several times if fail to clone the data causing by 
> the network problem.


2. Shift to htscodecs directory and run the following commands: 

```bash

$ cd htslib/htscodecs
$ autoreconf -i
$ ./configure
$ make

```



3. Go back to the upper directory and install main htslib by running the commands below:

```bash

$ cd htslib

$ autoreconf -i
$ ./configure
$ make

```

ignore this error.

```bash
test/test_khash.c: In function 'write_stats_str2int':
test/test_khash.c:53:9: warning: implicit declaration of function 'kh_stats' [-Wimplicit-function-declaration]
   53 |     if (kh_stats(str2int, h, &empty, &deleted, &hist_size, &hist) == 0) {
      |         ^~~~~~~~
test/test_khash.c:53:18: error: 'str2int' undeclared (first use in this function)
   53 |     if (kh_stats(str2int, h, &empty, &deleted, &hist_size, &hist) == 0) {
      |                  ^~~~~~~
test/test_khash.c:53:18: note: each undeclared identifier is reported only once for each function it appears in
make: *** [test/test_khash.o] Error 1
```


Compute resource
----------------

3G-4G for each threah if set -B 200.







imputation methods
------------------

- [2021 Efficient phasing and imputation of low-coverage sequencing data using large reference panels](https://www.nature.com/articles/s41588-020-00756-0)
- [2022 Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes](https://www.nature.com/articles/s41588-023-01438-3)
- [2023 Accurate rare variant phasing of whole-genome and whole-exome sequencing data in the UK Biobank](https://www.nature.com/articles/s41588-023-01415-w)
- <https://odelaneau.github.io/GLIMPSE>
