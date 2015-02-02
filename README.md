# MetaLREC

## Overview
* A program designed for correcting long reads with high error rates, such as PacBio reads, with the help from high identity short reads (such as Illumina or CCS reads) from the same sample.
* MetaLREC is written in C++ with MPI support, it's recommended now to use MPI instead of running single thread. Later threaded version might be supported.
* A MPI wrapper to run BLASR in parallel as a preprocessing step for error correction is also included.
* Version 1.0

## BUILD
* System requirement. The program was tested on Linux and MacOS with GNU compiler.
* To compile mpi version of MetaLREC only:
	type `make -f Make_mpi`
* To compile MPI wrapper for BLASR:
	type `make -f Make_align`
* Change your compiler in the corresponding makefiles as needed.
* Type `./mpi_align -h` or `./mpi_metalrec -h` for help to run the programs.

### Contributor 
* For questions and suggestions, contact JJ Chai, jjchai01@gmail.com
