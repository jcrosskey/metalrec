# MetaLREC

## OVERVIEW
* A program designed for correcting long reads with high error rates, such as PacBio reads, with the help from high identity short reads (such as Illumina or CCS reads) from the same sample.
* MetaLREC is written in C++ with MPI support, it's recommended now to use MPI instead of running single thread. Later threaded version might be supported.
* A MPI wrapper to run BLASR in parallel as a preprocessing step for error correction is also included.
* Version 1.0.1

## BUILD and RUN
- System requirement: 
    - The program was tested on Linux and MacOS with GNU compiler.
    - samtools is required
    - BLASR is required for the mpi-wrapper of BLASR
- To compile serial version of MetaLREC only:
	type `make`
- To compile MPI version of MetaLREC only:
	type `make mpi`
- To compile both versions of MetaLREC only:
	type `make all`
- To compile MPI wrapper for BLASR:
	type `make -f Make_align`
- Change your compiler in the corresponding makefiles as needed.
- Type `./mpi_align -h` or `./mpi_metalrec -h` for help to run the programs.
- Example for running MPI version to correct all the long reads in the .bam files:
	mpirun -np total_CPU_number -c config_file -od EC
- BAM files need to be sorted and indexed before being used as input to MetaLREC

## Contributor 
* For questions and suggestions, contact JJ Crosskey, jjchai01@gmail.com
