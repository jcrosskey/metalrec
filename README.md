# MetaLREC

## OVERVIEW
* MetaLREC is a program designed for correcting metagenomic long reads with high error rates (around 20%), such as PacBio reads, using high accuracy short reads (such as Illumina or CCS reads) from the same sample.

* MetaLREC is written in C++ with MPI support, it's recommended now to use MPI instead of running single thread. Later threaded version might be supported.

* Before error correction can be performed, mapping results from accurate reads onto the noisy reads needs to be done. We recommend using [BLASR](https://github.com/PacificBiosciences/blasr) for this step. A MPI wrapper to run BLASR in parallel as a preprocessing step for error correction is also included in MetaLREC.

## Current version
* Version 1.0.3

## Setup
- System requirement: 
    - The program was tested on Linux and MacOS with GNU compiler.
    - For MPI support, openMPI needs to be installed.

- Dependencies:
    - samtools is required
    - BLASR is required for the mpi-wrapper of BLASR

- To compile, change directory to MetaLREC/cpp/Debug:
    - To compile serial version of MetaLREC only: type `make`
    - To compile MPI version of MetaLREC only: type `make mpi`
    - To compile both versions of MetaLREC: type `make all`
    - To compile MPI wrapper for BLASR: type `make mpi_align`
    - Change your compiler in the corresponding makefiles as needed.


## A more detailed explanation on how to use MetaLREC and the preprocessing steps
Suppose we have two datasets, one from Illumina and one from PacBio.

* Step 1: Trim, correct, (and merge for tight insert paired end data), remove duplicate and contained reads from Illumina data as one would do for assembly preparation. A sample protocol can be found [here](https://bitbucket.org/omicsbio/omega2).
* Step 2: For the PacBio data, filter out the short reads (e.g. < 1000bps) which won't contribute much to improving assembly results.
* Step 3: 
    - If your data set is really big you might want to split the Illumina reads and PacBio reads into smaller chunks, so that the mapping can be done on these chunks simutaneously. It's better to split the PacBio data set if possible. If the Illumina data set is really big, it's also fine to split them too. After splitting is done, save the split results for Illumina and PacBio reads into two different folders IlluminaDir, and PacbioDir, go to step 4.
    - Otherwise simply use BLASR (or your preferred mapping tools that can accommadate indels and high error rate) to map the Illumina reads to PacBio reads and go to step 5.
    - BAM files need to be sorted and indexed before being used as input to MetaLREC. Full read sequence with soft clip is required, as well as CIGAR string and NM tag in the alignment result. 
```shell
blasr -noRefineAlign -advanceHalf -noSplitSubreads -minMatch 10 \
-sdpTupleSize 7 -minPctIdentity 70 -bestn 10 \
-sam -clipping soft -header \
$Illumina_reads $Pacbio_reads -out $samfile

# This is an exmaple to use 40 threads and 1.5G for each thread to convert sam file to bam file
samtools view -F 4 -bT ${ref} -@ 40 $samfile |  samtools sort -o ${bamfile} -T sorted -@ 40 -m 1500M
samtools index ${bamfile}
```

* Step 4: Run BLASR in parallel using the MPI wrapper, a sample align.config is provided in the package.

```shell
# In this case, the bam files are automatically generated and sorted
mpi_align -np $total_CPU_number -c $align_config -od $out_dir
```

* Step 5: Run error correction, a sample configuration file metalrec.config is provided in the package. Please change the file locations, you can also change the parameter values, though these have been tuned to provide good error correction performance.

```shell
mpi_metalrec -od $out_dir -c $ec_config -log ERROR
metalrec -od $out_dir -c $ec_config -log ERROR
```

## Contributor 
* For questions and suggestions, contact [JJ Crosskey](mailto:crosskey.jj@gmail.com)
