# default apple compiler
CC=	g++ 

# mpi compiler
MPICC=	mpic++ 

# openmp integrated GCC compiler
MPCC= /usr/local/bin/g++ -fopenmp 

# compile only, no linking
CFLAGS= -c 

all: metalrec_mpi

metalrec_mpi: metalrec_mpi.cpp utils.cpp directoryStructure.cpp
	$(MPICC) -o metalrec_mpi metalrec_mpi.cpp utils.cpp directoryStructure.cpp

clean:
	rm metalrec_mpi *o
