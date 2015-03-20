-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir.mk
-include objects.mk


#  -c : not to run the linker
#  -W -Wall -pedantic : maximal sensitivity of compiler against possible problems
#  -Wno-variadic-macros : do not warn variadic macros usage (from C90 was turned on by default)
#  -fmessage-length=0 : no line wrapping for error messages, each error message appears on a single line
#  -MP : add a phony target for each dependency other than the main file, causing each to depend on nothing
#  -MMD : dependency and header file options

# compiler
CC = g++
MPICC = mpic++

# 
ifneq ($(MAKECMSGOALS), clean)
	include $(sources:.cpp=.d)
endif
# compile time flags
CFLAGS = -O0 -g3 -W -Wno-long-long -pedantic -Wno-variadic-macros -fmessage-length=0 -MMD -MP -std=c++11
CFLAGS+=-Wall

# include dirs containing header files
#PATH_TO_CORE=/Users/cjg/Software/seqan/core
#PATH_TO_EXTRAS=/Users/cjg/Software/seqan/extras
#INCLUDES+=-I${PATH_TO_CORE}/include
#INCLUDES+=-I${PATH_TO_EXTRAS}/include
#INCLUDES+=-I/usr/local/include
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/6.0/include
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include

# lib paths in addition to /usr/lib
LFLAGS=

# libs to link into executable:
LIBS=

# source files and object files
# directive for source files
vpath %.cpp ../
vpath %.h ../
#SRCS := $(wildcard ../*.cpp)
SRCS := Utils.cpp Read.cpp Dataset.cpp Edge.cpp OverlapGraph.cpp HashTable.cpp ec.cpp

# use macro to get object files
OBJS = $(SRCS:.cpp=.o)

# program names
PROG = metalrec 
MPI_PROG = mpi_metalrec

.PHONY: default single mpi all clean

default: single

single:	$(PROG)
mpi:	$(MPI_PROG)

all: $(PROG) $(MPI_PROG)
	@echo Done

$(MPI_PROG): $(OBJS) mpi_metalrec.o
	$(MPICC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(MPI_PROG) $(OBJS) $(LIBS) mpi_metalrec.o
	@echo $(MPI_PROG) compiled

$(PROG): $(OBJS) metalrec.o
	$(CC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LIBS) metalrec.o
	@echo $(PROG) compiled

metalrec.o: metalrec.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

mpi_metalrec.o: mpi_metalrec.cpp
	$(MPICC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Either way will work, first one uses old-Fashioned suffix rules; second one uses static pattern rules
#.cpp.o:
$(OBJS): %.o:%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	$(RM) $(OBJS) $(PROG) $(MPI_PROG)

