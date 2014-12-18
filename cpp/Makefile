#  -c : not to run the linker
#  -W -Wall -pedantic : maximal sensitivity of compiler against possible problems
#  -Wno-variadic-macros : do not warn variadic macros usage (from C90 was turned on by default)
#  -fmessage-length=0 : no line wrapping for error messages, each error message appears on a single line
#  -MP : add a phony target for each dependency other than the main file, causing each to depend on nothing
#  -MMD : dependency and header file options

# compiler
CC = g++

# compile time flags
CFLAGS = -O0 -g3 -W -Wno-long-long -pedantic -Wno-variadic-macros -fmessage-length=0 -MMD -MP
CFLAGS+=-Wall

# include dirs containing header files
PATH_TO_CORE=/Users/cjg/Software/seqan/core
PATH_TO_EXTRAS=/Users/cjg/Software/seqan/extras
INCLUDES+=-I${PATH_TO_CORE}/include
INCLUDES+=-I${PATH_TO_EXTRAS}/include
#INCLUDES+=-I/usr/local/include
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/../lib/clang/6.0/include
#INCLUDES+=-I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include

# lib paths in addition to /usr/lib
LFLAGS=

# libs to link into executable:
LIBS=

# source files and object files
SRCS+=Utils.cpp Read.cpp Dataset.cpp Edge.cpp OverlapGraph.cpp HashTable.cpp main.cpp 
# use macro to get object files
OBJS = $(SRCS:.cpp=.o)

# program name
PROG = metalrec

.PHONY: default all clean

default: all

all: $(PROG)
	@echo $(PROG) has been compiled

$(PROG): $(OBJS)
	@echo $(CC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LIBS)
	$(CC) $(LFLAGS) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(PROG)

depend:	$(SRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE THIS LINE -- make depend needs it

Utils.o: Utils.h Common.h logcpp/log.h
Read.o: Read.h Common.h Utils.h logcpp/log.h
Dataset.o: Dataset.h Common.h Utils.h logcpp/log.h Read.h
Edge.o: Common.h Utils.h logcpp/log.h Edge.h Read.h
OverlapGraph.o: OverlapGraph.h Common.h Utils.h logcpp/log.h Dataset.h Read.h
OverlapGraph.o: Edge.h CS2/cs2.h CS2/types_cs2.h CS2/parser_cs2.h CS2/timer.h
HashTable.o: Common.h Utils.h logcpp/log.h HashTable.h Dataset.h Read.h
main.o: Common.h Utils.h logcpp/log.h Read.h Dataset.h Edge.h OverlapGraph.h
