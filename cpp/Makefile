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

Utils.o: Utils.h Common.h /usr/include/stdlib.h /usr/include/Availability.h
Utils.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
Utils.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
Utils.o: /usr/include/sys/_symbol_aliasing.h
Utils.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Utils.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
Utils.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Utils.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Utils.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Utils.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Utils.o: /usr/include/mach/i386/_structs.h
Utils.o: /usr/include/sys/_types/_sigaltstack.h
Utils.o: /usr/include/sys/_types/_ucontext.h
Utils.o: /usr/include/sys/_types/_pthread_attr_t.h
Utils.o: /usr/include/sys/_types/_sigset_t.h
Utils.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_uid_t.h
Utils.o: /usr/include/sys/resource.h /usr/include/stdint.h
Utils.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
Utils.o: /usr/include/sys/_types/_int32_t.h
Utils.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
Utils.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
Utils.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
Utils.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
Utils.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
Utils.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Utils.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Utils.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Utils.o: /usr/include/sys/_types/_ct_rune_t.h
Utils.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Utils.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
Utils.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
Utils.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Utils.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Utils.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
Utils.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Utils.o: /usr/include/time.h /usr/include/_structs.h
Utils.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
Utils.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
Utils.o: /usr/include/math.h /usr/include/sys/stat.h
Utils.o: /usr/include/sys/_types/_blkcnt_t.h
Utils.o: /usr/include/sys/_types/_blksize_t.h
Utils.o: /usr/include/sys/_types/_ino_t.h /usr/include/sys/_types/_ino64_t.h
Utils.o: /usr/include/sys/_types/_nlink_t.h /usr/include/sys/_types/_gid_t.h
Utils.o: /usr/include/sys/_types/_s_ifmt.h
Utils.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
Utils.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
Utils.o: logcpp/log.h /usr/include/sys/time.h
Utils.o: /usr/include/sys/_types/_fd_def.h
Utils.o: /usr/include/sys/_types/_suseconds_t.h
Utils.o: /usr/include/sys/_types/_fd_setsize.h
Utils.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Utils.o: /usr/include/sys/_types/_fd_isset.h
Utils.o: /usr/include/sys/_types/_fd_zero.h
Utils.o: /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h
Read.o: Read.h Common.h /usr/include/stdlib.h /usr/include/Availability.h
Read.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
Read.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
Read.o: /usr/include/sys/_symbol_aliasing.h
Read.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Read.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
Read.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Read.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Read.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Read.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Read.o: /usr/include/mach/i386/_structs.h
Read.o: /usr/include/sys/_types/_sigaltstack.h
Read.o: /usr/include/sys/_types/_ucontext.h
Read.o: /usr/include/sys/_types/_pthread_attr_t.h
Read.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h
Read.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
Read.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
Read.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
Read.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
Read.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
Read.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
Read.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
Read.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
Read.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Read.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Read.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Read.o: /usr/include/sys/_types/_ct_rune_t.h
Read.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Read.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
Read.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
Read.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Read.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Read.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
Read.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Read.o: /usr/include/time.h /usr/include/_structs.h
Read.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
Read.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
Read.o: /usr/include/math.h /usr/include/sys/stat.h
Read.o: /usr/include/sys/_types/_blkcnt_t.h
Read.o: /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_ino_t.h
Read.o: /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_nlink_t.h
Read.o: /usr/include/sys/_types/_gid_t.h /usr/include/sys/_types/_s_ifmt.h
Read.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
Read.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h Utils.h
Read.o: logcpp/log.h /usr/include/sys/time.h
Read.o: /usr/include/sys/_types/_fd_def.h
Read.o: /usr/include/sys/_types/_suseconds_t.h
Read.o: /usr/include/sys/_types/_fd_setsize.h
Read.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Read.o: /usr/include/sys/_types/_fd_isset.h
Read.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
Read.o: /usr/include/sys/_select.h
Dataset.o: Dataset.h Common.h /usr/include/stdlib.h
Dataset.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
Dataset.o: /usr/include/_types.h /usr/include/sys/_types.h
Dataset.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
Dataset.o: /usr/include/sys/_posix_availability.h
Dataset.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
Dataset.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
Dataset.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
Dataset.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
Dataset.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
Dataset.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
Dataset.o: /usr/include/sys/_types/_sigaltstack.h
Dataset.o: /usr/include/sys/_types/_ucontext.h
Dataset.o: /usr/include/sys/_types/_pthread_attr_t.h
Dataset.o: /usr/include/sys/_types/_sigset_t.h
Dataset.o: /usr/include/sys/_types/_size_t.h /usr/include/sys/_types/_uid_t.h
Dataset.o: /usr/include/sys/resource.h /usr/include/stdint.h
Dataset.o: /usr/include/sys/_types/_int8_t.h
Dataset.o: /usr/include/sys/_types/_int16_t.h
Dataset.o: /usr/include/sys/_types/_int32_t.h
Dataset.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
Dataset.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
Dataset.o: /usr/include/_types/_uint64_t.h
Dataset.o: /usr/include/sys/_types/_intptr_t.h
Dataset.o: /usr/include/sys/_types/_uintptr_t.h
Dataset.o: /usr/include/_types/_intmax_t.h /usr/include/_types/_uintmax_t.h
Dataset.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
Dataset.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
Dataset.o: /usr/include/libkern/_OSByteOrder.h
Dataset.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Dataset.o: /usr/include/sys/_types/_ct_rune_t.h
Dataset.o: /usr/include/sys/_types/_rune_t.h
Dataset.o: /usr/include/sys/_types/_wchar_t.h /usr/include/sys/_types/_null.h
Dataset.o: /usr/include/machine/types.h /usr/include/i386/types.h
Dataset.o: /usr/include/sys/_types/___offsetof.h
Dataset.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Dataset.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Dataset.o: /usr/include/sys/_types/_off_t.h
Dataset.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
Dataset.o: /usr/include/secure/_common.h /usr/include/time.h
Dataset.o: /usr/include/_structs.h /usr/include/sys/_structs.h
Dataset.o: /usr/include/sys/_types/_timespec.h
Dataset.o: /usr/include/sys/_types/_clock_t.h
Dataset.o: /usr/include/sys/_types/_time_t.h /usr/include/math.h
Dataset.o: /usr/include/sys/stat.h /usr/include/sys/_types/_blkcnt_t.h
Dataset.o: /usr/include/sys/_types/_blksize_t.h
Dataset.o: /usr/include/sys/_types/_ino_t.h
Dataset.o: /usr/include/sys/_types/_ino64_t.h
Dataset.o: /usr/include/sys/_types/_nlink_t.h
Dataset.o: /usr/include/sys/_types/_gid_t.h /usr/include/sys/_types/_s_ifmt.h
Dataset.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
Dataset.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h Utils.h
Dataset.o: logcpp/log.h /usr/include/sys/time.h
Dataset.o: /usr/include/sys/_types/_fd_def.h
Dataset.o: /usr/include/sys/_types/_suseconds_t.h
Dataset.o: /usr/include/sys/_types/_fd_setsize.h
Dataset.o: /usr/include/sys/_types/_fd_set.h
Dataset.o: /usr/include/sys/_types/_fd_clr.h
Dataset.o: /usr/include/sys/_types/_fd_isset.h
Dataset.o: /usr/include/sys/_types/_fd_zero.h
Dataset.o: /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h
Dataset.o: Read.h
Edge.o: Common.h /usr/include/stdlib.h /usr/include/Availability.h
Edge.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
Edge.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
Edge.o: /usr/include/sys/_symbol_aliasing.h
Edge.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Edge.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
Edge.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Edge.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Edge.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Edge.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Edge.o: /usr/include/mach/i386/_structs.h
Edge.o: /usr/include/sys/_types/_sigaltstack.h
Edge.o: /usr/include/sys/_types/_ucontext.h
Edge.o: /usr/include/sys/_types/_pthread_attr_t.h
Edge.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h
Edge.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
Edge.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
Edge.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
Edge.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
Edge.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
Edge.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
Edge.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
Edge.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
Edge.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Edge.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Edge.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Edge.o: /usr/include/sys/_types/_ct_rune_t.h
Edge.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Edge.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
Edge.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
Edge.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Edge.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Edge.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
Edge.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Edge.o: /usr/include/time.h /usr/include/_structs.h
Edge.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
Edge.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
Edge.o: /usr/include/math.h /usr/include/sys/stat.h
Edge.o: /usr/include/sys/_types/_blkcnt_t.h
Edge.o: /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_ino_t.h
Edge.o: /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_nlink_t.h
Edge.o: /usr/include/sys/_types/_gid_t.h /usr/include/sys/_types/_s_ifmt.h
Edge.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
Edge.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h Utils.h
Edge.o: logcpp/log.h /usr/include/sys/time.h
Edge.o: /usr/include/sys/_types/_fd_def.h
Edge.o: /usr/include/sys/_types/_suseconds_t.h
Edge.o: /usr/include/sys/_types/_fd_setsize.h
Edge.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Edge.o: /usr/include/sys/_types/_fd_isset.h
Edge.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
Edge.o: /usr/include/sys/_select.h Edge.h Read.h
OverlapGraph.o: OverlapGraph.h Common.h /usr/include/stdlib.h
OverlapGraph.o: /usr/include/Availability.h
OverlapGraph.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
OverlapGraph.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
OverlapGraph.o: /usr/include/sys/_symbol_aliasing.h
OverlapGraph.o: /usr/include/sys/_posix_availability.h
OverlapGraph.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
OverlapGraph.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
OverlapGraph.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
OverlapGraph.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
OverlapGraph.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
OverlapGraph.o: /usr/include/i386/_mcontext.h
OverlapGraph.o: /usr/include/mach/i386/_structs.h
OverlapGraph.o: /usr/include/sys/_types/_sigaltstack.h
OverlapGraph.o: /usr/include/sys/_types/_ucontext.h
OverlapGraph.o: /usr/include/sys/_types/_pthread_attr_t.h
OverlapGraph.o: /usr/include/sys/_types/_sigset_t.h
OverlapGraph.o: /usr/include/sys/_types/_size_t.h
OverlapGraph.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
OverlapGraph.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
OverlapGraph.o: /usr/include/sys/_types/_int16_t.h
OverlapGraph.o: /usr/include/sys/_types/_int32_t.h
OverlapGraph.o: /usr/include/sys/_types/_int64_t.h
OverlapGraph.o: /usr/include/_types/_uint8_t.h
OverlapGraph.o: /usr/include/_types/_uint16_t.h
OverlapGraph.o: /usr/include/_types/_uint32_t.h
OverlapGraph.o: /usr/include/_types/_uint64_t.h
OverlapGraph.o: /usr/include/sys/_types/_intptr_t.h
OverlapGraph.o: /usr/include/sys/_types/_uintptr_t.h
OverlapGraph.o: /usr/include/_types/_intmax_t.h
OverlapGraph.o: /usr/include/_types/_uintmax_t.h
OverlapGraph.o: /usr/include/sys/_types/_timeval.h
OverlapGraph.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
OverlapGraph.o: /usr/include/sys/_endian.h
OverlapGraph.o: /usr/include/libkern/_OSByteOrder.h
OverlapGraph.o: /usr/include/libkern/i386/_OSByteOrder.h
OverlapGraph.o: /usr/include/alloca.h /usr/include/sys/_types/_ct_rune_t.h
OverlapGraph.o: /usr/include/sys/_types/_rune_t.h
OverlapGraph.o: /usr/include/sys/_types/_wchar_t.h
OverlapGraph.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
OverlapGraph.o: /usr/include/i386/types.h
OverlapGraph.o: /usr/include/sys/_types/___offsetof.h
OverlapGraph.o: /usr/include/sys/_types/_dev_t.h
OverlapGraph.o: /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h
OverlapGraph.o: /usr/include/sys/_types/_va_list.h
OverlapGraph.o: /usr/include/sys/_types/_off_t.h
OverlapGraph.o: /usr/include/sys/_types/_ssize_t.h
OverlapGraph.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
OverlapGraph.o: /usr/include/time.h /usr/include/_structs.h
OverlapGraph.o: /usr/include/sys/_structs.h
OverlapGraph.o: /usr/include/sys/_types/_timespec.h
OverlapGraph.o: /usr/include/sys/_types/_clock_t.h
OverlapGraph.o: /usr/include/sys/_types/_time_t.h /usr/include/math.h
OverlapGraph.o: /usr/include/sys/stat.h /usr/include/sys/_types/_blkcnt_t.h
OverlapGraph.o: /usr/include/sys/_types/_blksize_t.h
OverlapGraph.o: /usr/include/sys/_types/_ino_t.h
OverlapGraph.o: /usr/include/sys/_types/_ino64_t.h
OverlapGraph.o: /usr/include/sys/_types/_nlink_t.h
OverlapGraph.o: /usr/include/sys/_types/_gid_t.h
OverlapGraph.o: /usr/include/sys/_types/_s_ifmt.h
OverlapGraph.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
OverlapGraph.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
OverlapGraph.o: Utils.h logcpp/log.h /usr/include/sys/time.h
OverlapGraph.o: /usr/include/sys/_types/_fd_def.h
OverlapGraph.o: /usr/include/sys/_types/_suseconds_t.h
OverlapGraph.o: /usr/include/sys/_types/_fd_setsize.h
OverlapGraph.o: /usr/include/sys/_types/_fd_set.h
OverlapGraph.o: /usr/include/sys/_types/_fd_clr.h
OverlapGraph.o: /usr/include/sys/_types/_fd_isset.h
OverlapGraph.o: /usr/include/sys/_types/_fd_zero.h
OverlapGraph.o: /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h
OverlapGraph.o: Dataset.h Read.h HashTable.h Edge.h CS2/cs2.h
OverlapGraph.o: /usr/include/string.h /usr/include/sys/_types/_rsize_t.h
OverlapGraph.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
OverlapGraph.o: /usr/include/secure/_string.h /usr/include/assert.h
OverlapGraph.o: CS2/types_cs2.h CS2/parser_cs2.h CS2/timer.h
OverlapGraph.o: /usr/include/unistd.h /usr/include/sys/unistd.h
OverlapGraph.o: /usr/include/sys/_types/_posix_vdisable.h
OverlapGraph.o: /usr/include/sys/_types/_seek_set.h
OverlapGraph.o: /usr/include/sys/_types/_useconds_t.h
OverlapGraph.o: /usr/include/sys/select.h /usr/include/sys/_types/_uuid_t.h
OverlapGraph.o: /usr/include/gethostuuid.h
HashTable.o: Common.h /usr/include/stdlib.h /usr/include/Availability.h
HashTable.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
HashTable.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
HashTable.o: /usr/include/sys/_symbol_aliasing.h
HashTable.o: /usr/include/sys/_posix_availability.h
HashTable.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
HashTable.o: /usr/include/sys/wait.h /usr/include/sys/_types/_pid_t.h
HashTable.o: /usr/include/sys/_types/_id_t.h /usr/include/sys/signal.h
HashTable.o: /usr/include/sys/appleapiopts.h /usr/include/machine/signal.h
HashTable.o: /usr/include/i386/signal.h /usr/include/machine/_mcontext.h
HashTable.o: /usr/include/i386/_mcontext.h /usr/include/mach/i386/_structs.h
HashTable.o: /usr/include/sys/_types/_sigaltstack.h
HashTable.o: /usr/include/sys/_types/_ucontext.h
HashTable.o: /usr/include/sys/_types/_pthread_attr_t.h
HashTable.o: /usr/include/sys/_types/_sigset_t.h
HashTable.o: /usr/include/sys/_types/_size_t.h
HashTable.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
HashTable.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
HashTable.o: /usr/include/sys/_types/_int16_t.h
HashTable.o: /usr/include/sys/_types/_int32_t.h
HashTable.o: /usr/include/sys/_types/_int64_t.h
HashTable.o: /usr/include/_types/_uint8_t.h /usr/include/_types/_uint16_t.h
HashTable.o: /usr/include/_types/_uint32_t.h /usr/include/_types/_uint64_t.h
HashTable.o: /usr/include/sys/_types/_intptr_t.h
HashTable.o: /usr/include/sys/_types/_uintptr_t.h
HashTable.o: /usr/include/_types/_intmax_t.h /usr/include/_types/_uintmax_t.h
HashTable.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
HashTable.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
HashTable.o: /usr/include/libkern/_OSByteOrder.h
HashTable.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
HashTable.o: /usr/include/sys/_types/_ct_rune_t.h
HashTable.o: /usr/include/sys/_types/_rune_t.h
HashTable.o: /usr/include/sys/_types/_wchar_t.h
HashTable.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
HashTable.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
HashTable.o: /usr/include/sys/_types/_dev_t.h
HashTable.o: /usr/include/sys/_types/_mode_t.h /usr/include/stdio.h
HashTable.o: /usr/include/sys/_types/_va_list.h
HashTable.o: /usr/include/sys/_types/_off_t.h
HashTable.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
HashTable.o: /usr/include/secure/_common.h /usr/include/time.h
HashTable.o: /usr/include/_structs.h /usr/include/sys/_structs.h
HashTable.o: /usr/include/sys/_types/_timespec.h
HashTable.o: /usr/include/sys/_types/_clock_t.h
HashTable.o: /usr/include/sys/_types/_time_t.h /usr/include/math.h
HashTable.o: /usr/include/sys/stat.h /usr/include/sys/_types/_blkcnt_t.h
HashTable.o: /usr/include/sys/_types/_blksize_t.h
HashTable.o: /usr/include/sys/_types/_ino_t.h
HashTable.o: /usr/include/sys/_types/_ino64_t.h
HashTable.o: /usr/include/sys/_types/_nlink_t.h
HashTable.o: /usr/include/sys/_types/_gid_t.h
HashTable.o: /usr/include/sys/_types/_s_ifmt.h
HashTable.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
HashTable.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
HashTable.o: Utils.h logcpp/log.h /usr/include/sys/time.h
HashTable.o: /usr/include/sys/_types/_fd_def.h
HashTable.o: /usr/include/sys/_types/_suseconds_t.h
HashTable.o: /usr/include/sys/_types/_fd_setsize.h
HashTable.o: /usr/include/sys/_types/_fd_set.h
HashTable.o: /usr/include/sys/_types/_fd_clr.h
HashTable.o: /usr/include/sys/_types/_fd_isset.h
HashTable.o: /usr/include/sys/_types/_fd_zero.h
HashTable.o: /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h
HashTable.o: HashTable.h Dataset.h Read.h
main.o: Common.h /usr/include/stdlib.h /usr/include/Availability.h
main.o: /usr/include/AvailabilityInternal.h /usr/include/_types.h
main.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
main.o: /usr/include/sys/_symbol_aliasing.h
main.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
main.o: /usr/include/i386/_types.h /usr/include/sys/wait.h
main.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
main.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
main.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
main.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
main.o: /usr/include/mach/i386/_structs.h
main.o: /usr/include/sys/_types/_sigaltstack.h
main.o: /usr/include/sys/_types/_ucontext.h
main.o: /usr/include/sys/_types/_pthread_attr_t.h
main.o: /usr/include/sys/_types/_sigset_t.h /usr/include/sys/_types/_size_t.h
main.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
main.o: /usr/include/stdint.h /usr/include/sys/_types/_int8_t.h
main.o: /usr/include/sys/_types/_int16_t.h /usr/include/sys/_types/_int32_t.h
main.o: /usr/include/sys/_types/_int64_t.h /usr/include/_types/_uint8_t.h
main.o: /usr/include/_types/_uint16_t.h /usr/include/_types/_uint32_t.h
main.o: /usr/include/_types/_uint64_t.h /usr/include/sys/_types/_intptr_t.h
main.o: /usr/include/sys/_types/_uintptr_t.h /usr/include/_types/_intmax_t.h
main.o: /usr/include/_types/_uintmax_t.h /usr/include/sys/_types/_timeval.h
main.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
main.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
main.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
main.o: /usr/include/sys/_types/_ct_rune_t.h
main.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
main.o: /usr/include/sys/_types/_null.h /usr/include/machine/types.h
main.o: /usr/include/i386/types.h /usr/include/sys/_types/___offsetof.h
main.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
main.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
main.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
main.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
main.o: /usr/include/time.h /usr/include/_structs.h
main.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
main.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
main.o: /usr/include/math.h /usr/include/sys/stat.h
main.o: /usr/include/sys/_types/_blkcnt_t.h
main.o: /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_ino_t.h
main.o: /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_nlink_t.h
main.o: /usr/include/sys/_types/_gid_t.h /usr/include/sys/_types/_s_ifmt.h
main.o: /usr/include/sys/_types/_filesec_t.h /usr/include/ctype.h
main.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h Utils.h
main.o: logcpp/log.h /usr/include/sys/time.h
main.o: /usr/include/sys/_types/_fd_def.h
main.o: /usr/include/sys/_types/_suseconds_t.h
main.o: /usr/include/sys/_types/_fd_setsize.h
main.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
main.o: /usr/include/sys/_types/_fd_isset.h
main.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
main.o: /usr/include/sys/_select.h Read.h Dataset.h Edge.h HashTable.h
main.o: OverlapGraph.h /usr/include/unistd.h /usr/include/sys/unistd.h
main.o: /usr/include/sys/_types/_posix_vdisable.h
main.o: /usr/include/sys/_types/_seek_set.h
main.o: /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h
main.o: /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h
