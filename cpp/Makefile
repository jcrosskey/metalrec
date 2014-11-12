#  -c : not to run the linker
#  -W -Wall -pedantic : maximal sensitivity of compiler against possible problems
#  -Wno-variadic-macros : do not warn variadic macros usage (from C90 was turned on by default)
#  -fmessage-length=0 : no line wrapping for error messages, each error message appears on a single line
#  -MP : add a phony target for each dependency other than the main file, causing each to depend on nothing
#  -MMD : dependency and header file options

# compiler
CC = g++

# compile time flags
CFLAGS = -O0 -g3 -W -Wall -Wno-long-long -pedantic -Wno-variadic-macros -fmessage-length=0 -MMD -MP

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
SRCS+=Utils.cpp Read.cpp Dataset.cpp main.cpp 
# use macro to get object files
OBJS = $(SRCS:.cpp=.o)

# program name
PROG = metalrec

.PHONY: default all clean

default: all

all: $(PROG)
	@echo $(PROG) has been compiled

$(PROG): $(OBJS)
	@echo $(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(PROG) $(OBJS) $(LFLAGS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(PROG)

depend:	$(SRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE THIS LINE -- make depend needs it

Utils.o: Utils.h Common.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iomanip
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/istream
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ostream
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ios
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iosfwd
Utils.o: /usr/include/wchar.h /usr/include/_types.h /usr/include/sys/_types.h
Utils.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
Utils.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Utils.o: /usr/include/i386/_types.h /usr/include/Availability.h
Utils.o: /usr/include/AvailabilityInternal.h /usr/include/sys/_types/_null.h
Utils.o: /usr/include/sys/_types/_size_t.h
Utils.o: /usr/include/sys/_types/_mbstate_t.h
Utils.o: /usr/include/sys/_types/_ct_rune_t.h
Utils.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdarg.h
Utils.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Utils.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
Utils.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Utils.o: /usr/include/time.h /usr/include/_structs.h
Utils.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
Utils.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
Utils.o: /usr/include/_wctype.h /usr/include/sys/_types/_wint_t.h
Utils.o: /usr/include/_types/_wctype_t.h /usr/include/ctype.h
Utils.o: /usr/include/runetype.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__locale
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/string
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstring
Utils.o: /usr/include/string.h /usr/include/sys/_types/_rsize_t.h
Utils.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
Utils.o: /usr/include/secure/_string.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdio
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwchar
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwctype
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cctype
Utils.o: /usr/include/wctype.h /usr/include/_types/_wctrans_t.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/algorithm
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/initializer_list
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstddef
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stddef.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/type_traits
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/utility
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple_03
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/memory
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/typeinfo
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/exception
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdint
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdint.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/new
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/limits
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__undef_min_max
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iterator
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base_03
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/tuple
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cassert
Utils.o: /usr/include/assert.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/atomic
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/stdexcept
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/mutex
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__mutex_base
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/chrono
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ctime
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ratio
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/climits
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/limits.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/system_error
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cerrno
Utils.o: /usr/include/errno.h /usr/include/sys/errno.h /usr/include/pthread.h
Utils.o: /usr/include/pthread_impl.h /usr/include/sched.h
Utils.o: /usr/include/sys/_types/_pthread_attr_t.h
Utils.o: /usr/include/sys/_types/_pthread_cond_t.h
Utils.o: /usr/include/sys/_types/_pthread_condattr_t.h
Utils.o: /usr/include/sys/_types/_pthread_key_t.h
Utils.o: /usr/include/sys/_types/_pthread_mutex_t.h
Utils.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
Utils.o: /usr/include/sys/_types/_pthread_once_t.h
Utils.o: /usr/include/sys/_types/_pthread_rwlock_t.h
Utils.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
Utils.o: /usr/include/sys/_types/_pthread_t.h
Utils.o: /usr/include/sys/_types/_mach_port_t.h
Utils.o: /usr/include/sys/_types/_sigset_t.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/functional
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_03
Utils.o: /usr/include/locale.h /usr/include/_locale.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/streambuf
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/locale
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdarg
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdlib
Utils.o: /usr/include/stdlib.h /usr/include/sys/wait.h
Utils.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Utils.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Utils.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Utils.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Utils.o: /usr/include/mach/i386/_structs.h
Utils.o: /usr/include/sys/_types/_sigaltstack.h
Utils.o: /usr/include/sys/_types/_ucontext.h /usr/include/sys/_types/_uid_t.h
Utils.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
Utils.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Utils.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Utils.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Utils.o: /usr/include/machine/types.h /usr/include/i386/types.h
Utils.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
Utils.o: /usr/include/sys/_types/_int32_t.h
Utils.o: /usr/include/sys/_types/_int64_t.h
Utils.o: /usr/include/sys/_types/_intptr_t.h
Utils.o: /usr/include/sys/_types/_uintptr_t.h
Utils.o: /usr/include/sys/_types/___offsetof.h
Utils.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Utils.o: /usr/include/nl_types.h /usr/include/sys/types.h
Utils.o: /usr/include/sys/_types/_blkcnt_t.h
Utils.o: /usr/include/sys/_types/_blksize_t.h
Utils.o: /usr/include/sys/_types/_gid_t.h
Utils.o: /usr/include/sys/_types/_in_addr_t.h
Utils.o: /usr/include/sys/_types/_in_port_t.h
Utils.o: /usr/include/sys/_types/_ino_t.h /usr/include/sys/_types/_ino64_t.h
Utils.o: /usr/include/sys/_types/_key_t.h /usr/include/sys/_types/_nlink_t.h
Utils.o: /usr/include/sys/_types/_useconds_t.h
Utils.o: /usr/include/sys/_types/_suseconds_t.h
Utils.o: /usr/include/sys/_types/_fd_def.h
Utils.o: /usr/include/sys/_types/_fd_setsize.h
Utils.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Utils.o: /usr/include/sys/_types/_fd_zero.h
Utils.o: /usr/include/sys/_types/_fd_isset.h
Utils.o: /usr/include/sys/_types/_fd_copy.h
Utils.o: /usr/include/sys/_types/_fsblkcnt_t.h
Utils.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/_types/_nl_item.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/bitset
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__bit_reference
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iostream
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/sstream
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/fstream
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/vector
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__split_buffer
Utils.o: /usr/include/math.h /usr/include/sys/stat.h
Utils.o: /usr/include/sys/_types/_s_ifmt.h
Utils.o: /usr/include/sys/_types/_filesec_t.h
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/map
Utils.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tree
Utils.o: logcpp/log.h /usr/include/sys/time.h /usr/include/sys/_select.h
Read.o: Common.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iomanip
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/istream
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ostream
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ios
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iosfwd
Read.o: /usr/include/wchar.h /usr/include/_types.h /usr/include/sys/_types.h
Read.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
Read.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
Read.o: /usr/include/i386/_types.h /usr/include/Availability.h
Read.o: /usr/include/AvailabilityInternal.h /usr/include/sys/_types/_null.h
Read.o: /usr/include/sys/_types/_size_t.h
Read.o: /usr/include/sys/_types/_mbstate_t.h
Read.o: /usr/include/sys/_types/_ct_rune_t.h
Read.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdarg.h
Read.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Read.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
Read.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
Read.o: /usr/include/time.h /usr/include/_structs.h
Read.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
Read.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
Read.o: /usr/include/_wctype.h /usr/include/sys/_types/_wint_t.h
Read.o: /usr/include/_types/_wctype_t.h /usr/include/ctype.h
Read.o: /usr/include/runetype.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__locale
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/string
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstring
Read.o: /usr/include/string.h /usr/include/sys/_types/_rsize_t.h
Read.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
Read.o: /usr/include/secure/_string.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdio
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwchar
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwctype
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cctype
Read.o: /usr/include/wctype.h /usr/include/_types/_wctrans_t.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/algorithm
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/initializer_list
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstddef
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stddef.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/type_traits
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/utility
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple_03
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/memory
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/typeinfo
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/exception
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdint
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdint.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/new
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/limits
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__undef_min_max
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iterator
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base_03
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/tuple
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cassert
Read.o: /usr/include/assert.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/atomic
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/stdexcept
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/mutex
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__mutex_base
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/chrono
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ctime
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ratio
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/climits
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/limits.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/system_error
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cerrno
Read.o: /usr/include/errno.h /usr/include/sys/errno.h /usr/include/pthread.h
Read.o: /usr/include/pthread_impl.h /usr/include/sched.h
Read.o: /usr/include/sys/_types/_pthread_attr_t.h
Read.o: /usr/include/sys/_types/_pthread_cond_t.h
Read.o: /usr/include/sys/_types/_pthread_condattr_t.h
Read.o: /usr/include/sys/_types/_pthread_key_t.h
Read.o: /usr/include/sys/_types/_pthread_mutex_t.h
Read.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
Read.o: /usr/include/sys/_types/_pthread_once_t.h
Read.o: /usr/include/sys/_types/_pthread_rwlock_t.h
Read.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
Read.o: /usr/include/sys/_types/_pthread_t.h
Read.o: /usr/include/sys/_types/_mach_port_t.h
Read.o: /usr/include/sys/_types/_sigset_t.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/functional
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_03
Read.o: /usr/include/locale.h /usr/include/_locale.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/streambuf
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/locale
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdarg
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdlib
Read.o: /usr/include/stdlib.h /usr/include/sys/wait.h
Read.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Read.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Read.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Read.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Read.o: /usr/include/mach/i386/_structs.h
Read.o: /usr/include/sys/_types/_sigaltstack.h
Read.o: /usr/include/sys/_types/_ucontext.h /usr/include/sys/_types/_uid_t.h
Read.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
Read.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
Read.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
Read.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Read.o: /usr/include/machine/types.h /usr/include/i386/types.h
Read.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
Read.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
Read.o: /usr/include/sys/_types/_intptr_t.h
Read.o: /usr/include/sys/_types/_uintptr_t.h
Read.o: /usr/include/sys/_types/___offsetof.h
Read.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Read.o: /usr/include/nl_types.h /usr/include/sys/types.h
Read.o: /usr/include/sys/_types/_blkcnt_t.h
Read.o: /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_gid_t.h
Read.o: /usr/include/sys/_types/_in_addr_t.h
Read.o: /usr/include/sys/_types/_in_port_t.h /usr/include/sys/_types/_ino_t.h
Read.o: /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_key_t.h
Read.o: /usr/include/sys/_types/_nlink_t.h
Read.o: /usr/include/sys/_types/_useconds_t.h
Read.o: /usr/include/sys/_types/_suseconds_t.h
Read.o: /usr/include/sys/_types/_fd_def.h
Read.o: /usr/include/sys/_types/_fd_setsize.h
Read.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Read.o: /usr/include/sys/_types/_fd_zero.h
Read.o: /usr/include/sys/_types/_fd_isset.h
Read.o: /usr/include/sys/_types/_fd_copy.h
Read.o: /usr/include/sys/_types/_fsblkcnt_t.h
Read.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/_types/_nl_item.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/bitset
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__bit_reference
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iostream
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/sstream
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/fstream
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/vector
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__split_buffer
Read.o: /usr/include/math.h /usr/include/sys/stat.h
Read.o: /usr/include/sys/_types/_s_ifmt.h
Read.o: /usr/include/sys/_types/_filesec_t.h
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/map
Read.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tree
Read.o: logcpp/log.h /usr/include/sys/time.h /usr/include/sys/_select.h
Read.o: Read.h Utils.h
Dataset.o: Dataset.h Common.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iomanip
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/istream
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ostream
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ios
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iosfwd
Dataset.o: /usr/include/wchar.h /usr/include/_types.h
Dataset.o: /usr/include/sys/_types.h /usr/include/sys/cdefs.h
Dataset.o: /usr/include/sys/_symbol_aliasing.h
Dataset.o: /usr/include/sys/_posix_availability.h
Dataset.o: /usr/include/machine/_types.h /usr/include/i386/_types.h
Dataset.o: /usr/include/Availability.h /usr/include/AvailabilityInternal.h
Dataset.o: /usr/include/sys/_types/_null.h /usr/include/sys/_types/_size_t.h
Dataset.o: /usr/include/sys/_types/_mbstate_t.h
Dataset.o: /usr/include/sys/_types/_ct_rune_t.h
Dataset.o: /usr/include/sys/_types/_rune_t.h
Dataset.o: /usr/include/sys/_types/_wchar_t.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdarg.h
Dataset.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
Dataset.o: /usr/include/sys/_types/_off_t.h
Dataset.o: /usr/include/sys/_types/_ssize_t.h /usr/include/secure/_stdio.h
Dataset.o: /usr/include/secure/_common.h /usr/include/time.h
Dataset.o: /usr/include/_structs.h /usr/include/sys/_structs.h
Dataset.o: /usr/include/sys/_types/_timespec.h
Dataset.o: /usr/include/sys/_types/_clock_t.h
Dataset.o: /usr/include/sys/_types/_time_t.h /usr/include/_wctype.h
Dataset.o: /usr/include/sys/_types/_wint_t.h /usr/include/_types/_wctype_t.h
Dataset.o: /usr/include/ctype.h /usr/include/runetype.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__locale
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/string
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstring
Dataset.o: /usr/include/string.h /usr/include/sys/_types/_rsize_t.h
Dataset.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
Dataset.o: /usr/include/secure/_string.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdio
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwchar
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwctype
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cctype
Dataset.o: /usr/include/wctype.h /usr/include/_types/_wctrans_t.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/algorithm
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/initializer_list
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstddef
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stddef.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/type_traits
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/utility
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple_03
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/memory
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/typeinfo
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/exception
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdint
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdint.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/new
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/limits
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__undef_min_max
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iterator
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base_03
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/tuple
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cassert
Dataset.o: /usr/include/assert.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/atomic
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/stdexcept
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/mutex
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__mutex_base
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/chrono
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ctime
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ratio
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/climits
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/limits.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/system_error
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cerrno
Dataset.o: /usr/include/errno.h /usr/include/sys/errno.h
Dataset.o: /usr/include/pthread.h /usr/include/pthread_impl.h
Dataset.o: /usr/include/sched.h /usr/include/sys/_types/_pthread_attr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_cond_t.h
Dataset.o: /usr/include/sys/_types/_pthread_condattr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_key_t.h
Dataset.o: /usr/include/sys/_types/_pthread_mutex_t.h
Dataset.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_once_t.h
Dataset.o: /usr/include/sys/_types/_pthread_rwlock_t.h
Dataset.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_t.h
Dataset.o: /usr/include/sys/_types/_mach_port_t.h
Dataset.o: /usr/include/sys/_types/_sigset_t.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/functional
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_03
Dataset.o: /usr/include/locale.h /usr/include/_locale.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/streambuf
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/locale
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdarg
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdlib
Dataset.o: /usr/include/stdlib.h /usr/include/sys/wait.h
Dataset.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
Dataset.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
Dataset.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
Dataset.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
Dataset.o: /usr/include/mach/i386/_structs.h
Dataset.o: /usr/include/sys/_types/_sigaltstack.h
Dataset.o: /usr/include/sys/_types/_ucontext.h
Dataset.o: /usr/include/sys/_types/_uid_t.h /usr/include/sys/resource.h
Dataset.o: /usr/include/sys/_types/_timeval.h /usr/include/machine/endian.h
Dataset.o: /usr/include/i386/endian.h /usr/include/sys/_endian.h
Dataset.o: /usr/include/libkern/_OSByteOrder.h
Dataset.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
Dataset.o: /usr/include/machine/types.h /usr/include/i386/types.h
Dataset.o: /usr/include/sys/_types/_int8_t.h
Dataset.o: /usr/include/sys/_types/_int16_t.h
Dataset.o: /usr/include/sys/_types/_int32_t.h
Dataset.o: /usr/include/sys/_types/_int64_t.h
Dataset.o: /usr/include/sys/_types/_intptr_t.h
Dataset.o: /usr/include/sys/_types/_uintptr_t.h
Dataset.o: /usr/include/sys/_types/___offsetof.h
Dataset.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
Dataset.o: /usr/include/nl_types.h /usr/include/sys/types.h
Dataset.o: /usr/include/sys/_types/_blkcnt_t.h
Dataset.o: /usr/include/sys/_types/_blksize_t.h
Dataset.o: /usr/include/sys/_types/_gid_t.h
Dataset.o: /usr/include/sys/_types/_in_addr_t.h
Dataset.o: /usr/include/sys/_types/_in_port_t.h
Dataset.o: /usr/include/sys/_types/_ino_t.h
Dataset.o: /usr/include/sys/_types/_ino64_t.h
Dataset.o: /usr/include/sys/_types/_key_t.h
Dataset.o: /usr/include/sys/_types/_nlink_t.h
Dataset.o: /usr/include/sys/_types/_useconds_t.h
Dataset.o: /usr/include/sys/_types/_suseconds_t.h
Dataset.o: /usr/include/sys/_types/_fd_def.h
Dataset.o: /usr/include/sys/_types/_fd_setsize.h
Dataset.o: /usr/include/sys/_types/_fd_set.h
Dataset.o: /usr/include/sys/_types/_fd_clr.h
Dataset.o: /usr/include/sys/_types/_fd_zero.h
Dataset.o: /usr/include/sys/_types/_fd_isset.h
Dataset.o: /usr/include/sys/_types/_fd_copy.h
Dataset.o: /usr/include/sys/_types/_fsblkcnt_t.h
Dataset.o: /usr/include/sys/_types/_fsfilcnt_t.h
Dataset.o: /usr/include/_types/_nl_item.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/bitset
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__bit_reference
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iostream
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/sstream
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/fstream
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/vector
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__split_buffer
Dataset.o: /usr/include/math.h /usr/include/sys/stat.h
Dataset.o: /usr/include/sys/_types/_s_ifmt.h
Dataset.o: /usr/include/sys/_types/_filesec_t.h
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/map
Dataset.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tree
Dataset.o: logcpp/log.h /usr/include/sys/time.h /usr/include/sys/_select.h
Dataset.o: Read.h Utils.h
main.o: Common.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iomanip
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__config
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/istream
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ostream
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ios
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iosfwd
main.o: /usr/include/wchar.h /usr/include/_types.h /usr/include/sys/_types.h
main.o: /usr/include/sys/cdefs.h /usr/include/sys/_symbol_aliasing.h
main.o: /usr/include/sys/_posix_availability.h /usr/include/machine/_types.h
main.o: /usr/include/i386/_types.h /usr/include/Availability.h
main.o: /usr/include/AvailabilityInternal.h /usr/include/sys/_types/_null.h
main.o: /usr/include/sys/_types/_size_t.h
main.o: /usr/include/sys/_types/_mbstate_t.h
main.o: /usr/include/sys/_types/_ct_rune_t.h
main.o: /usr/include/sys/_types/_rune_t.h /usr/include/sys/_types/_wchar_t.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdarg.h
main.o: /usr/include/stdio.h /usr/include/sys/_types/_va_list.h
main.o: /usr/include/sys/_types/_off_t.h /usr/include/sys/_types/_ssize_t.h
main.o: /usr/include/secure/_stdio.h /usr/include/secure/_common.h
main.o: /usr/include/time.h /usr/include/_structs.h
main.o: /usr/include/sys/_structs.h /usr/include/sys/_types/_timespec.h
main.o: /usr/include/sys/_types/_clock_t.h /usr/include/sys/_types/_time_t.h
main.o: /usr/include/_wctype.h /usr/include/sys/_types/_wint_t.h
main.o: /usr/include/_types/_wctype_t.h /usr/include/ctype.h
main.o: /usr/include/runetype.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__locale
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/string
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstring
main.o: /usr/include/string.h /usr/include/sys/_types/_rsize_t.h
main.o: /usr/include/sys/_types/_errno_t.h /usr/include/strings.h
main.o: /usr/include/secure/_string.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdio
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwchar
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cwctype
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cctype
main.o: /usr/include/wctype.h /usr/include/_types/_wctrans_t.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/algorithm
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/initializer_list
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstddef
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stddef.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/type_traits
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/utility
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tuple_03
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/memory
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/typeinfo
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/exception
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdint
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/stdint.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/new
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/limits
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__undef_min_max
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iterator
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_base_03
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/tuple
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cassert
main.o: /usr/include/assert.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/atomic
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/stdexcept
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/mutex
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__mutex_base
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/chrono
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ctime
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/ratio
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/climits
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/lib/clang/6.0/include/limits.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/system_error
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cerrno
main.o: /usr/include/errno.h /usr/include/sys/errno.h /usr/include/pthread.h
main.o: /usr/include/pthread_impl.h /usr/include/sched.h
main.o: /usr/include/sys/_types/_pthread_attr_t.h
main.o: /usr/include/sys/_types/_pthread_cond_t.h
main.o: /usr/include/sys/_types/_pthread_condattr_t.h
main.o: /usr/include/sys/_types/_pthread_key_t.h
main.o: /usr/include/sys/_types/_pthread_mutex_t.h
main.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
main.o: /usr/include/sys/_types/_pthread_once_t.h
main.o: /usr/include/sys/_types/_pthread_rwlock_t.h
main.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
main.o: /usr/include/sys/_types/_pthread_t.h
main.o: /usr/include/sys/_types/_mach_port_t.h
main.o: /usr/include/sys/_types/_sigset_t.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/functional
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__functional_03
main.o: /usr/include/locale.h /usr/include/_locale.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/streambuf
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/locale
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdarg
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/cstdlib
main.o: /usr/include/stdlib.h /usr/include/sys/wait.h
main.o: /usr/include/sys/_types/_pid_t.h /usr/include/sys/_types/_id_t.h
main.o: /usr/include/sys/signal.h /usr/include/sys/appleapiopts.h
main.o: /usr/include/machine/signal.h /usr/include/i386/signal.h
main.o: /usr/include/machine/_mcontext.h /usr/include/i386/_mcontext.h
main.o: /usr/include/mach/i386/_structs.h
main.o: /usr/include/sys/_types/_sigaltstack.h
main.o: /usr/include/sys/_types/_ucontext.h /usr/include/sys/_types/_uid_t.h
main.o: /usr/include/sys/resource.h /usr/include/sys/_types/_timeval.h
main.o: /usr/include/machine/endian.h /usr/include/i386/endian.h
main.o: /usr/include/sys/_endian.h /usr/include/libkern/_OSByteOrder.h
main.o: /usr/include/libkern/i386/_OSByteOrder.h /usr/include/alloca.h
main.o: /usr/include/machine/types.h /usr/include/i386/types.h
main.o: /usr/include/sys/_types/_int8_t.h /usr/include/sys/_types/_int16_t.h
main.o: /usr/include/sys/_types/_int32_t.h /usr/include/sys/_types/_int64_t.h
main.o: /usr/include/sys/_types/_intptr_t.h
main.o: /usr/include/sys/_types/_uintptr_t.h
main.o: /usr/include/sys/_types/___offsetof.h
main.o: /usr/include/sys/_types/_dev_t.h /usr/include/sys/_types/_mode_t.h
main.o: /usr/include/nl_types.h /usr/include/sys/types.h
main.o: /usr/include/sys/_types/_blkcnt_t.h
main.o: /usr/include/sys/_types/_blksize_t.h /usr/include/sys/_types/_gid_t.h
main.o: /usr/include/sys/_types/_in_addr_t.h
main.o: /usr/include/sys/_types/_in_port_t.h /usr/include/sys/_types/_ino_t.h
main.o: /usr/include/sys/_types/_ino64_t.h /usr/include/sys/_types/_key_t.h
main.o: /usr/include/sys/_types/_nlink_t.h
main.o: /usr/include/sys/_types/_useconds_t.h
main.o: /usr/include/sys/_types/_suseconds_t.h
main.o: /usr/include/sys/_types/_fd_def.h
main.o: /usr/include/sys/_types/_fd_setsize.h
main.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
main.o: /usr/include/sys/_types/_fd_zero.h
main.o: /usr/include/sys/_types/_fd_isset.h
main.o: /usr/include/sys/_types/_fd_copy.h
main.o: /usr/include/sys/_types/_fsblkcnt_t.h
main.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/_types/_nl_item.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/bitset
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__bit_reference
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/iostream
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/sstream
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/fstream
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/vector
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__split_buffer
main.o: /usr/include/math.h /usr/include/sys/stat.h
main.o: /usr/include/sys/_types/_s_ifmt.h
main.o: /usr/include/sys/_types/_filesec_t.h
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/map
main.o: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/__tree
main.o: logcpp/log.h /usr/include/sys/time.h /usr/include/sys/_select.h
main.o: Read.h Utils.h Dataset.h
