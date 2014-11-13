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
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/platform.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_gcc.h
Read.o: /usr/include/inttypes.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_nvcc.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_debug.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_test_system.h
Read.o: /usr/include/unistd.h /usr/include/sys/unistd.h
Read.o: /usr/include/sys/_types/_posix_vdisable.h
Read.o: /usr/include/sys/_types/_seek_set.h
Read.o: /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h
Read.o: /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h
Read.o: /usr/include/dirent.h /usr/include/sys/dirent.h /usr/include/cxxabi.h
Read.o: /usr/include/stddef.h /usr/include/sys/_types/_ptrdiff_t.h
Read.o: /usr/include/sys/_types/_rsize_t.h /usr/include/signal.h
Read.o: /usr/include/sys/_types/_pthread_t.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_helper.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/profiling.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_exception.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_metaprogramming.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_logic.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_control.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_math.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_type.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_enable_if.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_fundamental.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/macro_deprecated.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/builtin_functions.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_metafunctions.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_tags.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_transport.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_comparison.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_conversion.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/array_construct_destruct.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/hosted_type_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/test_system.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_concept.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/boost_preprocessor_subset.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/concept_checking.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_concepts.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_alphabet.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_concept.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_math.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_adapt_builtins.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_bio.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_qualities.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_storage.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_simple_type.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_tabs.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_funcs.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_profile.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_memset.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_aggregate.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/aggregate_concept.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_packed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_bit_compressed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_packed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_bit_compressed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_allocator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_smart_pointer.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_type.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_simple.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_tristate.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/volatile_ptr.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_simple.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_singlepool.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_multipool.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_chunkpool.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_to_std.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_parallelism.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_math.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_functions.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_log_space_value.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_iterator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_concept.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adaptor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_position.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adapt_std.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_iterator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_proxy.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_container.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/container_concept.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_tangle.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_forwards.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_lexical.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_infix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_suffix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_prefix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_array.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_alloc.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_cstyle.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_block.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_packed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_shortcuts.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/iter_concat_virtual.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_concatenator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_concat_direct.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_tight.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_generous.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_owner.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_list.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_string.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_vector.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_array_pointer.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_utils.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_tags.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_macros.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_primitives.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_misc.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_splitting.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_algorithms.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet_expansion.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_iterator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_string.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_functors.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_view.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_reverse.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_shortcuts.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_forwards.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/cstream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/chunk_collector.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_cstyle.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system.h
Read.o: /usr/include/pthread.h /usr/include/pthread_impl.h
Read.o: /usr/include/sched.h /usr/include/sys/_types/_pthread_cond_t.h
Read.o: /usr/include/sys/_types/_pthread_condattr_t.h
Read.o: /usr/include/sys/_types/_pthread_key_t.h
Read.o: /usr/include/sys/_types/_pthread_mutex_t.h
Read.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
Read.o: /usr/include/sys/_types/_pthread_once_t.h
Read.o: /usr/include/sys/_types/_pthread_rwlock_t.h
Read.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
Read.o: /usr/include/sys/_types/_mach_port_t.h /usr/include/errno.h
Read.o: /usr/include/sys/errno.h /usr/include/sys/_types/_errno_t.h
Read.o: /usr/include/semaphore.h /usr/include/sys/types.h
Read.o: /usr/include/sys/_types/_in_addr_t.h
Read.o: /usr/include/sys/_types/_in_port_t.h /usr/include/sys/_types/_key_t.h
Read.o: /usr/include/sys/_types/_fsblkcnt_t.h
Read.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/sys/fcntl.h
Read.o: /usr/include/sys/_types/_o_sync.h /usr/include/sys/_types/_o_dsync.h
Read.o: /usr/include/sys/semaphore.h /usr/include/aio.h
Read.o: /usr/include/sys/aio.h /usr/include/sys/mman.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_forwards.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_forwards.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_mutex.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_sema.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_event.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_thread.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_sync.h
Read.o: /usr/include/fcntl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_async.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_directory.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereaderiterator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereader.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream_algorithms.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/concept_stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_fstream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_iostream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_cstdio.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_sstream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_raw.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_embl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_genbank.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta_align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_cgviz.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_page.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_mapping.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_mmap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_external.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_mmap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_put.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_char_array.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_mmap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/file_stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single_mmap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double_mmap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/read_auto_format.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/is.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/tokenize.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/stream/lexical_cast.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_edit.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_io.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_data.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_simple.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_idmanager.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_edgestump.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_interface.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_directed.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_undirected.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_automaton.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_wordgraph.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_tree.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_fragment.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_hmm.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_vertex.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_outedge.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_adjacency.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_edge.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_property.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_oracle.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_trie.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_bfs.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_dfs.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_drawing.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_parsing.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_utility_parsing.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_mt19937.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/ext_MersenneTwister.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_uniform.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_normal.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_lognormal.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_geometric.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta_kfunc.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_rng_functor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_shuffle.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_union_find.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_heap_tree.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_hmm.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_lis_his.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align_adapt.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_interval_tree.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_scoring.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_fragment.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_aligngraph.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_exact_iterative.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_inexact.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_annotation.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_heap.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gapped_value_type.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_config.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_traceback.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/matrix_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_array.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_array.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gap_anchor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_anchor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_anchor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_metafunctions.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_cols.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_algorithm_tags.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_profile.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_band.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_scout.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_linear.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_affine.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_trace_segment.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_adaptor.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_setup.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_linear.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_affine.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_meta_info.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_sparse.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix_sparse.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_trace_matrix.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_algorithm_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_hirschberg_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_hirschberg_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_waterman_eggert_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded_waterman_eggert_impl.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_unbanded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_banded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_unbanded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_unbanded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_banded.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_specialized.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_operations.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/align/evaluate_alignment.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_name_store_cache.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_io_context.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/cigar.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_header_record.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_sam_conversion.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_tags_dict.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/read_sam.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/write_sam.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record_util.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_reader.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_reader.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_writer.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_writer.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_context.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_stream.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_base.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_read.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_matepair.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_library.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_contig.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_annotation.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_all.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align_intervals.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_intervaltree.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_sam.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_gff.h
Read.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_ucsc.h
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
Dataset.o: Read.h /Users/cjg/Software/seqan/core/include/seqan/sequence.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/platform.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_gcc.h
Dataset.o: /usr/include/inttypes.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_nvcc.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_debug.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_test_system.h
Dataset.o: /usr/include/unistd.h /usr/include/sys/unistd.h
Dataset.o: /usr/include/sys/_types/_posix_vdisable.h
Dataset.o: /usr/include/sys/_types/_seek_set.h
Dataset.o: /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h
Dataset.o: /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h
Dataset.o: /usr/include/dirent.h /usr/include/sys/dirent.h
Dataset.o: /usr/include/cxxabi.h /usr/include/stddef.h
Dataset.o: /usr/include/sys/_types/_ptrdiff_t.h
Dataset.o: /usr/include/sys/_types/_rsize_t.h /usr/include/signal.h
Dataset.o: /usr/include/sys/_types/_pthread_t.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_helper.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/profiling.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_exception.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_metaprogramming.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_logic.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_control.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_math.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_type.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_enable_if.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_fundamental.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/macro_deprecated.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/builtin_functions.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_metafunctions.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_tags.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_transport.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_comparison.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_conversion.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/array_construct_destruct.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/hosted_type_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/test_system.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_concept.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/boost_preprocessor_subset.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/concept_checking.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_concepts.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_alphabet.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_concept.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_math.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_adapt_builtins.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_bio.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_qualities.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_storage.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_simple_type.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_tabs.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_funcs.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_profile.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_memset.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_aggregate.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/aggregate_concept.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_packed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_bit_compressed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_packed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_bit_compressed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_allocator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_smart_pointer.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_type.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_simple.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_tristate.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/volatile_ptr.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_simple.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_singlepool.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_multipool.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_chunkpool.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_to_std.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_parallelism.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_math.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_functions.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_log_space_value.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_iterator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_concept.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adaptor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_position.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adapt_std.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_iterator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_proxy.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_container.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/container_concept.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_tangle.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_forwards.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_lexical.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_infix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_suffix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_prefix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_array.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_alloc.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_cstyle.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_block.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_packed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_shortcuts.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/iter_concat_virtual.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_concatenator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_concat_direct.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_tight.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_generous.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_owner.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_list.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_string.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_vector.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_array_pointer.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_utils.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_tags.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_macros.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_primitives.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_misc.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_splitting.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_algorithms.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet_expansion.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_iterator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_string.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_functors.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_view.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_reverse.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_shortcuts.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_forwards.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/cstream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/chunk_collector.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_cstyle.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system.h
Dataset.o: /usr/include/pthread.h /usr/include/pthread_impl.h
Dataset.o: /usr/include/sched.h /usr/include/sys/_types/_pthread_cond_t.h
Dataset.o: /usr/include/sys/_types/_pthread_condattr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_key_t.h
Dataset.o: /usr/include/sys/_types/_pthread_mutex_t.h
Dataset.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
Dataset.o: /usr/include/sys/_types/_pthread_once_t.h
Dataset.o: /usr/include/sys/_types/_pthread_rwlock_t.h
Dataset.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
Dataset.o: /usr/include/sys/_types/_mach_port_t.h /usr/include/errno.h
Dataset.o: /usr/include/sys/errno.h /usr/include/sys/_types/_errno_t.h
Dataset.o: /usr/include/semaphore.h /usr/include/sys/types.h
Dataset.o: /usr/include/sys/_types/_in_addr_t.h
Dataset.o: /usr/include/sys/_types/_in_port_t.h
Dataset.o: /usr/include/sys/_types/_key_t.h
Dataset.o: /usr/include/sys/_types/_fsblkcnt_t.h
Dataset.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/sys/fcntl.h
Dataset.o: /usr/include/sys/_types/_o_sync.h
Dataset.o: /usr/include/sys/_types/_o_dsync.h /usr/include/sys/semaphore.h
Dataset.o: /usr/include/aio.h /usr/include/sys/aio.h /usr/include/sys/mman.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_forwards.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_forwards.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_mutex.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_sema.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_event.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_thread.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_sync.h
Dataset.o: /usr/include/fcntl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_async.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_directory.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereaderiterator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereader.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream_algorithms.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/concept_stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_fstream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_iostream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_cstdio.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_sstream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_raw.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_embl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_genbank.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta_align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_cgviz.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_page.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_mapping.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_mmap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_external.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_mmap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_put.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_char_array.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_mmap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/file_stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single_mmap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double_mmap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/read_auto_format.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/is.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/tokenize.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/stream/lexical_cast.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_edit.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_io.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_data.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_simple.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_idmanager.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_edgestump.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_interface.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_directed.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_undirected.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_automaton.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_wordgraph.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_tree.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_fragment.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_hmm.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_vertex.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_outedge.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_adjacency.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_edge.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_property.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_oracle.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_trie.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_bfs.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_dfs.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_drawing.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_parsing.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_utility_parsing.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_mt19937.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/ext_MersenneTwister.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_uniform.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_normal.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_lognormal.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_geometric.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta_kfunc.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_rng_functor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_shuffle.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_union_find.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_heap_tree.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_hmm.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_lis_his.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align_adapt.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_interval_tree.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_scoring.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_fragment.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_aligngraph.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_exact_iterative.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_inexact.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_annotation.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_heap.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gapped_value_type.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_config.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_traceback.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/matrix_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_array.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_array.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gap_anchor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_anchor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_anchor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_metafunctions.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_cols.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_algorithm_tags.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_profile.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_band.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_scout.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_linear.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_affine.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_trace_segment.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_adaptor.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_setup.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_linear.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_affine.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_meta_info.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_sparse.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix_sparse.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_trace_matrix.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_algorithm_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_hirschberg_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_hirschberg_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_waterman_eggert_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded_waterman_eggert_impl.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_unbanded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_banded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_unbanded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_unbanded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_banded.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_specialized.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_operations.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/align/evaluate_alignment.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_name_store_cache.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_io_context.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/cigar.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_header_record.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_sam_conversion.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_tags_dict.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/read_sam.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/write_sam.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record_util.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_reader.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_reader.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_writer.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_writer.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_context.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_stream.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_base.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_read.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_matepair.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_library.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_contig.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_annotation.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_all.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align_intervals.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_intervaltree.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_sam.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_gff.h
Dataset.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_ucsc.h
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
main.o: /usr/include/sys/_select.h Read.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/platform.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_gcc.h
main.o: /usr/include/inttypes.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/platform/platform_nvcc.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_debug.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_test_system.h
main.o: /usr/include/unistd.h /usr/include/sys/unistd.h
main.o: /usr/include/sys/_types/_posix_vdisable.h
main.o: /usr/include/sys/_types/_seek_set.h
main.o: /usr/include/sys/_types/_useconds_t.h /usr/include/sys/select.h
main.o: /usr/include/sys/_types/_uuid_t.h /usr/include/gethostuuid.h
main.o: /usr/include/dirent.h /usr/include/sys/dirent.h /usr/include/cxxabi.h
main.o: /usr/include/stddef.h /usr/include/sys/_types/_ptrdiff_t.h
main.o: /usr/include/sys/_types/_rsize_t.h /usr/include/signal.h
main.o: /usr/include/sys/_types/_pthread_t.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/debug_helper.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/profiling.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_exception.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_metaprogramming.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_logic.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_control.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_math.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_type.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/metaprogramming_enable_if.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_fundamental.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/macro_deprecated.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/builtin_functions.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_metafunctions.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_tags.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_transport.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_comparison.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_conversion.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/array_construct_destruct.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/hosted_type_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/test_system.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_concept.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/boost_preprocessor_subset.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/concept_checking.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/fundamental_concepts.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_alphabet.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_concept.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_math.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_adapt_builtins.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_bio.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_qualities.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_storage.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_simple_type.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_tabs.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue_funcs.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_residue.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/alphabet_profile.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_memset.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_aggregate.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/aggregate_concept.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_packed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/pair_bit_compressed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/triple_packed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/tuple_bit_compressed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_allocator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_smart_pointer.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_type.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_simple.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/holder_tristate.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/volatile_ptr.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_simple.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_singlepool.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_multipool.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_chunkpool.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/allocator_to_std.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_parallelism.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_math.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_functions.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/math_log_space_value.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_iterator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_concept.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adaptor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_position.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/iterator_adapt_std.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/proxy_iterator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_proxy.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_container.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/container_concept.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/basic/basic_tangle.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_forwards.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_lexical.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_infix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_suffix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_prefix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_array.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_alloc.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_cstyle.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_block.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_packed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_shortcuts.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/iter_concat_virtual.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/sequence_concatenator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_concat_direct.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_tight.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_dependent_generous.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/string_set_owner.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_list.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_string.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_std_vector.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/adapt_array_pointer.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/sequence/segment_utils.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_tags.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_macros.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_primitives.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_atomic_misc.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_splitting.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/parallel/parallel_algorithms.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_alphabet_expansion.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_iterator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_string.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_functors.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_view.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_reverse.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/modifier/modifier_shortcuts.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_forwards.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/cstream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/chunk_collector.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_cstyle.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system.h
main.o: /usr/include/pthread.h /usr/include/pthread_impl.h
main.o: /usr/include/sched.h /usr/include/sys/_types/_pthread_cond_t.h
main.o: /usr/include/sys/_types/_pthread_condattr_t.h
main.o: /usr/include/sys/_types/_pthread_key_t.h
main.o: /usr/include/sys/_types/_pthread_mutex_t.h
main.o: /usr/include/sys/_types/_pthread_mutexattr_t.h
main.o: /usr/include/sys/_types/_pthread_once_t.h
main.o: /usr/include/sys/_types/_pthread_rwlock_t.h
main.o: /usr/include/sys/_types/_pthread_rwlockattr_t.h
main.o: /usr/include/sys/_types/_mach_port_t.h /usr/include/errno.h
main.o: /usr/include/sys/errno.h /usr/include/sys/_types/_errno_t.h
main.o: /usr/include/semaphore.h /usr/include/sys/types.h
main.o: /usr/include/sys/_types/_in_addr_t.h
main.o: /usr/include/sys/_types/_in_port_t.h /usr/include/sys/_types/_key_t.h
main.o: /usr/include/sys/_types/_fsblkcnt_t.h
main.o: /usr/include/sys/_types/_fsfilcnt_t.h /usr/include/sys/fcntl.h
main.o: /usr/include/sys/_types/_o_sync.h /usr/include/sys/_types/_o_dsync.h
main.o: /usr/include/sys/semaphore.h /usr/include/aio.h
main.o: /usr/include/sys/aio.h /usr/include/sys/mman.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_forwards.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_forwards.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_mutex.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_sema.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_event.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/system_thread.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_sync.h
main.o: /usr/include/fcntl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_async.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/system/file_directory.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereaderiterator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_filereader.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/stream_algorithms.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/concept_stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_fstream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_iostream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_cstdio.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_sstream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_raw.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_embl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_genbank.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_fasta_align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_cgviz.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_page.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_mapping.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_mmap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/string_external.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/file/file_format_mmap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_put.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/stream_char_array.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/adapt_mmap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/file_stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_single_mmap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/record_reader_double_mmap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/read_auto_format.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/is.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/tokenize.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/stream/lexical_cast.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_edit.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_io.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_matrix_data.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/score/score_simple.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_idmanager.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_edgestump.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_interface.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_directed.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_undirected.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_automaton.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_wordgraph.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_tree.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_fragment.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_hmm.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_vertex.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_outedge.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_adjacency.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_edge.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_property.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_oracle.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_impl_trie.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_bfs.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_iterator_dfs.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_drawing.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_parsing.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_types/graph_utility_parsing.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_mt19937.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/ext_MersenneTwister.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_uniform.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_normal.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_lognormal.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_geometric.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_beta_kfunc.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_rng_functor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/random/random_shuffle.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_union_find.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_heap_tree.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_hmm.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_algorithms/graph_algorithm_lis_his.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_impl_align_adapt.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_interval_tree.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_scoring.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_fragment.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_aligngraph.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_exact_iterative.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_inexact.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/graph_align/graph_algorithm_refine_annotation.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/priority_type_heap.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gapped_value_type.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_config.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_traceback.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/matrix_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_array.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_array.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gap_anchor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_anchor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/gaps_iterator_anchor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_metafunctions.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_cols.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/align_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_algorithm_tags.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_profile.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_band.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_scout.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_linear.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_cell_affine.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_trace_segment.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_adaptor.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_setup.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_linear.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_formula_affine.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_meta_info.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_sparse.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_score_matrix_sparse.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_matrix_navigator_trace_matrix.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_traceback_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/dp_algorithm_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_hirschberg_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_myers_hirschberg_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_waterman_eggert_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded_waterman_eggert_impl.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_unbanded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_banded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_unbanded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_banded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_unbanded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/local_alignment_enumeration_banded.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/global_alignment_specialized.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/alignment_operations.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/align/evaluate_alignment.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/misc/misc_name_store_cache.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_io_context.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/cigar.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_header_record.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_sam_conversion.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_tags_dict.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/read_sam.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/write_sam.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_alignment_record_util.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_reader.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_reader.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/xam_writer.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/sam_writer.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/bam_io/bam_stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_context.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_io_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/gff_io/gff_stream.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_base.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_read.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_matepair.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_library.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_contig.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_annotation.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_all.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_align_intervals.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_intervaltree.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_sam.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_gff.h
main.o: /Users/cjg/Software/seqan/core/include/seqan/store/store_io_ucsc.h
main.o: Dataset.h
