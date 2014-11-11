# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS +=  Utils.cpp Read.cpp Dataset.cpp main.cpp 

OBJS += Utils.o Read.o Dataset.o main.o 

CC = g++
PROG = metalrec
CLFAGS = -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP

all: metalrec
metalrec: $(OBJS)
		$(CC) $(CFLAGS) -o $(PROG) $(OBJS)
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
clean:
	rm -f $(OBJS) $(PROG)

depend:
	makedepend -- $(CFLAGS) -- $(CPP_SRCS)
# DO NOT DELETE

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
Read.o: Common.h /usr/include/stdlib.h /usr/include/Availability.h
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
Read.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
Read.o: logcpp/log.h /usr/include/sys/time.h
Read.o: /usr/include/sys/_types/_fd_def.h
Read.o: /usr/include/sys/_types/_suseconds_t.h
Read.o: /usr/include/sys/_types/_fd_setsize.h
Read.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
Read.o: /usr/include/sys/_types/_fd_isset.h
Read.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
Read.o: /usr/include/sys/_select.h Read.h Utils.h
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
Dataset.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
Dataset.o: logcpp/log.h /usr/include/sys/time.h
Dataset.o: /usr/include/sys/_types/_fd_def.h
Dataset.o: /usr/include/sys/_types/_suseconds_t.h
Dataset.o: /usr/include/sys/_types/_fd_setsize.h
Dataset.o: /usr/include/sys/_types/_fd_set.h
Dataset.o: /usr/include/sys/_types/_fd_clr.h
Dataset.o: /usr/include/sys/_types/_fd_isset.h
Dataset.o: /usr/include/sys/_types/_fd_zero.h
Dataset.o: /usr/include/sys/_types/_fd_copy.h /usr/include/sys/_select.h
Dataset.o: Read.h Utils.h
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
main.o: /usr/include/runetype.h /usr/include/sys/_types/_wint_t.h
main.o: logcpp/log.h /usr/include/sys/time.h
main.o: /usr/include/sys/_types/_fd_def.h
main.o: /usr/include/sys/_types/_suseconds_t.h
main.o: /usr/include/sys/_types/_fd_setsize.h
main.o: /usr/include/sys/_types/_fd_set.h /usr/include/sys/_types/_fd_clr.h
main.o: /usr/include/sys/_types/_fd_isset.h
main.o: /usr/include/sys/_types/_fd_zero.h /usr/include/sys/_types/_fd_copy.h
main.o: /usr/include/sys/_select.h Read.h Utils.h Dataset.h
