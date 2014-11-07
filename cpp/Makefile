# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Dataset.cpp \
../Edge.cpp \
../HashTable.cpp \
../OverlapGraph.cpp \
../Read.cpp \
../main.cpp 

OBJS += \
./Dataset.o \
./Edge.o \
./HashTable.o \
./OverlapGraph.o \
./Read.o \
./main.o 

CC = g++
PROG = ec
CLFAGS = -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP

all: ec
omega: $(OBJS)
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
