################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/exafmm.cpp \
../src/key.cpp \
../src/main.cpp \
../src/new.cpp \
../src/node_client.cpp \
../src/node_server.cpp \
../src/silo_output.cpp 

OBJS += \
./src/exafmm.o \
./src/key.o \
./src/main.o \
./src/new.o \
./src/node_client.o \
./src/node_server.o \
./src/silo_output.o 

CPP_DEPS += \
./src/exafmm.d \
./src/key.d \
./src/main.d \
./src/new.d \
./src/node_client.d \
./src/node_server.d \
./src/silo_output.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -g -O3 -tbb -ftrapuv -check-pointers=rw `pkg-config --cflags hpx_application` -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


