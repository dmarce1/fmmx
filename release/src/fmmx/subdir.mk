################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/fmmx/exafmm.cpp \
../src/fmmx/key.cpp \
../src/fmmx/main.cpp \
../src/fmmx/new.cpp \
../src/fmmx/node_client.cpp \
../src/fmmx/node_server.cpp \
../src/fmmx/silo_output.cpp 

OBJS += \
./src/fmmx/exafmm.o \
./src/fmmx/key.o \
./src/fmmx/main.o \
./src/fmmx/new.o \
./src/fmmx/node_client.o \
./src/fmmx/node_server.o \
./src/fmmx/silo_output.o 

CPP_DEPS += \
./src/fmmx/exafmm.d \
./src/fmmx/key.d \
./src/fmmx/main.d \
./src/fmmx/new.d \
./src/fmmx/node_client.d \
./src/fmmx/node_server.d \
./src/fmmx/silo_output.d 


# Each subdirectory must supply rules for building sources it contributes
src/fmmx/%.o: ../src/fmmx/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -O3 -ipo -inline-level=2 -tbb -DNDEBUG -std=c++0x -xHost `pkg-config --cflags hpx_application` -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


