# Compiler and Flags
CC = mpicc
CFLAGS = -Wall -Isrc -Isrc/simulation_parameters $(HDF5_FLAGS)  -O2 -std=c99 -I/usr/local/include/mpi
LIBS = $(HDF5_LIBS) -lm

NP ?= 30  # default number of processes

# HDF5 library location
#HDF5_FLAGS = -I/opt/homebrew/opt/hdf5/include
#HDF5_LIBS = -L/opt/homebrew/opt/hdf5/lib -lhdf5

HDF5_FLAGS = -I/astro/local/hdf5/1.14.3/intel_oneapi/include
HDF5_LIBS = -L/astro/local/hdf5/1.14.3/intel_oneapi/lib -lhdf5_hl -lhdf5

# Source Files
SRC = $(wildcard src/*.c src/*/*.c src/*/*/*.c src/*/*/*/*.c)

# Targets
TARGET = main.o

compile: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRC) $(LIBS) -o $@

execute:
	mpiexec -n $(NP) ./$(TARGET)

all: compile execute

new: clean all

clean:
	rm -f $(TARGET)