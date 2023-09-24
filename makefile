# Compiler and Flags
CC = gcc
CFLAGS = -Wall -Isrc -Isrc/shared_files $(HDF5_FLAGS) -O3
LIBS = $(HDF5_LIBS)

# HDF5 library location
HDF5_FLAGS = -I/opt/homebrew/opt/hdf5/include
HDF5_LIBS = -L/opt/homebrew/opt/hdf5/lib -lhdf5

# Source Files
SRC = $(wildcard src/*.c src/*/*.c src/*/*/*.c src/*/*/*/*.c)

# Targets
TARGET = main.o

compile: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRC) $(LIBS) -o $@

execute:
	./$(TARGET)

all: compile execute

new: clean all

clean:
	rm -f $(TARGET)
