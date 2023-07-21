CC = gcc
CFLAGS = -Isrc
LIBS = 

ARRAY_MEMORY_MANAGEMENT_DIR = ./shared_files/array_memory_management/
SNAPSHOT_IO_OPERATIONS_DIR = ./shared_files/snapshot_io_operations/
DERIVATIVES_DIR = ./shared_files/derivatives/

hdf5 =  -I/opt/homebrew/Cellar/hdf5/1.14.1/include -L/opt/homebrew/Cellar/hdf5/1.14.1/lib -lhdf5 # Location of my local hdf5 library
MEMORY_FILES = ${ARRAY_MEMORY_MANAGEMENT_DIR}array_allocation/*.c ${ARRAY_MEMORY_MANAGEMENT_DIR}array_deallocation/*.c
SNAPSHOT_FILES = ${SNAPSHOT_IO_OPERATIONS_DIR}snapshot_saving/*.c ${SNAPSHOT_IO_OPERATIONS_DIR}snapshot_loading/*.c
FILES = main.c #initialize_solar_s.c
DERIVATIVES = ${DERIVATIVES_DIR}spacial_derivatives/*.c ${DERIVATIVES_DIR}time_derivatives/*.c

compile:
	${CC} ${CFLAGS} ${hdf5} -o main.o ${FILES} ${SNAPSHOT_FILES} ${MEMORY_FILES} ${DERIVATIVES} ${LIBS}
execute:
	./main.o

all: compile execute
