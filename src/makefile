CC = gcc
CFLAGS = -Isrc
LIBS = 

ARRAY_MEMORY_MANAGEMENT_DIR = ./shared_files/array_memory_management/
SNAPSHOT_IO_OPERATIONS_DIR = ./shared_files/snapshot_io_operations/
DERIVATIVES_DIR = ./shared_files/derivatives/

hdf5 =  -I/opt/homebrew/Cellar/hdf5/1.14.1/include -L/opt/homebrew/Cellar/hdf5/1.14.1/lib -lhdf5 # Location of my local hdf5 library

# shared files
MEMORY_FILES = ${ARRAY_MEMORY_MANAGEMENT_DIR}array_allocation/*.c ${ARRAY_MEMORY_MANAGEMENT_DIR}array_deallocation/*.c
SNAPSHOT_FILES = ${SNAPSHOT_IO_OPERATIONS_DIR}snapshot_saving/*.c ${SNAPSHOT_IO_OPERATIONS_DIR}snapshot_loading/*.c
SOLAR_S_INITIALIZATION = ./shared_files/solar_s_initialization/*.c
DERIVATIVES = ${DERIVATIVES_DIR}spacial_derivatives/*.c ${DERIVATIVES_DIR}time_derivatives/*.c
INTERPOLATION = ./shared_files/interpolation/*.c
EXTRAPOLATION = 

SHARED_FILES = ${MEMORY_FILES} ${SNAPSHOT_FILES} ${DERIVATIVES} ${SOLAR_S_INITIALIZATION} ${INTERPOLATION} ${EXTRAPOLATION}

2D_HD_SOLVER_DIR = ./hd_solver/2D/
HD_SOLVER_FILES = ${2D_HD_SOLVER_DIR}*.c ${2D_HD_SOLVER_DIR}rhs_functions/*.c ${2D_HD_SOLVER_DIR}boundaries/*.c ${2D_HD_SOLVER_DIR}initialization/*.c

FILES = main.c ${HD_SOLVER_FILES} ./shared_files/structs/*.c

compile:
	${CC} ${CFLAGS} ${hdf5} -o main.o ${FILES} ${SHARED_FILES} ${LIBS}
execute:
	./main.o

all: compile execute
