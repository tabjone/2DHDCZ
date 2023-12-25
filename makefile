# Compiler and Flags
CC = mpicc
CFLAGS = -Wall -Isrc -Isrc/shared_files $(HDF5_FLAGS) -O3 -std=c99 -I/usr/local/include/mpi
LIBS = $(HDF5_LIBS) -lm

NP = 3  # default number of processes

# HDF5 library location
HDF5_FLAGS = -I/opt/homebrew/opt/hdf5/include
HDF5_LIBS = -L/opt/homebrew/opt/hdf5/lib -lhdf5

# Source Files
SRC = $(filter-out src/solver_3D/%, $(wildcard src/*.c src/*/*.c src/*/*/*.c src/*/*/*/*.c))

# Targets
TARGET = main.o

compile: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(SRC) $(LIBS) -o $@

execute_mpi:
	mpiexec -n $(NP) ./$(TARGET)
execute:
	./$(TARGET)

all_mpi: compile execute_mpi
all: compile execute

new: clean all
new_mpi: clean all_mpi

clean:
	rm -f $(TARGET)

# Test Source Files (including main_tests.c and other dependencies)
TEST_SRC = src/main_tests.c $(filter-out src/main.c, $(SRC))

# Test Target
TARGET_TESTS = main_tests

# Test Compile Target
test_compile: $(TARGET_TESTS)

$(TARGET_TESTS):
	$(CC) $(CFLAGS) $(TEST_SRC) $(LIBS) -o $@

# Test Execution Target
test_execute:
	./$(TARGET_TESTS)

test: test_compile test_execute