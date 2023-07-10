CC = gcc
CFLAGS = -Isrc
LIBS = 

hdf5 =  -I/opt/homebrew/Cellar/hdf5/1.14.1/include -L/opt/homebrew/Cellar/hdf5/1.14.1/lib -lhdf5 # Location of my local hdf5 library
MEMORY_FILES = ./array_allocation/*.c ./array_deallocation/*.c
SNAPSHOT_FILES = #load_snapshot.c save_snapshot.c
FILES = main.c #initialize_solar_s.c

compile:
	${CC} ${CFLAGS} ${hdf5} -o main.o ${FILES} ${SNAPSHOT_FILES} ${MEMORY_FILES} ${LIBS}
execute:
	./main.o

all: compile execute
