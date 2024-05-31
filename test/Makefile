# Compiler & compiler flags
CC = g++
CFLAGS = -O2 -Wall -lm

ROOT_DIR := $(shell dirname $(PWD))

# Include folders
INCLUDE = -I$(ROOT_DIR)/include -I$(ROOT_DIR)/test

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -L$(ROOT_DIR) -L$(ROOT_DIR)/test -lxatu -larmadillo -lopenblas -llapack -larpack -fopenmp -lgfortran

# Conditional flags for compilation
ifeq ($(DEBUG), 1)
	CFLAGS = -Wall -lm -g
endif

dummy_bin_folder := $(shell mkdir -p bin)

tests: tests.cpp
	$(CC) -o bin/$@.x $< $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm -f bin/*.x
