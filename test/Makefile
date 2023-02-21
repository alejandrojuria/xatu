# Compiler & compiler flags
CC = g++
CFLAGS = -O2 -Wall -lm

ROOT_DIR := $(shell dirname $(PWD))

# Include folders
INCLUDE = -I$(ROOT_DIR)/include

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -L$(ROOT_DIR) -lxatu -larmadillo -lopenblas -llapack -fopenmp

all: hbn_base hbn_davidson hbn_spin

hbn_base: hbn_base.cpp
	$(CC) -o ./$@.x $< $(CFLAGS) $(INCLUDE) $(LIBS)

hbn_davidson: hbn_davidson.cpp
	$(CC) -o ./$@.x $< $(CFLAGS) $(INCLUDE) $(LIBS)

hbn_spin: hbn_spin.cpp
	$(CC) -o ./$@.x $< $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm -f ./*.x
