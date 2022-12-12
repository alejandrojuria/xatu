# Compiler & compiler flags
CC = g++
CF = gfortran
CFLAGS = -O2 -Wall -lm
FFLAGS = -O2 -Wall -Wno-tabs -lm

# Include folders
INCLUDE = -I$(PWD)/include

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -llapack -lopenblas -larmadillo -fopenmp -lgfortran

# Compilation targets
CC_SRC_FILES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(CC_SRC_FILES))
FC_SRC_FILES := $(wildcard src/*.f90)
OBJECTS_FC := $(patsubst src/%.f90, build/%.o, $(FC_SRC_FILES))
OBJECTS += $(OBJECTS_FC)

# Create folders
dummy_build_folder := $(shell mkdir -p build)
dummy_bin_folder := $(shell mkdir -p bin)

build:	$(OBJECTS)
	
xatu: 	$(OBJECTS) main/xatu.cpp
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LIBS) $(INCLUDE)

exciton: $(OBJECTS) main/gexciton.cpp
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LIBS) $(INCLUDE)

test: $(OBJECTS) main/test.cpp
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LIBS) $(INCLUDE)

# Compilation steps
# $< refers to first prerequisite and $@ to the target
build/%.o: src/%.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(LIBS) $(INCLUDE)

build/%.o: src/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(LIBS) $(INCLUDE) -Wno-tabs

clean:
	rm -f build/*.o bin/*
