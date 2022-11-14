# Compiler & compiler flags
CC = g++
CFLAGS = -O2 -Wall -lm

# Include folders
INCLUDE = -I$(PWD)/include

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -llapack -lopenblas -larmadillo -fopenmp

# Compilation targets
SRC_FILES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(SRC_FILES))

# Create folders
dummy_build_folder := $(shell mkdir -p build)
dummy_bin_folder := $(shell mkdir -p bin)

build:	$(OBJECTS)
	
xatu: 	$(OBJECTS) main/xatu.cpp
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LIBS) $(INCLUDE)

exciton: $(OBJECTS) main/gexciton.cpp
	$(CC) -o bin/$@ $^ $(CFLAGS) $(LIBS) $(INCLUDE)

# Compilation steps
# $< refers to first prerequisite and $@ to the target
build/%.o: src/%.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(LIBS) $(INCLUDE)

clean:
	rm -f build/*.o bin/*
