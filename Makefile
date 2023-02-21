# Compiler & compiler flags
CC = g++
CFLAGS = -O2 -Wall -lm

# Include folders
INCLUDE = -I$(PWD)/include

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -L$(PWD) -lxatu -larmadillo -lopenblas -llapack -fopenmp

# Compilation targets
SRC_FILES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(SRC_FILES))

# Create folders
dummy_build_folder := $(shell mkdir -p build)
dummy_bin_folder := $(shell mkdir -p bin)

build:	$(OBJECTS)
	ar rcs libxatu.a $(OBJECTS)
	
xatu: main/xatu.cpp $(OBJECTS) 
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

ribbon_offset: main/ribbon_offset.cpp $(OBJECTS) 
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

exciton: main/exciton.cpp $(OBJECTS)
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

transition_rate: main/transition.cpp $(OBJECTS)
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

transition_conv: main/transition_conv.cpp $(OBJECTS)
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

# Compilation steps
# $< refers to first prerequisite and $@ to the target
build/%.o: src/%.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) $(LIBS) 

clean:
	rm -f build/*.o bin/* libxatu.a
