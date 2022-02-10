CC = g++
CFLAGS = -O2 -Wall -lm
# SRC_DIR = $(PWD)/lib
INCLUDE = -I/usr/lib/armadillo/include -I$(PWD)/include
LIBS = -DARMA_DONT_USE_WRAPPER -lopenblas  -llapack -fopenmp
SRC_FILES := $(wildcard src/*.cpp)
OBJECTS := $(patsubst src/%.cpp, build/%.o, $(SRC_FILES))

build: $(OBJECTS)

system: $(OBJECTS) main/system.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBS) -o bin/$@ $^

exciton: $(OBJECTS) main/gexciton.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBS) -o bin/$@ $^

# Compilation steps
# $< refers to first prerequisite and $@ to the target
build/%.o: src/%.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBS) -c $< -o $@

clean:
	rm -f build/*.o bin/*