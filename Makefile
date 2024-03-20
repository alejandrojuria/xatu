# Compiler & compiler flags
CC = g++
FC = gfortran
CFLAGS = -O2 -Wall -lm
FFLAGS = -O2 -Wall -Wno-tabs -lm

# Include folders
INCLUDE = -I$(PWD)/include

# Libraries
LIBS = -DARMA_DONT_USE_WRAPPER -L$(PWD) -lxatu -larmadillo -lopenblas -llapack -larpack -fopenmp -lgfortran

# Conditional flags for compilation
ifeq ($(DEBUG), 1)
	CFLAGS = -Wall -lm -g
	FFLAGS = -Wall -Wno-tabs -lm -g
endif
ifeq ($(HDF5), 1)
	CFLAGS += -DARMA_USE_HDF5
	LIBS += -lhdf5
endif

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
	ar rcs libxatu.a $(OBJECTS) 
	
xatu: main/xatu.cpp $(OBJECTS) 
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

%: main/%.cpp $(OBJECTS)
	$(CC) -o bin/$@ $< $(CFLAGS) $(INCLUDE) $(LIBS)

# Compilation steps
# $< refers to first prerequisite and $@ to the target
build/%.o: src/%.cpp
	$(CC) -c $< -o $@ $(CFLAGS) $(INCLUDE) $(LIBS) 

build/%.o: src/%.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(LIBS) $(INCLUDE)

clean:
	rm -f build/*.o bin/* libxatu.a
