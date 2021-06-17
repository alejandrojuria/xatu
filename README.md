# Setup guide
## Ubuntu 20.04 LTS native and WSL
Install the required libraries:
```sudo apt install libopenblas-dev, liblapack-dev```

Clone the Armadillo library repository:
```git clone https://gitlab.com/conradsnicta/armadillo-code.git```

To install Armadillo run:
```
cd armadillo-code
cmake .
make install
```

While running cmake, if the libraries were installed correctly they should be detected in the setup. Then, to link the libraries during compilation:
```
g++ test.cpp -o test -fopenmp -larmadillo -lopenblas -llapack
```

## General OS
In case that OpenBLAS or LAPACK are not available through repositories, it is always possible to manually download and compile those libraries.
Clone OpenBLAS: ```git clone https://github.com/xianyi/OpenBLAS.git```

To link the libraries we have to specify the directories where they are installed:
```
INCLUDE = -I/dir/armadillo/include -I/another_dir/OpenBLAS/include/
LIBS = -L/another_dir/OpenBLAS/lib
CFLAGS = -fopenmp -O2 -Wall
g++ test.cpp -o test $(INCLUDE) $(LIBS) $(CFLAGS)
```