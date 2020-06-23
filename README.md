# Setup guide
## Ubuntu 20.04 LTS native and WSL
Install the required libraries:
```sudo apt install libopenblas-dev, liblapack-dev```

Clone the Armadillo library repository:
```git clone https://gitlab.com/conradsnicta/armadillo-code.git```

To install run:
```cd armadillo-code
cmake .
make install```

While running cmake, if the libraries were installed correctly they should be detected in the setup.
