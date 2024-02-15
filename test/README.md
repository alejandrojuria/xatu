We provide a battery of regression tests to ease the development of the code. As regression tests, they ensure that the simulations done are
reproducible, i.e. repeating the same simulation with the same parameters returns the same physical values. To this end, the tests are done
using the model files under ```/models```, meaning that any change in these files will results in the tests failing. The exciton configuration files,
on the other hand, are not used by the tests and can be modified freely. 

To use the tests correctly, first build the library with the changes to be tested from the root of the repository:
```
make build
```

Then build and run the test from the ```test/``` folder:
```
cd test/
make tests
./bin/tests.x
```

You will get a string for each test with its description, and whether it was succesful or not, and a final summary with the total amount of passed and failed tests.
The testing is based on the Catch2 library, so we refer to its documentation for further information on how to deploy additional tests. Lastly, the library provides
a very useful ```help``` flag for the executable:
```
./bin/tests.x -h
```

which will list parameters that can be passed to the test executable to tune the testing process (for instance, to run only one test if you know which one is failing).