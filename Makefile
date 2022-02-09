CC = g++
CFLAGS = -O2 -Wall -lm
INCLUDE = -I/usr/lib/armadillo/include -I$(PWD)/lib
LIBS = -DARMA_DONT_USE_WRAPPER -lopenblas  -llapack -fopenmp -L$(PWD)/lib
PARSE = lib/ConfigurationBase.cpp lib/SystemConfiguration.cpp
#INCLUDE = -I/opt/OpenBLAS/include/  
#LIBS = -L/opt/OpenBLAS/lib
#INCLUDE = -L/usr/lib/x86_x64-linux-gnu
#LIBS = -fopenmp -llapack -lopenblas -larmadillo -larpack

exciton: lib/Exciton.cpp lib/Zigzag.cpp main/exciton.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp main/exciton.cpp -o bin/exciton $(CFLAGS) $(INCLUDE) $(LIBS)

transition: lib/Exciton.cpp lib/Zigzag.cpp main/transition.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp main/transition.cpp -o bin/transition $(CFLAGS) $(INCLUDE) $(LIBS)

zigzag: lib/Zigzag.cpp main/zigzag.cpp
	$(CC) lib/Zigzag.cpp main/zigzag.cpp -o bin/zigzag $(CFLAGS) $(INCLUDE) $(LIBS)

width: lib/Zigzag.cpp main/width.cpp
	$(CC) lib/Zigzag.cpp main/width.cpp -o bin/width $(CFLAGS) $(INCLUDE) $(LIBS)

dos: lib/Exciton.cpp lib/Zigzag.cpp main/main_dos.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp main/main_dos.cpp -o bin/dos $(CFLAGS) $(INCLUDE) $(LIBS)

realspace_wf: lib/System.cpp lib/GExciton.cpp lib/libwavefunction.cpp main/rswf.cpp
	$(CC) lib/libwavefunction.cpp main/rswf.cpp lib/GExciton.cpp lib/System.cpp -o bin/realspace_wf $(CFLAGS) $(INCLUDE) $(LIBS)

old_realspace_wf: lib/Zigzag.cpp lib/Exciton.cpp lib/libwavefunction_old.cpp main/rswf_old.cpp
	$(CC) lib/libwavefunction_old.cpp main/rswf_old.cpp lib/Exciton.cpp lib/Zigzag.cpp -o bin/realspace_wf_old $(CFLAGS) $(INCLUDE) $(LIBS)

wfs: utils/calculate_rs_wf.cpp
	$(CC) utils/calculate_rs_wf.cpp -o bin/wfs $(CFLAGS) $(INCLUDE) $(LIBS)

system: lib/Crystal.cpp lib/System.cpp main/system.cpp
	$(CC) lib/ConfigurationBase.cpp lib/SystemConfiguration.cpp lib/Crystal.cpp lib/System.cpp main/system.cpp -o bin/system $(CFLAGS) $(INCLUDE) $(LIBS)

gexciton: lib/System.cpp lib/GExciton.cpp main/gexciton.cpp
	  $(CC) lib/System.cpp lib/GExciton.cpp main/gexciton.cpp -o bin/gexciton $(CFLAGS) $(INCLUDE) $(LIBS)

entanglement: lib/System.cpp lib/GExciton.cpp lib/libwavefunction.cpp main/entanglement.cpp
	$(CC) lib/libwavefunction.cpp main/entanglement.cpp lib/GExciton.cpp lib/System.cpp -o bin/entanglement $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm zigzag excitons test_file test_file_excitons

clean_spenctrum:
	rm spectrum*

clean_bands:
	rm bands*
