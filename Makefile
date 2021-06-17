CC = g++
CFLAGS = -O2 -Wall -lm
INCLUDE = -I/usr/local/lib/armadillo-9.800.2/include
LIBS = -DARMA_DONT_USE_WRAPPER -lopenblas  -llapack -fopenmp
#INCLUDE = -I/opt/OpenBLAS/include/  
#LIBS = -L/opt/OpenBLAS/lib
#INCLUDE = -L/usr/lib/x86_x64-linux-gnu
#LIBS = -fopenmp -llapack -lopenblas -larmadillo -larpack

exciton: lib/Exciton.cpp lib/Zigzag.cpp lib/main_exciton.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp lib/main_exciton.cpp -o exciton $(CFLAGS) $(INCLUDE) $(LIBS)

transition: lib/Exciton.cpp lib/Zigzag.cpp lib/main_transition.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp lib/main_transition.cpp -o transition $(CFLAGS) $(INCLUDE) $(LIBS)

zigzag: lib/Zigzag.cpp lib/main_zigzag.cpp
	$(CC) lib/Zigzag.cpp lib/main_zigzag.cpp -o zigzag $(CFLAGS) $(INCLUDE) $(LIBS)

width: lib/Zigzag.cpp lib/main_width.cpp
	$(CC) lib/Zigzag.cpp lib/main_width.cpp -o width $(CFLAGS) $(INCLUDE) $(LIBS)

dos: lib/Exciton.cpp lib/Zigzag.cpp lib/main_dos.cpp
	$(CC) lib/Exciton.cpp lib/Zigzag.cpp lib/main_dos.cpp -o dos $(CFLAGS) $(INCLUDE) $(LIBS)

realspace_wf: lib/System.cpp lib/GExciton.cpp lib/libwavefunction.cpp lib/main_rswf.cpp
	$(CC) lib/libwavefunction.cpp lib/main_rswf.cpp lib/GExciton.cpp lib/System.cpp -o realspace_wf $(CFLAGS) $(INCLUDE) $(LIBS)

old_realspace_wf: lib/Zigzag.cpp lib/Exciton.cpp lib/libwavefunction_old.cpp lib/main_rswf_old.cpp
	$(CC) lib/libwavefunction_old.cpp lib/main_rswf_old.cpp lib/Exciton.cpp lib/Zigzag.cpp -o realspace_wf_old $(CFLAGS) $(INCLUDE) $(LIBS)

wfs: utils/calculate_rs_wf.cpp
	$(CC) utils/calculate_rs_wf.cpp -o wfs $(CFLAGS) $(INCLUDE) $(LIBS)

system: lib/System.cpp lib/main_system.cpp
	$(CC) lib/System.cpp lib/main_system.cpp -o system $(CFLAGS) $(INCLUDE) $(LIBS)

gexciton: lib/System.cpp lib/GExciton.cpp lib/main_gexciton.cpp
	  $(CC) lib/System.cpp lib/GExciton.cpp lib/main_gexciton.cpp -o gexciton $(CFLAGS) $(INCLUDE) $(LIBS)

entanglement: lib/System.cpp lib/GExciton.cpp lib/libwavefunction.cpp lib/main_entanglement.cpp
	$(CC) lib/libwavefunction.cpp lib/main_entanglement.cpp lib/GExciton.cpp lib/System.cpp -o entanglement $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm zigzag excitons test_file test_file_excitons

clean_spenctrum:
	rm spectrum*

clean_bands:
	rm bands*
