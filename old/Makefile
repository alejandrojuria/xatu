CC = g++
CFLAGS = -O2 -Wall -lm
#INCLUDE = -I/usr/local/lib/armadillo-9.800.2/include -I/usr/lib/OpenBLAS/include/
#LIBS = -L/usr/lib/OpenBLAS/lib -lopenblas -lpthread -llapack -fopenmp
#INCLUDE = -I/opt/OpenBLAS/include/  
#LIBS = -L/opt/OpenBLAS/lib
LIBS = -fopenmp -llapack -lopenblas -larmadillo

zigzag: zigzag.cpp
	$(CC) zigzag.cpp -o zigzag $(CFLAGS) $(INCLUDE) $(LIBS)

excitons: libzigzag.cpp excitons.cpp
	$(CC) excitons.cpp libzigzag.cpp -o excitons $(CFLAGS) $(INCLUDE) $(LIBS)

realspace_wf: libzigzag.cpp libexcitons.cpp libwavefunction.cpp
	$(CC) libwavefunction.cpp main_rswf.cpp libexcitons.cpp libzigzag.cpp -o wavefunction $(CFLAGS) $(INCLUDE) $(LIBS)

realspace_wf: libzigzag.cpp libexcitons.cpp libwavefunction.cpp
	$(CC) libwavefunction.cpp main_rswf.cpp libexcitons.cpp libzigzag.cpp -o wavefunction $(CFLAGS) $(INCLUDE) $(LIBS)

bands: lib/libzigzag.cpp lib/libexcitons.cpp lib/main_bands.cpp
	$(CC) lib/main_bands.cpp  lib/libexcitons.cpp lib/libzigzag.cpp -o bands $(CFLAGS) $(INCLUDE) $(LIBS)

invariant: libinvariant.cpp libzigzag.cpp main_invariant.cpp
	$(CC) main_invariant.cpp libinvariant.cpp libzigzag.cpp -o invariant $(CFLAGS) $(INCLUDE) $(LIBS)

invariant_sparse: libinvariant_sparse.cpp libzigzag.cpp main_invariant_sparse.cpp
	$(CC) main_invariant_sparse.cpp libinvariant_sparse.cpp libzigzag.cpp -o invariant_sparse $(CFLAGS) $(INCLUDE) $(LIBS)

invariant_dense: libinvariant_dense.cpp libzigzag.cpp main_invariant_dense.cpp
	$(CC) main_invariant_dense.cpp libinvariant_dense.cpp libzigzag.cpp -o invariant_dense $(CFLAGS) $(INCLUDE) $(LIBS)

all: libzigzag.cpp libexcitons.cpp main.cpp
	$(CC) main.cpp libexcitons.cpp libzigzag.cpp -o main $(CFLAGS) $(INCLUDE) $(LIBS)

spectrum: lib/libzigzag.cpp lib/libexcitons.cpp lib/main_spectrum.cpp
	$(CC) lib/main_spectrum.cpp lib/libexcitons.cpp lib/libzigzag.cpp -o spectrum $(CFLAGS) $(INCLUDE) $(LIBS)

transition: lib/libzigzag.cpp lib/libexcitons.cpp lib/main_transition.cpp
	$(CC) lib/main_transition.cpp lib/libexcitons.cpp lib/libzigzag.cpp -o transition $(CFLAGS) $(INCLUDE) $(LIBS)

hartree: hartree-fock_check.cpp libexcitons.cpp
	$(CC) hartree-fock_check.cpp libexcitons.cpp libzigzag.cpp -o hartree $(CFLAGS) $(INCLUDE) $(LIBS)

class: lib/Exciton_class.cpp lib/Zigzag_class.cpp lib/main_class_exciton.cpp
	$(CC) lib/Exciton_class.cpp lib/Zigzag_class.cpp lib/main_class_exciton.cpp -o class $(CFLAGS) $(INCLUDE) $(LIBS)

class_transition: lib/Exciton_class.cpp lib/Zigzag_class.cpp lib/main_class_transition.cpp
	$(CC) lib/Exciton_class.cpp lib/Zigzag_class.cpp lib/main_class_transition.cpp -o class_transition $(CFLAGS) $(INCLUDE) $(LIBS)

class_zz: lib/Zigzag_class.cpp lib/main_class_zigzag.cpp
	$(CC) lib/Zigzag_class.cpp lib/main_class_zigzag.cpp -o class_zz $(CFLAGS) $(INCLUDE) $(LIBS)

clean:
	rm zigzag excitons test_file test_file_excitons

clean_spenctrum:
	rm spectrum*

clean_bands:
	rm bands*
