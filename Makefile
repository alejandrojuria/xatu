CC = g++
CFLAGS = -O2 -Wall -lm
#INCLUDE = -I/usr/local/lib/armadillo-9.800.2/include -I/usr/lib/OpenBLAS/include/
#LIBS = -L/usr/lib/OpenBLAS/lib -lopenblas -lpthread -llapack -fopenmp
#INCLUDE = -I/opt/OpenBLAS/include/  
#LIBS = -L/opt/OpenBLAS/lib
#INCLUDE = -L/usr/lib/x86_x64-linux-gnu
LIBS = -fopenmp -llapack -lopenblas -larmadillo -larpack

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

clean:
	rm zigzag excitons test_file test_file_excitons

clean_spenctrum:
	rm spectrum*

clean_bands:
	rm bands*
