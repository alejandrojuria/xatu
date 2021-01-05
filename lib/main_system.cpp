#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "System.hpp"
#include "GExciton.hpp"


#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){

    auto start = high_resolution_clock::now();

	// ------------------------ Model output into files ------------------------

	FILE* textfile = fopen("bands_from_file", "w");
	bool writeBands = true;

	// ------------------------ Initialization ------------------------

    std::string filename = "model_bi";
    System system = System(filename);

	int Ncell = 100;
    arma::mat kpoints = system.brillouin_zone_mesh(Ncell);
	double Q = 0.0;

	// ---------------------------- Main loop ----------------------------
	for (int i = 0; i < Ncell; i++) {

        arma::vec eigenval;
        arma::cx_mat eigenvec;
		arma::cx_mat h = system.hamiltonian(kpoints.row(i));
		arma::eig_sym(eigenval, eigenvec, h);

		if(writeBands){
            fprintf(textfile, "%lf\t", kpoints.row(i)(1));
            for(int j = 0; j < system.basisdim; j++){
                fprintf(textfile, "%lf\t", eigenval(j));
            };
            fprintf(textfile, "\n");
		};
	};

	// Close all text files
	fclose(textfile);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};