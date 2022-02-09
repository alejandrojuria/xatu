#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "System.hpp"

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

    std::string filename = "./models/2band_insulator_model.txt";
    System system = System(filename);
	arma::mat kpoints;

	int Ncell = 100;
	int nk = 50;
	if (system.ndim == 1){
		kpoints = system.brillouinZoneMesh(Ncell);
	}
	else if (system.ndim == 2){
		
		double norm = arma::norm(system.reciprocalLattice.row(0));
		arma::rowvec M = system.reciprocalLattice.row(0)/2.;
		arma::rowvec K = (system.reciprocalLattice.row(0) + 
						  system.reciprocalLattice.row(1))/2.;
		arma::rowvec G = {0., 0., 0.};
		kpoints = arma::zeros(3*nk + 1, 3);
		for(int i = 0; i < nk; i++){
			kpoints.row(i) = G*(nk-i)/nk + K*i/nk;
		}
		for(int i = 0; i < nk; i++){
			kpoints.row(i + nk) = K*(nk-i)/nk + M*i/nk;
		}
		for(int i = 0; i < nk; i++){
			kpoints.row(i + 2*nk) = M*(nk-i)/nk + G*i/nk;
		}
		kpoints.row(3*nk) = G;
		/*
		int nkdiv = 100;
		kpoints = arma::zeros(3*nkdiv + 1, 3);
		for(int i = 0; i < nkdiv; i++){
			double alpha = (double)i/(nkdiv+1);
			kpoints.row(i) = G*(1 - alpha) + M*alpha;
		}
		for(int i = 0; i < nkdiv; i++){
			double alpha = (double)i/(nkdiv+1);
			kpoints.row(i + nkdiv) = M*(1 - alpha) + K*alpha;
		}
		for(int i = 0; i < nkdiv; i++){
			double alpha = (double)i/(nkdiv+1);
			kpoints.row(i + nkdiv*2) = K*(1 - alpha) + G*alpha;
		}
		kpoints.row(3*nkdiv) = G;

		kpoints = arma::zeros(3, 3);
		kpoints.row(0) = K/5;
		kpoints.row(1) = system.rotateC3(K/5);
		kpoints.row(2) = system.rotateC3(system.rotateC3(K/5));
		*/

	};	

	double Q = 0.0;

	// ---------------------------- Main loop ----------------------------
	cx_vec stored_eigvec;
	for (int i = 0; i < kpoints.n_rows; i++) {

        arma::vec eigenval;
        arma::cx_mat eigenvec;
		arma::cx_mat h = system.hamiltonian(kpoints.row(i), true);
		arma::eig_sym(eigenval, eigenvec, h);
		arma::vec eigval = real(eigenval);
		eigval = arma::sort(eigval);

		if(writeBands){
            fprintf(textfile, "%d\t", i);
            for(int j = 0; j < system.basisdim; j++){
                fprintf(textfile, "%lf\t", eigval(j));
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