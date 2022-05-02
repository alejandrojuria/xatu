#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "GExciton.hpp"
#include "System.hpp"
#include "Crystal.hpp"
#include "Result.hpp"
#include "utils.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    std::cout << argc << std::endl;

    if (argc != 2){
        throw std::invalid_argument("Error: One input file is required (system config.)");
    };

    auto start = high_resolution_clock::now();
    
    int ncell = 60;
    double filling = 1./2;
    arma::rowvec Q = {0., 0., 0.};
    arma::rowvec parameters = {1., 1., 10};
    std::string modelfile = argv[1];    

    // ----------------- Model parameters & Output --------------------
    bool writeEigvals = false;
    std::string filename = "eigval.out";
    FILE* textfile_en = fopen(filename.c_str(), "w");

    bool writeStates = true;
    std::string filename_st = "states_tb.out";
    FILE* textfile_st = fopen(filename_st.c_str(), "w");

    // -------------------------- Main body ---------------------------

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "N. cells: " << ncell*ncell << endl;

    System model = System(modelfile);
	model.setFilling(filling);

	FILE* filecrystal = fopen("crystal.out", "w");
	arma::rowvec K = model.reciprocalLattice.row(0)*2./3. + model.reciprocalLattice.row(1)*1./3.;
	arma::rowvec rotatedK = model.rotateC3(K);
	arma::rowvec rotatedK2 = model.rotateC3(rotatedK);
	arma::rowvec M = model.reciprocalLattice.row(0)*1./2.;
	writeVectorToFile(K, filecrystal);
	writeVectorToFile(rotatedK, filecrystal);
	writeVectorToFile(rotatedK2, filecrystal);
	writeVectorToFile(M, filecrystal);
	fclose(filecrystal);
    
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    model.brillouinZoneMesh(ncell);

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

	arma::cx_mat h, eigvec;
	arma::vec eigval;
    for(int i = 0; i < model.nk; i++){
		h = model.hamiltonian(model.kpoints.row(i));
		arma::eig_sym(eigval, eigvec, h);

		if (writeEigvals){
			for (int j = 0; j < eigval.n_elem; j++){
				fprintf(textfile_en, "%f\t", eigval(j));
			}
		}

		if(writeStates){
			for (int j = 0; j < eigval.n_elem; j++){
				for(int n = 0; n < eigval.n_elem; n++){
					fprintf(textfile_st, "%f\t%f\t", real(eigvec.col(j)(n)), imag(eigvec.col(j)(n)));
				}
				fprintf(textfile_st, "\n");
			}
		}
	}
    
    fclose(textfile_en);
    fclose(textfile_st);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};