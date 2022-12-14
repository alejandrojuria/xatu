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
#include "CrystalDFTConfiguration.hpp"
#include "BiRibbon.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    auto start = high_resolution_clock::now();

    // ----------------- Model parameters & Output --------------------
    std::string filename = "eigval.out";
    FILE* textfile_en = fopen(filename.c_str(), "w");

    // -------------------------- Main body ---------------------------

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    int N = 19;
    BiRibbon model = BiRibbon(N);
    model.applyElectricField(0.003);
    
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    arma::rowvec G = {0., 0., 0.};
	arma::rowvec K = model.reciprocalLattice.row(0)/2.;

    int nk = 200;
    model.brillouinZoneMesh(nk);

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

	arma::cx_mat h, eigvec;
	arma::vec eigval;
    double left_edge_occupation, right_edge_occupation, spin;

    for(int i = 0; i < model.kpoints.n_rows; i++){
        arma::rowvec k = model.kpoints.row(i);
		model.solveBands(k, eigval, eigvec);

		for (int j = 0; j < eigval.n_elem; j++){
            left_edge_occupation = arma::norm(eigvec.col(j).subvec(0, 15));
            right_edge_occupation = arma::norm(eigvec.col(j).subvec(2*(N+1)*8 - 16, 2*(N+1)*8 - 1));
            spin = model.expectedSpinZValue(eigvec.col(j));
			fprintf(textfile_en, "%f\t%f\t%f\t%f\t", eigval(j), left_edge_occupation, right_edge_occupation, spin);
		}
        fprintf(textfile_en, "\n");

	}

    
    fclose(textfile_en);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};