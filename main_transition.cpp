/* Script to compute transition rate from ground state bulk exciton to
 edge non-interacting electron-hole pair according to Fermi Golden rule */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "Exciton.hpp"
#include "Zigzag.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){

    //omp_set_num_threads(24);
    
    // lambda = 0.0; // no SOC
    cx_mat eigvecX;
    vec eigvalX;

    std::string filename = "transition_ribbon_4_300_noonsite";
    FILE* textfile = fopen(filename.c_str(), "a");

    int N = 15;
    int nBulkBands = 4;
    int nEdgeBands = 1;
    double Q = 0.0;

    // ----------------------- Main body -----------------------

    vec NArray = arma::regspace(7, 4, 31);
    for(int n = 0; n < (int)NArray.n_elem; n++){

        auto start = high_resolution_clock::now();

        int N = NArray(n); 
        int Ncell = 300;
        int nk = 2*Ncell - 1;

        cout << "Ribbon width: " << N << " cells" << endl;
        cout << "N. cells: " << Ncell << endl;
        cout << "nk = " << nk << endl;
        cout << "Q = " << Q << endl;
        cout << "#bulk bands: " << nBulkBands << endl;
        cout << "Computing GS bulk exciton..." << endl;

        Exciton bulkExciton = Exciton(N, Ncell, Q, nBulkBands, 0);
        bulkExciton.BShamiltonian();
        cout << "Diagonalizing... " << std::flush;
        int eig_number = 10;
        //eigs_gen(eigvalX, eigvecX, bulkExciton.HBS, eig_number, "sr");
        eig_sym(eigvalX, eigvecX, bulkExciton.HBS);
        cout << "Done" << endl;

        vec energies = bulkExciton.computeEnergies(eigvecX.col(0));
        cout << "Eigenenergy: " << eigvalX(0) << " eV" << endl;
        cout << "Kinetic energy: " << energies(0) << "eV" << endl;
        cout << "Potential energy: " << energies(1) << "eV" << endl;

        double initialEnergy = eigvalX(0);
        cx_vec initialCoefs = eigvecX.col(0);
        initialCoefs /= sqrt(bulkExciton.kpoints(1) - bulkExciton.kpoints(0));

        cout << "Computing transition rate..." << endl;

        Exciton exciton = Exciton(N, Ncell, Q, nBulkBands, nEdgeBands, vec{0, 1});
        double transitionRate = exciton.fermiGoldenRule(initialCoefs, initialEnergy);

        cout << "Done" << endl;
        cout << "Transition rate is: " << transitionRate << "\n" << endl;

        fprintf(textfile, "%d\t%.10e\n", N, transitionRate);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
	    std::cout << "Elapsed time during iteration: " << duration.count() << " ms" << std::endl;
    };

    fclose(textfile);

	return 0;
};
