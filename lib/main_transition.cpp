/* Script to compute transition rate from ground state bulk exciton to
 edge non-interacting electron-hole pair according to 
 Fermi Golden rule */

#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "excitons.hpp"
#include "zigzag.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

// Shared variable declaration
double a, c;
double Es, Ep, Vsss, Vsps, Vpps, Vppp;
double lambda, zeeman, onsiteEdge;
vec a1, a2, tau;
vec n1, n2, n3;
vec Gamma, K, M;
mat M0, M1, M2p, M2m, Mzeeman;
cx_mat Mso;
cx_mat H0, Ha, Hsoc, Hzeeman;

cx_mat HBS;
mat HK;
mat states;

int main(){

    auto start = high_resolution_clock::now();
    
    initializeConstants();
    // lambda = 0.0; // no SOC
    cx_mat eigvecX;
    vec eigvalX;

    std::string filename = "transition_dumb";
    FILE* textfile = fopen(filename.c_str(), "w");

    int N = 15;

    // ----------------------- Main body -----------------------

    initializeBlockMatrices();
    prepareHamiltonian(N);

    int nBulkBands = 2; // Bulk bands
    vec NcellArray = {100, 200, 300, 400, 500, 600};
    for(int n = 0; n < (int)NcellArray.n_elem; n++){

        int Ncell = NcellArray(n);

        int nk = 2*Ncell;
        double epsk = 0.001;
        vec kpoints = arma::linspace(0.0, 2*PI/a, nk);
        kpoints = kpoints(arma::span(0, nk-2));
        nk -= 1;
        int Qindex = 0;
        double Q = kpoints(Qindex);

        cout << "N. cells: " << Ncell << endl;
        cout << "nk = " << nk << endl;
        cout << "Q = " << Q << "(" << Qindex << "/" << nk/2 << ")" << endl;
        cout << "#bulk bands: " << nBulkBands << endl;
        cout << "Computing GS bulk exciton...\n" << endl;

        mat states = createBasis(N, Q, kpoints, nBulkBands, 0);

        BShamiltonian(N, Ncell, states, kpoints);
        arma::eig_sym(eigvalX, eigvecX, HBS);
        cout << "Done" << endl;

        vec energies = computeEnergies(eigvecX.col(0), HBS, HK);
        cout << "Eigenenergy: " << eigvalX(0) << " eV" << endl;
        cout << "Kinetic energy: " << energies(0) << "eV" << endl;
        cout << "Potential energy: " << energies(1) << "eV" << endl;

        double initialEnergy = eigvalX(0);
        cx_vec initialCoefs = eigvecX.col(0);
        initialCoefs /= sqrt(kpoints(1) - kpoints(0));

        cout << "Computing transition rate...\n" << endl;
        
        double transitionRate = fermiGoldenRule(initialCoefs, initialEnergy,
                nBulkBands, Q, N, Ncell, kpoints);

        cout << "Done" << endl;
        cout << "Transition rate is: " << transitionRate << endl;

        fprintf(textfile, "%d\t%.10e\n", nk, transitionRate);
    }

    fclose(textfile);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};