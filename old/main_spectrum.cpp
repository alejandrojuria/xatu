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
    //lambda = 0.0; // no SOC
    cx_mat eigvecX;
    vec eigvalX;

    // ----------------- Model parameters & Output --------------------
    int N = 15;

    std::string filename = "spectrum_bulk_energy_evolution_centeres_BZ";
    FILE* textfile = fopen(filename.c_str(), "w");
    bool writeEigvals = false;

    std::string filename_wf = "spectrum_bulk_wf_2bands";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
    bool printWF = false;

    bool calculateSpin = true;
    std::string filename_spin = "spectrum_bulk_spin_evolution_centeres_BZ";
    FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
    bool writeSpin = true;

    // ----------------------- Main body -----------------------

    initializeBlockMatrices();
    prepareHamiltonian(N);

    //double zeeCenter = 0.693268 - PI/a; // According to main_bands calculations
    int nBulkBands = 2;       // Bulk bands
    vec zeemanArray = {0.01, 0.008, 0.005, 0.003, 0.001, 0.0008, 0.0005, 0.0003, 0.0001};
<<<<<<< Updated upstream:old/main_spectrum.cpp
    //vec cellArray = arma::regspace(10, 150);
    vec cellArray = {150};
=======
    vec cellArray = {10, 20, 50, 75, 100, 125, 150, 175, 200};
    //vec cellArray = {100};
>>>>>>> Stashed changes:lib/main_spectrum.cpp
    int Ncell = 100;
    int nk = 2*Ncell;
    double epsk = 0.001;
    vec kpoints = arma::linspace(0.0, 2*PI/a, nk);
    kpoints = kpoints(arma::span(0, nk-2));
    nk -= 1;

    int Qspacing = nk/40;
    vec QIndexArray = {0};
    int nEdgeStates = 0; // Incorporate edge states

    fprintf(textfile, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");

    for(unsigned int n = 0; n < (int)cellArray.n_elem; n++){

        Ncell = cellArray(n);
        nk = 2*Ncell;
        double epsk = 0.001;
        vec kpoints = arma::linspace(0.0, 2*PI/a, nk);
        kpoints = kpoints(arma::span(0, nk-2));
        nk -= 1;
        int nEdgeStates = 0;

<<<<<<< Updated upstream:old/main_spectrum.cpp
        double Q = 0;
=======
        double Q = 0.0;
>>>>>>> Stashed changes:lib/main_spectrum.cpp

        //zeeman = zeemanArray(n); // Overwrite zeeman term value

        cout << "Zeeman: " << zeeman << " eV" << endl;
        cout << "N. cells: " << Ncell << endl;

        cout << "nk = " << nk << endl;
        //cout << "Q = " << Q << "(" << Qindex << "/" << nk/2 << ")" << endl;
        cout << "#bulk bands: " << nBulkBands << endl;

        mat states = createBasis(N, Q, kpoints, nBulkBands, nEdgeStates);
        BShamiltonian(N, Ncell, states, kpoints);
        arma::eig_sym(eigvalX, eigvecX, HBS);

        vec energies = computeEnergies(eigvecX.col(0), HBS, HK);
        cout << energies(0) << energies(1) << endl;
        printf("%12.8lf\t%12.8lf\n", eigvalX(0), eigvalX(1));
        printf("%12.8lf\t%12.8lf\n", eigvalX(2), eigvalX(3));

        if(writeEigvals){
            //writeEigenvaluesToFile(textfile, eigvalX, Q);
            fprintf(textfile, "%d\t%12.9lf\n", nk, eigvalX(0));
        };

        if(printWF){
            int state = 0;
            int nbands = sqrt(states.n_rows/nk);
            int nbands2 = nbands*nbands;
            for (int i = 0; i < nk; i++){
                double coef = 0;
                for(int nband = 0; nband < nbands2; nband++){
                    coef += abs(eigvecX.col(state)(nbands2*i + nband))*abs(eigvecX.col(state)(nbands2*i + nband));
                    coef += abs(eigvecX.col(state + 1)(nbands2*i + nband))*abs(eigvecX.col(state + 1)(nbands2*i + nband));
                };
                coef /= (kpoints(1) - kpoints(0)); // L2 norm instead of l2
                fprintf(textfile_wf, "%11.8lf\t%11.8lf\n", kpoints(i), coef);
            };
            fprintf(textfile_wf, "#\n");
        };

        cx_vec spin;
        if(calculateSpin){
                //cout << "Computing " << i+1 << " exciton spin..." << endl;
                spin = spinX(eigvecX.col(0), states, kpoints, N);
                cout << "Hole Sz is: " << spin(0) << endl;
                cout << "Electron Sz is: " << spin(1) << "\n" << endl;
                cout << "Total spin is: " << spin(2) << "\n" << endl;
            
            fprintf(textfile_spin, "%d\t%13.10lf\n", nk, spin(2));
        };
    };

    fclose(textfile);
    fclose(textfile_wf);
    fclose(textfile_spin);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};