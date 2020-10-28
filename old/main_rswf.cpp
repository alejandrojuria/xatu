#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>
#include <omp.h>

#include "excitons.hpp"
#include "zigzag.hpp"
#include "libwavefunction.hpp"

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
double lambda, zeeman;
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

    // ---------------- Model parameters & Output ----------------
    int N = 15;
    double Q = 0.;

    std::string filename_en = "energies_placeholder";
    FILE* textfile_en = fopen(filename_en.c_str(), "w");
    bool print_en = false;

    std::string filename_k_wf = "first_wf_coefs_1001nk_2bands";
    FILE* textfile_k_wf = fopen(filename_k_wf.c_str(), "w");
    bool print_kWF = false;

    std::string filename_rs_wf = "rs_wf_first_301nk_4bands_noedge";
    FILE* textfile_rs_wf = fopen(filename_rs_wf.c_str(), "w");
    bool print_rsWF = true;

    // ----------------------- Main body -----------------------
    initializeBlockMatrices();
    prepareHamiltonian(N);    

    int Ncell = 50;
    int nk = 2*Ncell + 1;
    int nbands = 2;
    int nEdgeStates = nk; // full edge bands
    cout << "nk = " << nk << endl;
    cout << "nbands = " << nbands << endl;
    vec kpoints = arma::linspace(-PI/a, PI/a, nk);

    states = createBasis(N, Q, kpoints, nbands, nEdgeStates);
    BShamiltonian(N, Ncell, states, kpoints, true);
    arma::eig_sym(eigvalX, eigvecX, HBS);
    cout << "BSE diagonalized" << endl;

    cx_vec groundStateExciton = eigvecX.col(0);
    cx_vec degenerateGSExciton = eigvecX.col(1);
    cx_vec totalExciton = groundStateExciton;

    vec energies = computeEnergies(eigvecX.col(0), HBS, HK);
    cout << eigvalX(0) << eigvalX(1) << eigvalX(2) << eigvalX(3)
         << energies(0) << energies(1) << endl;

    // Print ground state kinetic, potential and total energies
    if(print_en == true){
        fprintf(textfile_en, "#bands\t#k points\tTotal energy\tTB energy\tV energy (eV)\n");
        fprintf(textfile_en, 
                "%d\t%d\t%10.7lf\t%10.7lf\t%10.7lf\n", 
                nbands, nk, eigvalX(0), energies(0), energies(1));
    };

    // Print k wavefunction
    if(print_kWF == true){
        for (int i = 0; i < nk; i++){
            double coef = abs(totalExciton(4*i))*abs(totalExciton(4*i)) + 
            abs(totalExciton(4*i+1))*abs(totalExciton(4*i+1)) + 
            abs(totalExciton(4*i+2))*abs(totalExciton(4*i+2)) + 
            abs(totalExciton(4*i+3))*abs(totalExciton(4*i+3));

            fprintf(textfile_k_wf, 
                    "%11.8lf\t%11.8lf\n",
                    kpoints(i), coef);
        };
    };

    // Print real-space wavefunction
    if(print_rsWF == true){

        cout << __LINE__ << endl;

        mat motiv = createMotiv(N);
        cout << __LINE__ << endl;

        int holeIndex = N; // Middle of nanoribbon
        int holeCell = 0; 
        rowvec holePosition = motiv.row(holeIndex) + 
                              rowvec({0, holeCell*a, 0});
                cout << __LINE__ << endl;


        cx_mat RScoefsGS = RScoefficientCalc(groundStateExciton, kpoints,
                            holeIndex, N, nEdgeStates);
        cout << __LINE__ << endl;

        cx_mat RScoefsDEG = RScoefficientCalc(degenerateGSExciton, kpoints,
                holeIndex, N, nEdgeStates);
        // Write hole position
        fprintf(textfile_rs_wf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                        holePosition(0), holePosition(1), 0.0);

        double coefSum = 0;
        vec coefs = arma::zeros(nk*2*(N+1));
        int it = 0;

        for(int cell = -Ncell; cell <= Ncell; cell++){
            for (int n = 0; n < 2*(N+1); n++){

                coefs(it) = 
                atomCoefficientSquared(n, cell, holeCell, RScoefsGS, kpoints);
                        
                coefSum += coefs(it);
                it++;
            };
        };
        cout << __LINE__ << endl;

        cout << coefSum << endl;
        it = 0;
        cout << "Writing w.f. to file" << endl;
        for(int cell = -Ncell; cell <= Ncell; cell++){
            for (int n = 0; n < 2*(N+1); n++){

                rowvec position = motiv.row(n) + rowvec({0, cell*a, 0});

                fprintf(textfile_rs_wf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                        position(0), position(1), coefs(it)/coefSum);
                it++;
            };
        };
    };

    fclose(textfile_en);
    fclose(textfile_k_wf);
    fclose(textfile_rs_wf);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};