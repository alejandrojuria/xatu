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

    // ----------------- Model parameters & Output --------------------
    int N = 15;
    double Q = 0.;

    std::string filename = "convergence_edges_vnot0_noExtremes";
    FILE* textfile = fopen(filename.c_str(), "w");

    std::string filename_wf = "first_wf_coefs_1001nk_2bands";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
    bool printWF = false;

    // ----------------------- Main body -----------------------
    initializeBlockMatrices();
    prepareHamiltonian(N);

    vec NcellArray = {150, 300, 450, 600, 750, 900, 1000};
    vec kArray = NcellArray*2 + 1;
    fprintf(textfile, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");

    for(unsigned int n = 0; n < kArray.n_elem; n++){      

        int Ncell = NcellArray(n);
        int nkAux = kArray(n);
        int nbands = 0;
        cout << "nk = " << nkAux << endl;
        vec kpointsAux = arma::linspace(-PI/a, PI/a, nkAux);
        vec kpoints = kpointsAux(arma::span(0, kpointsAux.n_elem-1));
        int nk = kpoints.n_elem;
        int nEdgeStates = nk;

        mat states = createBasis(N, Q, kpoints, nbands, nEdgeStates);

        BShamiltonian(N, Ncell, states, kpoints, true);
        arma::eig_sym(eigvalX, eigvecX, HBS);

        vec energies = computeEnergies(eigvecX.col(0), HBS, HK);
        cout << eigvalX(0) << energies(0) << energies(1) << endl;

        fprintf(textfile, 
                "%d\t%10.7lf\t%10.7lf\t%10.7lf\n", 
                nk, eigvalX(0), energies(0), energies(1));

        if(printWF == true && nk == 1001){
            for (int i = 0; i < nk; i++){
                double coef = abs(eigvecX.col(0)(4*i))*abs(eigvecX.col(0)(4*i)) + 
                abs(eigvecX.col(0)(4*i+1))*abs(eigvecX.col(0)(4*i+1)) + 
                abs(eigvecX.col(0)(4*i+2))*abs(eigvecX.col(0)(4*i+2)) + 
                abs(eigvecX.col(0)(4*i+3))*abs(eigvecX.col(0)(4*i+3));

                fprintf(textfile_wf, 
                        "%11.8lf\t%11.8lf\n", 
                        kpoints(i), coef);
            };
            fclose(textfile_wf);
        };

    };

    fclose(textfile);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};