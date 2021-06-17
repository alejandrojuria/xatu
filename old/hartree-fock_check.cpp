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

    initializeConstants();
    int N = 15;
    double Q = 0.;

    initializeBlockMatrices();
    prepareHamiltonian(N);

    vec NcellArray = {150, 300, 450, 600, 750, 900, 1050};
    vec nkArray = 2*NcellArray + 1;

    vec eigval;
    cx_mat eigvec;

    double valence = 2.*(N+1)*5-1;
    double conduction = 2.*(N+1)*5;

    cout << "nk\tEc\t\tEv\t\tD\t\tX\t\tTotal\n" << endl;

    for(int n = 0; n < (int)NcellArray.n_elem; n++){

        int Ncell = NcellArray(n);
        int nk = nkArray(n);

        int kindex = (int)(79*nk/80);
        vec kpoints = arma::linspace(-PI/a, PI/a, nk);

        vec state = {valence, conduction, kpoints(kindex), Q};

        arma::eig_sym(eigval, eigvec, hamiltonian(kpoints(kindex), H0, Ha, Hsoc));

        double eigc = eigval((int)conduction);
        double eigv = eigval((int)valence);
        arma::cx_vec coefC = eigvec.col((int)conduction);
        arma::cx_vec coefV = eigvec.col((int)valence);

        std::complex<double> fourierT = fourierTrans(0, Ncell);

        double D = real(tDirect(fourierT, coefV, coefC, coefV, coefC, Ncell));
        double X = real(tExchange(fourierT, coefV, coefC, coefV, coefC, Ncell));

        cout << nk << "\t" << eigc << "\t" << eigv << "\t" << D << "\t" << X << "\t" << (eigc-eigv-D+X) << endl;

    };


    return 0;
}