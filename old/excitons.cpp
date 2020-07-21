#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "zigzag.hpp"
#include "excitons.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

// Shared variable definition
double a, c;
double Es, Ep, Vsss, Vsps, Vpps, Vppp;
double lambda;
vec a1, a2, tau;
vec n1, n2, n3;
vec Gamma, K, M;
mat M0, M1, M2p, M2m;
cx_mat Mso;
cx_mat H0, Ha, Hsoc;

cx_mat HBS;
mat HK;

auto start = high_resolution_clock::now();

/*      =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x ò 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= 
Source: http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mstvh0_cpp.txt */
void STVH0(double X, double *SH0) {
        double A0,BY0,P0,Q0,R,S,T,T2,TA0;
	int K, KM;

        S=1.0;
        R=1.0;
        if (X <= 20.0) {
           A0=2.0*X/PI;
           for (K=1; K<61; K++) {
              R=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
    e15:       *SH0=A0*S;
        }
        else {
           KM=int(0.5*(X+1.0));
           if (X >= 50.0) KM=25;
           for (K=1; K<=KM; K++) {
              R=-R*pow((2.0*K-1.0)/X,2);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
    e25:       T=4.0/X;
           T2=T*T;
           P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
           Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
           TA0=X-0.25*PI;
           BY0=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
           *SH0=2.0/(PI*X)*S+BY0;
        }
}

/* Obtain eigenstates from tight-binding hamiltonian. 
   Input: double k. Output: complex matrix 2*(N+1)*8 x 2*(N+1)*8*/
cx_mat eigenstatesH0(double k){
    vec eigenval;
    cx_mat eigenvec;
    cx_mat H = hamiltonian(k, H0, Ha, Hsoc);
    arma::eig_sym(eigenval, eigenvec, H);

    return eigenvec;
};

/* Calculate value of interaction potential (Keldysh). Units are eV.
   Input: double k. Output: complex double */
double potential(double r){
    double eps = 1.;
    double eps1 = 1.;
    double eps2 = 5.06;
    double eps_bar = (eps1 + eps2)/2;
    double r0 = c*eps/(eps1 + eps2);
    double R = abs(r)/r0;
    double SH0;

    if(r == 0){
        STVH0(a, &SH0);
        return ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(a));
    }
    else{
        STVH0(R, &SH0);
        return ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(R));
    };
};

/* Calculate lattice Fourier transform of Keldsyh potential
   Input: double k, int Ncell. Output:  complex double. Vk */
std::complex<double> fourierTrans(double k, int Ncell = 200){
    std::complex<double> i(0,1);
    std::complex<double> Vk = potential(0);
    for(int n = 0; n < Ncell; n++){
        Vk += (potential((n+1) * a))*std::exp(i * k * (double)(n+1)*a) + 
              (potential(- (n+1) * a))*std::exp(- i * k * (double)(n+1)*a);
    };
    return Vk;
};

std::complex<double> tDirect(std::complex<double> Vk,
                             const arma::cx_vec& coefsK, 
                             const arma::cx_vec& coefsKQ,
                             const arma::cx_vec& coefsK2, 
                             const arma::cx_vec& coefsK2Q, 
                             int Ncell){
    
    std::complex<double> D = 1./(2*Ncell + 1);
    D  *= Vk*arma::cdot(coefsKQ, coefsK2Q)*arma::cdot(coefsK2, coefsK);

    return D;
};

std::complex<double> tExchange(std::complex<double> VQ, 
                               const arma::cx_vec& coefsK, 
                               const arma::cx_vec& coefsKQ,
                               const arma::cx_vec& coefsK2, 
                               const arma::cx_vec& coefsK2Q, 
                               int Ncell){
    
    std::complex<double> X = 1./(2*Ncell + 1);
    X *= VQ*cdot(coefsKQ, coefsK)*cdot(coefsK2, coefsK2Q);

    return X;
};

/* Initilise basis to be used in the construction of the BSE matrix.
   We consider only the latest valence band, and the first conduction band.
   Input: int N (cells along finite direction), 
   int nEdgeStates (0 by default, otherwise odd). Output: matrix */
arma::mat createBasis(int N, double Q, const vec& kpoints, int nEdgeStates){

    int nk = kpoints.n_elem;

    // Sanity check
    if (nEdgeStates > nk){
        throw std::invalid_argument("More edge states than k points ");
    };

    vec valence = arma::regspace(2*(N+1)*5 - 3, 2*(N+1)*5 - 3);
    vec conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2);
    vec edgeC = arma::regspace(2*(N+1)*5 + 0, 2*(N+1)*5 + 0);
    vec edgeV = arma::regspace(2*(N+1)*5 - 0, 2*(N+1)*5 - 2);

    int nv = valence.n_elem;
    int nc = conduction.n_elem;

    int basisDimension = nv*nc*nk + nEdgeStates*edgeC.n_elem*edgeV.n_elem;
    mat states = arma::zeros(basisDimension, 4);

    int it = 0;
    for (int i = 0; i < nk; i++){

        if ((nEdgeStates > 0) && 
            (abs(i - (int)(nk/2)) <= (int)(nEdgeStates/2))){

            valence = arma::join_cols(valence, edgeV);
            conduction = arma::join_cols(conduction, edgeC);
        };

        for (int j = 0; j < (int)conduction.n_elem; j++){
            for (int k = 0; k < (int)valence.n_elem; k++){
                mat state = { conduction(j), valence(k), kpoints(i), Q };
                states.row(it) = state;

                it++;
            };
        };
    };

    return states;
};

/* Routine to calculate the index i associated to k within the 
kpoints vector.
Input: double k, vec kpoints. Output: index i */
int determineKIndex(double k, const vec& kpoints){
    int ndiv = kpoints.n_elem - 1;
    return round((k + PI/a)*a*ndiv/(2*PI));
};

/* Initialize BSE hamiltonian matrix and kinetic matrix. Recursive approach:
Instead of calculating the energies and coeficients dinamically, which
is too expensive, instead we first calculate those for each k, save them
in the stack, and then call them consecutively as we build the matrix.
Analogously, we calculate the Fourier transform of the potential beforehand,
saving it in the stack so that it can be later called in the matrix element
calculation.
Input: int N (cells finite direction), vec states, int Ncells (periodic 
direction), int nEdgeStates. Output: None (updates previously declared matrices) */
void BShamiltonian(int N, int Ncell, const mat& states, const vec& kpoints){

    int basisDimBSE = states.n_rows;
    int basisDimTB  = 2*(N+1)*8;
    int nk = kpoints.n_elem;

    HBS = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK  = arma::zeros(basisDimBSE, basisDimBSE);

    cx_cube eigvecKStack(basisDimTB, basisDimTB, nk);
    mat eigvalKStack(basisDimTB, nk);
    cx_vec  ftStack = arma::zeros<cx_vec>(nk, 1);

    double k, k2;
    int k_index, kQ_index, k2_index, k2Q_index;
    double Q = states(0,3);
    int Q_index = determineKIndex(Q, kpoints);
    std::complex<double> D, X;

    // ---------------- Load data in the stack ----------------
    vec auxEigVal(basisDimTB);
    for (int i = 0;  i < nk; i++){

        arma::eig_sym(auxEigVal, eigvecKStack.slice(i), 
                      hamiltonian(kpoints(i), H0, Ha, Hsoc));
        
        eigvalKStack.col(i) = auxEigVal;
        ftStack(i) = fourierTrans(kpoints(i), Ncell);
    };

    // ----------- Calculate matrix elements of HBS, HK -----------
    for (int i = 0; i < basisDimBSE; i++){

        k = states(i, 2);
        int c = states(i, 0);
        int v = states(i, 1);

        k_index = determineKIndex(k, kpoints);
        if (Q != 0){
            kQ_index = 
            (k_index + Q_index < nk) ? k_index + Q_index : k_index + Q_index - nk + 1;
        }
        else{
            kQ_index = k_index;
        }

        for (int j = i; j < basisDimBSE; j++){

            k2 = states(j, 2);
            int c2 = states(j, 0);
            int v2 = states(j, 1);

            k2_index = determineKIndex(k2, kpoints);
            if (Q != 0){
                k2Q_index = 
                (k2_index + Q_index < nk) ? k2_index + Q_index : k2_index + Q_index - nk + 1;
            }
            else{
                k2Q_index = k2_index;
            };

            cx_vec coefsK = eigvecKStack.slice(k_index).col(v);
            cx_vec coefsKQ = eigvecKStack.slice(kQ_index).col(c);
            cx_vec coefsK2 = eigvecKStack.slice(k2_index).col(v2);
            cx_vec coefsK2Q = eigvecKStack.slice(k2Q_index).col(c2);

            D = tDirect(ftStack(abs(k_index - k2_index), 0),
                             coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);
            X = tExchange(ftStack(Q_index, 0),
                             coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);

            if (i == j){
                HBS(i, j) = (eigvalKStack(states(i, 0), kQ_index) - eigvalKStack(states(i, 2), k_index))/2;
                HK(i, j) = HBS(i, j).real();
            }
            else{
                HBS(i, j) = eigvalKStack(states(i, 0), kQ_index) - eigvalKStack(states(i, 2), k_index) -
                           (D - X);
            };
            HBS = HBS + HBS.t();
        };
    };
};

/* Compute expected value of tight-binding energy and potential term for 
an exciton eigenstate.
Input: cx_vec eigvec, cx_mat HBS, cx_mat HK. Output: Vector 2x1 (<T>, <V>) */
vec computeEnergies(const cx_vec& eigvec, const cx_mat& HBS, const mat& HK){

    std::complex<double> kineticEnergy = arma::cdot(eigvec, HK*eigvec);

    cx_mat HV = HBS - HK;
    std::complex<double> potentialEnergy = arma::cdot(eigvec, HV*eigvec);

    vec energies = {kineticEnergy.real(), potentialEnergy.real()};
    return energies;
};

int main(){
    
    initializeConstants();
    cx_mat eigvecX;
    vec eigvalX;

    // -------------------- Model parameters --------------------
    int N = 15;
    int Ncell = 200;
    int nk = 25;
    vec kpoints = arma::linspace(-PI/a, PI/a, nk);
    double Q = 0.;

    std::string filename = "test_file_excitons";
    FILE* textfile = fopen(filename.c_str(), "w");

    // ----------------------- Main body -----------------------

    initializeBlockMatrices();
    prepareHamiltonian(N);
    mat states = createBasis(N, Q, kpoints);
    BShamiltonian(N, Ncell, states, kpoints);
    arma::eig_sym(eigvalX, eigvecX, HBS);

    vec energies = computeEnergies(eigvecX.col(0), HBS, HK);

    fprintf(textfile, "TB energy\tV Energy\n");
    fprintf(textfile, "%10.7f\%10.7f", energies(0), energies(1));
	fclose(textfile);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};