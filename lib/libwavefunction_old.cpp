#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "Zigzag.hpp"
#include "Exciton.hpp"
#include "libwavefunction.hpp"

/* Routine to calculate the coefficients of the exciton real-space 
wavefunction, each one associated to each atom of the lattice.
Here we use the vectorized implementation.
NB: this function requires having the basis predefined, meaning
that we need valence and conduction arrays.
Input: cx_mat BSE coefficiens, vec kpoints,
int holePos (1D representation of positions, i.e. index),
int Ncell, int N, int nEdgeStates.
Output: cx_vec real-space coefficients; already summed over k */
arma::cx_mat RScoefficientCalc(Exciton& exciton, const arma::cx_vec& BSEcoefs, 
                               int holePos){


    arma::vec eigValTB;
    arma::cx_mat eigVecTB;
    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = 2*(exciton.N+1)*8;
    int nk = exciton.kpoints.n_elem;
    // Matrix dimension is: #orbitals*#motif*#orbitals
    arma::cx_mat RScoefs = arma::zeros<arma::cx_mat>(dimTB*8, nk);

    for (int i = 0; i < nk; i++){

        kBlockEnd += exciton.nBulkBands*exciton.nBulkBands + 
                     exciton.nEdgeBands*exciton.nEdgeBands + 
                     2*exciton.nBulkBands*exciton.nEdgeBands;
        
        eigVecTB = exciton.eigvecKStack.slice(i);
        arma::cx_mat c = exciton.eigvecKQStack.slice(i);

        arma::cx_mat v = eigVecTB.rows(8*holePos, 8*(holePos+1) - 1);
        v = v.cols(0, exciton.nBulkBands + exciton.nEdgeBands - 1);
        c = c.cols(exciton.nBulkBands + exciton.nEdgeBands, c.n_cols - 1);

        arma::cx_mat M = arma::kron(c, arma::conj(v));
        arma::cx_vec A = BSEcoefs(arma::span(kBlock, kBlockEnd - 1));

        arma::cx_vec coefsK = M*A;
        RScoefs.col(i) = coefsK;

        kBlock = kBlockEnd;
    };

    return RScoefs;
};

/* Instead of using the s hydrogenic wave functions, we are simply
going to calculate the coefficient corresponding to each atom, since
initially we assumed orthonormality (i.e. they are punctual).
It is already modulus squared. 
Input: int n (1D atom position representation). 
int cell (corresponding unit cell) 
cx_mat RScoefs (matrix of real space coefficients as calculated previously)
vec kpoints.
Output: coefficient probability on atom (n, cell)*/
double atomCoefficientSquared(int n, int cell, int hCell, const arma::cx_mat& RScoefs,
                              Exciton& exciton){

    std::complex<double> coef = 0;
    double tCoef = 0;
    std::complex<double> j(0,1);
    for(int orb = 0; orb < 64; orb++){
        coef = 0;
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += RScoefs.col(i)(64*n + orb)*exp(j * exciton.kpoints(i) * (double)(cell - hCell)* exciton.a);
        };
        tCoef += abs(coef)*abs(coef);
    };

    return tCoef;
};