#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "System.hpp"
#include "GExciton.hpp"
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
arma::cx_mat RScoefficientCalc(GExciton& exciton, const arma::cx_vec& BSEcoefs, 
                               int holePos){


    arma::vec eigValTB;
    arma::cx_mat eigVecTB;
    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = exciton.basisdim;
    int nk = exciton.nk;
    // Matrix dimension is: #orbitals*#motif*#orbitals
    arma::cx_mat RScoefs = arma::zeros<arma::cx_mat>(dimTB*exciton.norbitals, nk);

    for (int i = 0; i < nk; i++){

        kBlockEnd += (exciton.excitonbasisdim/nk);
        
        eigVecTB = exciton.eigvecKStack.slice(i);
        arma::cx_mat c = exciton.eigvecKQStack.slice(i);

        /* Precaution when picking the coeffs. associated to holePos, since
        now our basis splits H in two blocks for spin*/
        arma::cx_mat v = eigVecTB.rows(exciton.norbitals*holePos, 
                                        exciton.norbitals*(holePos+1) - 1);
        /*arma::cx_mat v2 = eigVecTB.rows(dimTB/2 + exciton.norbitals/2*holePos, 
                                        dimTB/2 + exciton.norbitals/2*(holePos+1) - 1);
        arma::cx_mat v = arma::join_cols(v1, v2);*/

        v = v.cols(0, exciton.valenceBands.n_elem - 1);
        c = c.cols(exciton.valenceBands.n_elem, c.n_cols - 1);

        //arma::cx_mat M = arma::kron(c, arma::conj(v));
        arma::cx_mat M = arma::kron(c, v);
        arma::cx_vec A = BSEcoefs(arma::span(kBlock, kBlockEnd - 1));
        //arma::cout << "M: " << M << arma::endl;
        //arma::cout << "A: " << A << arma::endl;

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
double atomCoefficientSquared(int n, const arma::rowvec& cell, const arma::rowvec& hCell, 
                              const arma::cx_mat& RScoefs, GExciton& exciton){

    std::complex<double> coef = 0;
    double tCoef = 0;
    std::complex<double> imag(0,1);
    int orbitalblock = exciton.norbitals*exciton.norbitals;

    // We need to separate in spin blocks
    for(int orb = 0; orb < orbitalblock; orb++){        

        coef = 0;
        // Spin up
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += RScoefs.col(i)(orbitalblock*n + orb)*exp(imag * arma::dot(exciton.kpoints.row(i), cell -  
            hCell));
        };
        tCoef += abs(coef)*abs(coef);
         
        // Spin down
        /*coef = 0;
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += (RScoefs.col(i)(RScoefs.n_rows/2 + orbitalblock*n + orb)*exp(imag * arma::dot(exciton.kpoints.row(i), cell - hCell)));
        };
        tCoef += abs(coef)*abs(coef); */
    };
    
    //tCoef = real(coef);
    // arma::cout << real(coef) << "---" << tCoef << arma::endl;

    return tCoef;
};

double fourierTransformExciton(const arma::cx_vec& BSEcoefs, 
                               GExciton& exciton, 
                               const arma::rowvec& electron_position, const arma::rowvec& hole_position){

    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = exciton.basisdim;
    int nk = exciton.nk;
    arma::cx_double ft;
    std::complex<double> imag(0,1);
    // Matrix dimension is: #orbitals*#motif*#orbitals

    for (int i = 0; i < nk; i++){

        kBlockEnd += (exciton.excitonbasisdim/nk);
        
        arma::cx_vec A = BSEcoefs(arma::span(kBlock, kBlockEnd - 1));
        //arma::cx_vec v = exciton.eigvecKStack.slice(i).col();
        //arma::cx_vec c = exciton.eigvecKStack.slice(i).col();
        double summed_coefs = pow(arma::norm(A), 2);

        ft += summed_coefs * std::exp(-imag*arma::dot(exciton.kpoints.row(i), electron_position - hole_position));

        kBlock = kBlockEnd;
    };
    double result = std::real(ft);

    return result;
}

double realSpaceWavefunction(GExciton& exciton, const arma::cx_vec& BSEcoefs,
                             int electronIndex, int holeIndex,
                             const arma::rowvec& eCell, const arma::rowvec& hCell){

    std::complex<double> imag(0, 1);
    std::complex<double> amplitude = 0;
    std::complex<double> totalAmplitude = 0;
    arma::cx_vec BSEsquared = arma::conj(BSEcoefs) % BSEcoefs;
    // std::cout << arma::max(BSEsquared) << std::endl;
    double maxval = std::real(arma::max(BSEsquared));

    /* for (auto i : BSEsquared){
        if (std::real(i) > maxval/10){
            //std::cout << i << std::endl;
        }
    }*/

    for (int i = 0; i < exciton.excitonbasisdim; i++){
        for (int j = 0; j < exciton.excitonbasisdim; j++){


            amplitude = 0;
            arma::irowvec basisState = exciton.basisStates.row(i);
            arma::irowvec basisState2 = exciton.basisStates.row(j);

            int v1 = exciton.bandToIndex[basisState(0)];
            int v2 = exciton.bandToIndex[basisState2(0)];
            int c1 = exciton.bandToIndex[basisState(1)];
            int c2 = exciton.bandToIndex[basisState2(1)];
            int kIndex = basisState(2);
            int k2Index = basisState2(2);
            arma::rowvec kpoint = exciton.kpoints.row(kIndex);
            arma::rowvec kpoint2 = exciton.kpoints.row(k2Index);

            amplitude += std::exp(imag*arma::dot(kpoint - kpoint2, eCell - hCell))*
                         std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(electronIndex))*
                         exciton.eigvecKStack.slice(kIndex).col(v1)(holeIndex)*
                         std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(holeIndex))*
                         exciton.eigvecKQStack.slice(k2Index).col(c2)(electronIndex);

            if (false){
                amplitude += std::exp(imag*arma::dot(exciton.Q, eCell - hCell))*
                         std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(electronIndex))*
                         exciton.eigvecKQStack.slice(k2Index).col(c2)(holeIndex)*
                         std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(holeIndex))*
                         exciton.eigvecKStack.slice(kIndex).col(v1)(electronIndex);
            
                amplitude += std::exp(-imag*arma::dot(exciton.Q, eCell - hCell))*
                            std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(electronIndex))*
                            exciton.eigvecKStack.slice(kIndex).col(v1)(holeIndex)*
                            std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(holeIndex))*
                            exciton.eigvecKQStack.slice(k2Index).col(c2)(electronIndex);
                
                amplitude += std::exp(-imag*arma::dot(kpoint - kpoint2, eCell - hCell))*
                            std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(electronIndex))*
                            exciton.eigvecKQStack.slice(k2Index).col(c2)(holeIndex)*
                            std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(holeIndex))*
                            exciton.eigvecKStack.slice(kIndex).col(v1)(electronIndex);
            }

            amplitude *= std::conj(BSEcoefs(i))*BSEcoefs(j);
            totalAmplitude += amplitude;     
        }
    }

    //arma::cout << totalAmplitude << "---" << std::abs(totalAmplitude) << arma::endl;
    if (std::imag(totalAmplitude) > 1E-5){
        arma::cout << "Warning: Imaginary part of amplitude is not" << arma::endl;
    }
    return std::abs(totalAmplitude);
};


// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
std::complex<double> densityMatrix(GExciton& exciton, const arma::cx_vec& BSEcoefs, 
                                    int eIndex, int hIndex){

    std::complex<double> rho;
    int nk = exciton.nk;
    int nc = exciton.conductionBands.n_elem;
    int nv = exciton.valenceBands.n_elem;

    for (int i = 0; i < exciton.basisStates.n_rows; i++){
        arma::irowvec state = exciton.basisStates.row(i);
        int v = state(0);
        int c = state(1);
        int vI = i%nv;
        int cI = (int)(i/nv)%nc;
        int kIndex = state(2);
        arma::rowvec k = exciton.kpoints.row(kIndex);

        // Conduction term
        for (int cIndex = 0; cIndex < nc; cIndex++){
            int c2 = exciton.conductionBands(cIndex);
            rho += std::conj(BSEcoefs(vI + cIndex*nv + kIndex*nc*nv))*BSEcoefs(i)*
                    std::conj(exciton.eigvecKQStack.slice(kIndex)(eIndex, c2))*
                    exciton.eigvecKQStack.slice(kIndex)(hIndex, c);
        }
        // Valence term
        for (int vIndex = 0; vIndex < nv; vIndex++){
            int v2 = exciton.valenceBands(vIndex);
            rho -= std::conj(BSEcoefs(vIndex + cI*nv + kIndex*nc*nv))*BSEcoefs(i)*
                    std::conj(exciton.eigvecKStack.slice(kIndex)(eIndex, v))*
                    exciton.eigvecKStack.slice(kIndex)(hIndex, v2);
        }
    }
    
    // Fermi sea term;
    for(int i = 0; i < exciton.kpoints.n_rows; i++){
        arma::cx_rowvec blochState = exciton.eigvecKStack.slice(i).row(eIndex);
        blochState = blochState(arma::span(0, exciton.fermiLevel));
        arma::cx_rowvec blochState2 = exciton.eigvecKStack.slice(i).row(hIndex);
        blochState2 = blochState2(arma::span(0, exciton.fermiLevel));

        rho += arma::cdot(blochState, blochState2);
    }

    return rho;
};

// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
std::complex<double> densityMatrixK(int kIndex, GExciton& exciton, const arma::cx_vec& BSEcoefs, 
                                    int eIndex, int hIndex){

    std::complex<double> rho;
    int nk = exciton.nk;
    int nc = exciton.conductionBands.n_elem;
    int nv = exciton.valenceBands.n_elem;

    for (int cIndex = 0; cIndex < nc; cIndex++){
        for (int vIndex = 0; vIndex < nv; vIndex++){
            int v = exciton.valenceBands(vIndex);
            int c = exciton.conductionBands(cIndex);
            arma::rowvec k = exciton.kpoints.row(kIndex);

            // Conduction term
            for (int cIndex2 = 0; cIndex2 < nc; cIndex2++){
                int c2 = exciton.conductionBands(cIndex2);
                rho += std::conj(BSEcoefs(vIndex + cIndex2*nv + kIndex*nc*nv))*
                        BSEcoefs(vIndex + cIndex*nv + kIndex*nc*nv)*
                        std::conj(exciton.eigvecKStack.slice(kIndex)(eIndex, c2))*
                        exciton.eigvecKStack.slice(kIndex)(hIndex, c);
            }
            // Valence term
            for (int vIndex2 = 0; vIndex2 < nv; vIndex2++){
                int v2 = exciton.valenceBands(vIndex2);
                rho -= std::conj(BSEcoefs(vIndex2 + cIndex*nv + kIndex*nc*nv))*
                        BSEcoefs(vIndex + cIndex*nv + kIndex*nc*nv)*
                        std::conj(exciton.eigvecKStack.slice(kIndex)(eIndex, v))*
                        exciton.eigvecKStack.slice(kIndex)(hIndex, v2);
            }
        }
    }

    // Fermi sea term;
    arma::cx_rowvec blochState = exciton.eigvecKStack.slice(kIndex).row(eIndex);
    blochState = blochState(arma::span(0, exciton.fermiLevel));
    arma::cx_rowvec blochState2 = exciton.eigvecKStack.slice(kIndex).row(hIndex);
    blochState2 = blochState2(arma::span(0, exciton.fermiLevel));

    rho += arma::cdot(blochState, blochState2);
    return rho;
};