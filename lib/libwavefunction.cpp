#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "Zigzag.hpp"
#include "Exciton.hpp"
#include "libwavefunction.hpp"

/* Implementation of the hydrogenic atom wavefunction for the state
|100> (exponential, radial symmetry) to represent the exciton wavefunction
at each atom. 
We pass a 3 dimensional position array, although effectively we only use
the coordinates in-plane (x,y).
Input: array[3] position, array[3] reference position.
Output: double wavefunction value. */
double sWaveFunction(const arma::vec& position, const arma::vec& position0){

    int Z = 1;
    double cte = 2*sqrt(Z*Z*Z)/(sqrt(4*PI));
    // double a0 = 0.529; // Bohr radius [A]

    double x2 = (position(0) - position0(0))*(position(0) - position0(0));
    double y2 = (position(1) - position0(1))*(position(1) - position0(1));
    double r = sqrt(x2 + y2);

    double wf = cte*exp(-Z*0.5*r);
    return wf;
};

/* Routine to collapse the tight-binding coefficient matrix (eigenvalues)
by summing over rows corresponding to the same atom, but different orbital,
i.e., we treat all orbitals in the same footing for simplicity when 
plotting the exciton w.f.
Input: mat tight-binding coefficients.
Output: mat (dimension previous/#orbitals) collapsed TB coefs. */
arma::cx_mat collapseTBcoefs(const arma::cx_mat& coefs){
    
    int dimTB = coefs.n_rows;
    arma::cx_mat reducedCoefs = arma::zeros<arma::cx_mat>(dimTB/8, dimTB); 
    for(int n = 0; n < dimTB; n++){
        int i = (int)n/8;
        reducedCoefs.row(i) += coefs.row(n);
        arma::cout << coefs.row(n)(0) << arma::endl;
        arma::cout << reducedCoefs.row(i)(0) << arma::endl;
        arma::cout << "--------------" << arma::endl;

    };
    
    return reducedCoefs;
};

/* Routine to calculate the coefficients of the exciton real-space 
wavefunction, each one associated to each atom of the lattice.
Here we use the vectorized implementation.
NB: this function requires having the basis predefined, meaning
that we need valence and conduction arrays.
Input: cx_mat BSE coefficiens, vec kpoints,
int holePos (1D representation of positions, i.e. index),
int Ncell, int N, int nEdgeStates.
Output: cx_vec real-space coefficients; already summed over k */
arma::cx_mat RScoefficientCalc(const arma::cx_vec& BSEcoefs, 
                               const arma::vec& kpoints, 
                               int holePos, int N, 
                               int nEdgeStates){


    /* Valence, conduction and edges were in the basis construction
    scope so we have to define them again here. 
    CAREFUL: 
    They must be defined the same as in the BSE basis routine */

    arma::uvec valence = arma::regspace<arma::uvec>(2*(N+1)*5 - 4, 2*(N+1)*5 - 3);
    arma::uvec conduction = arma::regspace<arma::uvec>(2*(N+1)*5 + 2, 2*(N+1)*5 + 3);
    arma::uvec edgeC = arma::regspace<arma::uvec>(2*(N+1)*5 + 0, 2*(N+1)*5 + 1);
    arma::uvec edgeV = arma::regspace<arma::uvec>(2*(N+1)*5 - 2, 2*(N+1)*5 - 1);

    arma::uvec tValence;
    arma::uvec tConduction;
    
    int nv = valence.n_elem;
    int nc = conduction.n_elem;
    int nk = kpoints.n_elem;
    int basisDimension = nv*nc*nk + nEdgeStates*edgeC.n_elem*edgeV.n_elem + 
                         nEdgeStates*edgeC.n_elem*nv + 
                         nEdgeStates*edgeV.n_elem*nc;

    // Sanity check
    if(basisDimension != (int)states.n_rows){
        std::invalid_argument("Valence and conduction bands do not match basis states for BSE calculation");
    };

    arma::vec eigValTB;
    arma::cx_mat eigVecTB;
    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = 2*(N+1)*8;
    // Matrix dimension is: #orbitals*#motif*#orbitals
    arma::cx_mat RScoefs = arma::zeros<arma::cx_mat>(dimTB*8, nk);

    for (int i = 0; i < nk; i++){

        arma::eig_sym(eigValTB, eigVecTB, 
                      hamiltonian(kpoints(i), H0, Ha, Hsoc));
                
        if ((nEdgeStates > 0) && 
            (abs(i - (int)(nk/2)) <= (int)(nEdgeStates/2))){

            tValence = arma::join_cols(valence, edgeV);
            tConduction = arma::join_cols(conduction, edgeC);

            kBlockEnd += nv*nc + edgeC.n_elem*edgeV.n_elem + 
                         edgeC.n_elem*nv + edgeV.n_elem*nc;
        }
        else{
            tValence = valence;
            tConduction = conduction;

            kBlockEnd += nv*nc;
        };

        arma::cx_mat v = eigVecTB.rows(8*holePos, 8*(holePos+1) - 1);
        v = conj(v.cols(tValence)); // .t() takes hermitian conjugate
        arma::cx_mat c = eigVecTB.cols(tConduction);

        arma::cx_mat M = arma::kron(c, v);
        arma::cx_vec A = BSEcoefs(arma::span(kBlock, kBlockEnd - 1));

        arma::cx_vec coefsK = M*A;
        RScoefs.col(i) = coefsK;

        kBlock = kBlockEnd;
    };

    return RScoefs;
};

/* Routine to calculate the probability of finding the electron (hole) 
at a certain position. We can specify how many unit cells we want to use
in the calculation (we only specify upwards, although it is for both 
downwards and upwards), since those far away will be negligible.
Input: vec position, cx_vec RS coefs, int N, 
int Ncell desired (in one direction), vec kpoints.
Ouput: double electronicDensity.*/
double electronicDensity(const arma::vec& position, 
                         const arma::cx_mat& RScoefs, 
                         int N, int Ncell, const arma::vec& kpoints){
    
    std::complex<double> prob = 0;
    std::complex<double> i(0,1);
    arma::mat motiv = createMotiv(N);

    int cell = (int)position(1)/a;

    for(int j = 0; j < (int)RScoefs.n_cols; j++){
        arma::cx_vec coefs = RScoefs.col(j);

        for (int n = 0; n < (int)coefs.n_elem; n++){

            arma::vec atomPosition = motiv.row(n).t();
            arma::vec atomPositionUp = atomPosition + 
                                       arma::vec({0, cell*a, 0});
            arma::vec atomPositionDown = atomPosition + 
                                         arma::vec({0, -cell*a, 0});

            prob += coefs(n)*sWaveFunction(position, atomPosition)*
                    exp(i * kpoints(j) * (double)cell*a) + 
                    coefs(n)*sWaveFunction(position, atomPositionUp)*
                    exp(i * kpoints(j) * (double)(cell+1)*a) + 
                    coefs(n)*sWaveFunction(position, atomPositionDown)*
                    exp(i * kpoints(j) * (double)(cell-1)*a) + 
                    coefs(n)*sWaveFunction(position, atomPositionUp)*
                    exp(i * kpoints(j) * (double)(cell+2)*a) + 
                    coefs(n)*sWaveFunction(position, atomPositionDown)*
                    exp(i * kpoints(j) * (double)(cell-2)*a);;
        };
    };

    return abs(prob)*abs(prob);
};

/* Instead of using the s hydrogenic wave functions, we are simply
going to calculate the coefficient corresponding to each atom, since
initially we assumed orthonormality (i.e. they are punctual).
It is already modulus squared. 
Input: int n (1D atom position representation). 
int cell (corresponding unit cell) 
cx_mat RScoefs (matrix of real space coefficients as calculated previously)
vec kpoints.
Output: coefficient probability on atom (n, Ncell)*/
double atomCoefficientSquared(int n, int cell, int hCell, const arma::cx_mat& RScoefs,
                              const arma::vec& kpoints){

    std::complex<double> coef = 0;
    double tCoef = 0;
    std::complex<double> j(0,1);
    for(int orb = 0; orb < 64; orb++){
        coef = 0;
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += RScoefs.col(i)(64*n + orb)*exp(j * kpoints(i) * (double)(cell - hCell)* a);
        };
        tCoef += abs(coef)*abs(coef);
    };

    return tCoef;
};