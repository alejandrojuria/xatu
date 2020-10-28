#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "zigzag.hpp"
#include "excitons.hpp"

using namespace arma;
using namespace std::chrono;

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
    double eps = 43.;
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


std::complex<double> realfourierTrans(double k, int Ncell){
    double eps = 1.;
    double eps1 = 1.;
    double eps2 = 5.06;
    double eps_bar = (eps1 + eps2)/2;
    double r0 = c*eps/(eps1 + eps2);
    double cte = ec/(2*eps0*eps_bar);
    double kappak = 1 + r0*k;

    if(k == 0){
        kappak = 1 + r0*0.001;
        return (cte*1E-10/(0.001*kappak));
    }
    else{
        return cte*1E-10/(k*kappak);
    };

}

/* Calculate lattice Fourier transform of Keldsyh potential
   Input: double k, int Ncell. Output:  complex double. Vk */
std::complex<double> fourierTrans(double k, int Ncell){
    std::complex<double> i(0,1);
    std::complex<double> Vk = potential(0);
    for(int n = 0; n < Ncell; n++){
        Vk += potential((n+1) * a)*std::exp(i * k * (double)(n+1)*a)
            + potential(- (n+1) * a)*std::exp(- i * k * (double)(n+1)*a);
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
arma::mat createBasis(int N, double Q, const vec& kpoints, int nbands, int nEdgeStates){

    int nk = kpoints.n_elem;
    vec valence;
    vec conduction;

    // Sanity check
    if (nEdgeStates > nk){
        throw std::invalid_argument("More edge states than k points");
    };
    if (nbands % 2 != 0 ){
        throw std::invalid_argument("Number of bulk bands has to be even");
    };

    if (nbands == 0){
        valence = {};
        conduction = {};
    }
    else{
        valence = arma::regspace(2*(N+1)*5 - 3 - nbands + 1, 2*(N+1)*5 - 3);
        conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2 + nbands - 1);
    };

    vec edgeC = arma::regspace(2*(N+1)*5 + 0, 2*(N+1)*5 + 0);
    vec edgeV = arma::regspace(2*(N+1)*5 - 2, 2*(N+1)*5 - 2);

    vec tValence;
    vec tConduction;

    int nv = valence.n_elem;
    int nc = conduction.n_elem;

    int basisDimension = nv*nc*nk + nEdgeStates*edgeC.n_elem*edgeV.n_elem + 
                         nEdgeStates*edgeC.n_elem*nv + 
                         nEdgeStates*edgeV.n_elem*nc;

    mat states = arma::zeros(basisDimension, 4);

    int it = 0;
    for (int i = 0; i < nk; i++){

        if ((nEdgeStates > 0) && 
            (abs(i - (int)(nk/2)) <= (int)(nEdgeStates/2))){
            
            tValence = arma::join_cols(valence, edgeV);
            tConduction = arma::join_cols(conduction, edgeC);

        }
        else{
            tValence = valence;
            tConduction = conduction;
        };

        for (int k = 0; k < (int)tConduction.n_elem; k++){
            for (int j = 0; j < (int)tValence.n_elem; j++){
            
                mat state = { tValence(j), tConduction(k), kpoints(i), Q };
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
    if(k > kpoints[ndiv]){
        k -= 2*PI/a;
    };
    return round((k - kpoints(0))*ndiv/(kpoints(ndiv) - kpoints(0)));
};

/* Routine to calculate the coefficients corresponding to wavefunctions
in the atomic gauge.
Input:  */
cx_cube atomicGCoefs(const cx_cube& coefs, const mat& motiv, const vec& kpoints, int N){
    int basisDimTB = motiv.n_rows*8;
    int nk = coefs.n_slices;
    arma::cx_mat formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
    std::complex<double> i(0,1);

    cx_cube atomicCoefsStack = zeros<cx_cube>(basisDimTB, basisDimTB, nk);
    for (int j = 0; j < nk; j++){
        
        // Initialise vector with form factors
        formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
        for (int k = 0; k < basisDimTB; k++){
            int kAtom = (int)k/8;
            rowvec atomPos = motiv.row(kAtom);
            formFactorArray(k) += std::exp(-i * kpoints(j) * atomPos(0));
        };
        cx_mat formFactorMatrix = arma::kron(arma::ones(1, basisDimTB), formFactorArray).st();
        // % stands for element-wise multiplication
        atomicCoefsStack.slice(j) = coefs.slice(j) % formFactorMatrix;
    };

    return atomicCoefsStack;
}

/* Initialize BSE hamiltonian matrix and kinetic matrix. Recursive approach:
Instead of calculating the energies and coeficients dinamically, which
is too expensive, instead we first calculate those for each k, save them
in the stack, and then call them consecutively as we build the matrix.
Analogously, we calculate the Fourier transform of the potential beforehand,
saving it in the stack so that it can be later called in the matrix element
calculation.
Input: int N (cells finite direction), vec states, int Ncells (periodic 
direction), int nEdgeStates. Output: None (updates previously declared matrices) 
BEWARE: Does not work for Q < 0 (Expected to use reflection symmetry)*/
void BShamiltonian(int N, int Ncell, const mat& states, const vec& kpoints, 
                   bool doApproximation){

    int basisDimBSE = states.n_rows;
    int basisDimTB  = 2*(N+1)*8;
    int nk = kpoints.n_elem;

    HBS = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK  = arma::zeros(basisDimBSE, basisDimBSE);

    cx_cube eigvecKStack(basisDimTB, basisDimTB, nk);
    mat eigvalKStack(basisDimTB, nk);
    cx_vec ftStack = arma::zeros<cx_vec>(nk);

    double k, k2;
    int k_index, kQ_index, k2_index, k2Q_index;
    double Q = states(0, 3);

    int Q_index = determineKIndex(Q, kpoints);
    std::complex<double> ft, D, X;
    std::complex<double> ftX = fourierTrans(Q, Ncell);

    // ---------------- Load data in the stack ----------------
    vec auxEigVal(basisDimTB);
    for (int i = 0; i < nk; i++){
        arma::eig_sym(auxEigVal, eigvecKStack.slice(i), hamiltonian(kpoints(i), H0, Ha, Hsoc));
 
        eigvalKStack.col(i) = auxEigVal; 
        //ftStack(i) = fourierTrans(kpoints(i), Ncell);
    };
    std::complex<double> ftZero = 0.0;

    mat motiv = createMotiv(N);
    cx_cube atomicGCoefsStack = atomicGCoefs(eigvecKStack, motiv, kpoints, N);
    // cx_cube atomicGCoefsStack = eigvecKStack;
    //int i0 = int(nk/2);
    //ftStack(i0) = 0.0;

    // ----------- Calculate matrix elements of HBS, HK -----------
    for (int i = 0; i < basisDimBSE; i++){

        k = states(i, 2);
        int v = states(i, 0);
        int c = states(i, 1);

        k_index = determineKIndex(k, kpoints);
        if (Q != 0){
            kQ_index = determineKIndex(k + Q, kpoints);
        }
        else{   
            kQ_index = k_index;
        }

        for (int j = i; j < basisDimBSE; j++){

            k2 = states(j, 2);
            int v2 = states(j, 0);
            int c2 = states(j, 1);

            k2_index = determineKIndex(k2, kpoints);
            if (Q != 0){
                k2Q_index = determineKIndex(k2 + Q, kpoints);
            }
            else{
                k2Q_index = k2_index;
            };

            // Using the atomic gauge
            cx_vec coefsK = atomicGCoefsStack.slice(k_index).col(v);
            cx_vec coefsKQ = atomicGCoefsStack.slice(kQ_index).col(c);
            cx_vec coefsK2 = atomicGCoefsStack.slice(k2_index).col(v2);
            cx_vec coefsK2Q = atomicGCoefsStack.slice(k2Q_index).col(c2);

            int kk2_index = determineKIndex(abs(k2 - k), kpoints);
            std::complex<double> ftD = fourierTrans(k2 - k, Ncell);
            
            if(k == k2){
                ftD = ftZero;
            }
            else{
                //ft = ftStack(kk2_index);
            };

            D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);

            if(Q == 0){
                X = 0;
            }
            else{
                X = tExchange(ftX,
                             coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);
            }
            if (i == j){
                HBS(i, j) = (eigvalKStack(states(i, 1), kQ_index) - eigvalKStack(states(i, 0), k_index))/2
                            - (D - X);
                HK(i, j) = eigvalKStack(states(i, 1), kQ_index) - eigvalKStack(states(i, 0), k_index);

            }
            else{
                HBS(i, j) =  - (D - X);
            };
            if(true){
                cout << "i = " << i << ", j = " << j << "\n" << endl;
                cout << k_index << "--" << kpoints(k_index) << "--" << k << endl;
                cout << k2_index << "--" << kpoints(k2_index) << "--" << k2 << endl;
                cout << kQ_index << "--" << kpoints(kQ_index) << "--" << k+Q << endl;
                cout << k2Q_index << "--" << kpoints(k2Q_index) << "--" << k2+Q << endl;
                //cout << kk2_index << "--" << kpoints(kk2_index) << "--" << abs(k-k2) << endl;
                //cout << Q_index << "--" << kpoints(Q_index) << "--" << Q << endl;
                if(i == 97 && j == 99){
                    cout << "Test FT" << endl;
                    cout << ftStack(k2Q_index) << endl;
                    cout << fourierTrans(k2 + Q, Ncell) << endl;
                };
            };
        };
    };
    HBS = HBS + HBS.t();
    cout << HBS.submat(0,0 ,nk/2-1, nk/2-1) << endl;;
    cout << HBS.submat(nk/2, nk/2, nk-1, nk-1) << endl;
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
