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
        STVH0(a/r0, &SH0);
        return ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(a));
        //return 0.0;
    }
    else{
        STVH0(R, &SH0);
        return ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(R));
        //return 0.0;
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
};

/* Calculate lattice Fourier transform of Keldsyh potential
   Input: double k, int Ncell. Output:  complex double. Vk */
std::complex<double> fourierTrans(double k, int Ncell){
    std::complex<double> i(0,1);
    std::complex<double> Vk = potential(0);
    for(int n = 0; n < Ncell; n++){
        Vk += 2*potential((n+1) * a)*cos(k * (double)(n+1)*a);
    };
    //potential((n+1) * a)*std::exp(i*k*(double)(n+1)*a)
    //2*potential((n+1) * a)*cos(k * (double)(n+1)*a)
    return Vk;
};

std::complex<double> tDirect(std::complex<double> Vk,
                             const arma::cx_vec& coefsK, 
                             const arma::cx_vec& coefsKQ,
                             const arma::cx_vec& coefsK2, 
                             const arma::cx_vec& coefsK2Q, 
                             int Ncell){
    
    std::complex<double> D = 1./(2*Ncell + 1);
    cx_double I_first_pair = arma::cdot(coefsKQ, coefsK2Q);
    cx_double I_second_pair = arma::cdot(coefsK2, coefsK);

    /*if (abs(I_first_pair) < 1E-15){
        I_first_pair = 0.0;
    }
    if (abs(I_second_pair) < 1E-15){
        I_second_pair = 0.0;
    }*/
    D  *= Vk*I_first_pair*I_second_pair;

    return D;
};

std::complex<double> tExchange(std::complex<double> VQ, 
                               const arma::cx_vec& coefsK, 
                               const arma::cx_vec& coefsKQ,
                               const arma::cx_vec& coefsK2, 
                               const arma::cx_vec& coefsK2Q, 
                               int Ncell){
    
    std::complex<double> X = 1./(2*Ncell + 1);
    cx_double I_first_pair = arma::cdot(coefsKQ, coefsK);
    cx_double I_second_pair = arma::cdot(coefsK2, coefsK2Q);

    /*if (abs(I_first_pair) < 1E-15){
        I_first_pair = 0.0;
    }
    if (abs(I_second_pair) < 1E-15){
        I_second_pair = 0.0;
    }*/
    X *= VQ*I_first_pair*I_second_pair;

    return X;
};

mat initializePotentialMatrix(int Ncell, const arma::mat& motif){

    int nAtom = motif.n_rows;
    int dimRows = nAtom*nAtom;
    int dimCols = (2*Ncell+1);

    arma::mat potentialMat = arma::zeros<mat>(dimRows, dimCols);
    arma::vec potentialVector = arma::zeros<vec>(nAtom, 1);

    rowvec cell = arma::zeros(2*Ncell + 1, 3);
    cell.col(1) = arma::regspace(0, 2*Ncell + 1)*a;

    vec ones = arma::zeros(nAtom, 1);
    arma::mat motif_combinations = arma::kron(motif, ones) - arma::kron(ones, motif);

    for(int i = 0; i < 2*Ncell + 1; i++){
        arma::mat position = cell.row(i) - motif_combinations;
        vec pos_module = arma::diagvec(position*position.t());
        for(int j = 0; j < nAtom; j++){
            potentialVector(j) = potential(pos_module(j));
        };
        potentialMat.col(i) = potentialVector;
    };

    cout << "Potential matrix computed" << endl;
    return potentialMat;
};

// Old implementation of potential matrix
/* mat initializePotentialMatrix(int Ncell, const arma::mat& motif){

    int nAtom = motif.n_rows;
    int dimRows = nAtom*nAtom;
    int dimCols = (2*Ncell+1)*(2*Ncell+1);
    mat potentialMat = arma::zeros<mat>(dimRows, dimCols);
    rowvec cellDifference;
    rowvec atomPosDiff;
    double distance;
    for(int i = -Ncell; i < Ncell + 1; i++){
        for(int j = -Ncell; j < Ncell + 1; j++){

            cellDifference = {0.0, (i - j)*a, 0.0};
            for(int n = 0; n < motif.n_rows; n++){
                for(int m = 0; m < motif.n_rows; m++){
                    atomPosDiff = cellDifference - (motif.row(n) - motif.row(m));
                    distance = sqrt(arma::dot(atomPosDiff, atomPosDiff));

                    potentialMat(n*nAtom + m, 
                                (i + Ncell)*(2*Ncell + 1) + (j + Ncell)) = 
                                potential(distance);
                };
            };
        };
    };
    cout << "Potential matrix computed" << endl;
    return arma::kron(potentialMat, arma::ones(64, 1));
};*/

/* Exact implementation of interaction term, valid for both direct and exchange */
std::complex<double> exactInteractionTerm(const arma::cx_vec& coefsK1, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsK3, 
                                     const arma::cx_vec& coefsK4, 
                                     const arma::vec& kArray,
                                     const arma::mat& motif,
                                     const arma::mat& potential,
                                     int Ncell){
                                        
    cx_vec firstCoefArray = conj(coefsK1) % coefsK3;
    cx_vec secondCoefArray = conj(coefsK2) % coefsK4;
    cx_vec coefVector = arma::kron(secondCoefArray, firstCoefArray);
    cout << coefVector.n_rows << "--" << coefVector.n_cols << endl;

    double k_diff = kArray(2)-kArray(0) + kArray(3) - kArray(1);
    std::complex<double> i(0,1);
    cx_vec expArray = arma::regspace<cx_vec>(0, 2*Ncell + 1)*a*k_diff*i;
    expArray = arma::exp(expArray);
    cout << __LINE__ << endl;

    cx_vec result = potential.st()*coefVector;
    result = result % expArray;
    
    cout << __LINE__ << endl;
    std::complex<double> term = arma::sum(result)/(2.*Ncell + 1.);
    return term;
};

/* Routine to calculate the self-energy within the interacting TB theory
to renormalize energies of diagonal elements in the BSE equation */
std::complex<double> selfEnergy(const arma::cx_cube& eigvecStack, 
                                const arma::rowvec& state,
                                const arma::mat& basis,
                                const arma::vec& kpoints,
                                int Ncell, int Qindex){

    // Only for edge states (the self-energy is calculated exclusively
    // between edge states)
    std::complex<double> selfe = 0.0;
    std::complex<double> Xconduction = 0.0;
    std::complex<double> Xvalence = 0.0;
    std::complex<double> Dconduction = 0.0;
    std::complex<double> Dvalence = 0.0;
    std::complex<double> U;
    double k = state(2);
    int nk = kpoints.n_elem;
    int kindex = determineKIndex(k, kpoints);
    double Q = state(3);
    int v = state(0);
    int c = state(1);
    int kQindex;
    if(kindex + Qindex > nk - 1){
        kQindex = kindex + Qindex - nk + 1;
    }
    else{
        kQindex = kindex + Qindex;
    };
    arma::cx_mat slice;
    slice = eigvecStack.slice(kQindex);
    cx_vec cCoefs = slice.col(c);

    slice = eigvecStack.slice(kindex);
    cx_vec vCoefs = slice.col(v);

    int seaState;
    double seaK;
    int basisMultiplicity = basis.n_rows/kpoints.n_elem;

    for(int i = 0; i < (int)basis.n_rows; i++){

        int ik = (int)i/basisMultiplicity;
        slice = eigvecStack.slice(ik);
        seaState = basis.row(i)(0);
        seaK = basis.row(i)(2);
        U = arma::cdot(cCoefs, slice.col(seaState));
        Xconduction = U*conj(U)*fourierTrans(k+Q-seaK, Ncell)/(2.*Ncell + 1);

        U = arma::cdot(vCoefs, slice.col(seaState));
        Xvalence = U*conj(U)*fourierTrans(k-seaK, Ncell)/(2.*Ncell + 1);

        selfe += -Xconduction + Xvalence;
    };

    cout << selfe << endl;
    return selfe;
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

    vec edgeC = arma::regspace(2*(N+1)*5 + 1, 2*(N+1)*5 + 1);
    vec edgeV = arma::regspace(2*(N+1)*5 - 2, 2*(N+1)*5 - 2);

    vec tValence;
    vec tConduction;
    vec auxBands;

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

        /*if(kpoints(i) + Q > PI/a){
            edgeC = arma::regspace(2*(N+1)*5 + 0, 2*(N+1)*5 + 0);

            if(kpoints(i) > PI/a){
                edgeV = arma::regspace(2*(N+1)*5 - 1, 2*(N+1)*5 - 1);
            }
            if(kpoints(i) + Q >= 2*PI/a){
                edgeC = arma::regspace(2*(N+1)*5 + 1, 2*(N+1)*5 + 1);
            };
        };
        tValence = arma::join_cols(valence, edgeV);
        tConduction = arma::join_cols(conduction, edgeC);*/
            

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


/* Compute the basis elements for the spinful exciton problem. Reorders basis
in blocks of defined spin (so that they are diagonal for later calculation of eigenstates of BSE)*/
mat createSOCBasis(int N, double Q, const vec& kpoints, int nbands){

    vec valence = arma::regspace(2*(N+1)*5 - 3 - nbands + 1, 2*(N+1)*5 - 3);
    vec conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2 + nbands - 1);

    int nv = valence.n_elem;
    int nc = conduction.n_elem;
    int nk = kpoints.n_elem;
    int basisDimension = nv*nc*nk;

    arma::mat states = arma::zeros(basisDimension, 4);

    int counter = 0;
    for (int vIndex = 0; vIndex < nv; vIndex++){
        for (int cIndex = 0; cIndex < nc; cIndex++){
            for(int i = 0; i < nk; i++){
                
                mat state = {valence(vIndex), conduction(cIndex), kpoints(i), Q};
                states.row(counter) = state;

                counter++;
            };
        };
    };
    return states;
};

/* Routine to fix the band crossing at the K=PI/a point, intented
to work ONLY for edge states. Does not work for band crossings happening
in other in other k points for bulk bands.
Input: vec eigenval, cx_mat eigenvec
Output: void (updates input references) */
void fixBandCrossing(vec& eigenval, cx_mat& eigenvec, int N){

    double auxEigenval;
    cx_vec auxEigenvec;
    int vband = 2*(N+1)*5 - 2;
    int cband = 2*(N+1)*5 + 1;

    // Valence band
    auxEigenval = eigenval(vband);
	eigenval(vband) = eigenval(vband + 1);
	eigenval(vband + 1) = auxEigenval;

	auxEigenvec = eigenvec.col(vband);
	eigenvec.col(vband) = eigenvec.col(vband + 1);
	eigenvec.col(vband + 1) = auxEigenvec;

	// Conduction band
	auxEigenval = eigenval(cband);
	eigenval(cband) = eigenval(cband - 1);
	eigenval(cband - 1) = auxEigenval;

	auxEigenvec = eigenvec.col(cband);
	eigenvec.col(cband) = eigenvec.col(cband - 1);
	eigenvec.col(cband - 1) = auxEigenvec;

    return;
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

    cx_cube atomicCoefsStack = arma::zeros<cx_cube>(basisDimTB, basisDimTB, nk);
    for (int j = 0; j < nk; j++){
        
        // Initialise vector with form factors
        formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
        for (int k = 0; k < basisDimTB; k++){
            int kAtom = (int)k/8;
            rowvec atomPos = motiv.row(kAtom);
            formFactorArray(k) += std::exp(-i * kpoints(j) * atomPos(1));
        };
        cx_mat formFactorMatrix = arma::kron(arma::ones(1, basisDimTB), formFactorArray);
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
void BShamiltonian(int N, int Ncell, const mat& states, const vec& kpoints, std::string ordering){

    int basisDimBSE = states.n_rows;
    int basisDimTB  = 2*(N+1)*8;
    int nk = kpoints.n_elem;

    int kMultplicity = basisDimBSE/nk;
    int nbands = (int)sqrt(kMultplicity);
    arma::uvec vbands = arma::zeros<arma::uvec>(nbands);
    arma::uvec cbands = arma::zeros<arma::uvec>(nbands);

    // Determine bands that participate
    // Depends on basis ordering
    if(ordering == "kpoints"){
        for(int i = 0; i < nbands; i++){
            vbands(i) = states.row(i)(0);
            cbands(i) = states.row(i*nbands)(1);
        }
    }
    else if(ordering == "spin"){
        for(int i = 0; i < nbands; i++){
            vbands(i) = states.row(i*nk*nbands)(0);
            cbands(i) = states.row(i*nk)(1);
        };
    }
    else{
        throw("Incorrect value of argument ordering, exiting...");
    };

    arma::uvec bands = arma::join_cols(vbands,cbands);
    int nTotalBands = bands.n_elem;

    // Create dictionary that maps bands to indices for storage
    std::map<int, int> bandToIndex;
    for(int i = 0; i < nbands; i++){
        bandToIndex[vbands(i)] = i;
        bandToIndex[cbands(i)] = i + nbands;
    };

    HBS = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK  = arma::zeros(basisDimBSE, basisDimBSE);

    cx_cube eigvecKStack(basisDimTB, nTotalBands, nk);
    cx_cube eigvecKQStack(basisDimTB, nTotalBands, nk);
    mat eigvalKStack(nTotalBands, nk);
    mat eigvalKQStack(nTotalBands, nk);
    cx_vec ftStack = arma::zeros<cx_vec>(nk);

    double k, k2;
    int k_index, kQ_index, k2_index, k2Q_index;
    double Q = states(0, 3);
    std::complex<double> ft, ftD, D, X;
    std::complex<double> ftX = fourierTrans(Q, Ncell);
    cout << "ftX: " << ftX << endl;
    bool doFixCrossing = false;

    // ---------------- Load data in the stack ----------------
    mat motif = createMotiv(N);
    vec auxEigVal(basisDimTB);
    cx_mat auxEigvec(basisDimTB, basisDimTB);
    cx_mat h;

    for (int i = 0; i < nk; i++){
		h = hamiltonian(kpoints(i), H0, Ha, Hsoc);
        arma::eig_sym(auxEigVal, auxEigvec, h);

        if(doFixCrossing && (kpoints(i) > PI/a) && (kpoints(i) < 2*PI/a)){
            fixBandCrossing(auxEigVal, auxEigvec, N);
        };
        
        eigvalKStack.col(i) = auxEigVal(bands);
        eigvecKStack.slice(i) = auxEigvec.cols(bands);

        if(Q != 0){
            h = hamiltonian(kpoints(i) + Q, H0, Ha, Hsoc);
            arma::eig_sym(auxEigVal, auxEigvec, h);

            if(doFixCrossing && (kpoints(i) + Q) > PI/a && (kpoints(i) + Q) < 2*PI/a){
                fixBandCrossing(auxEigVal, auxEigvec, N);
            };

            eigvalKQStack.col(i) = auxEigVal(bands);
            eigvecKQStack.slice(i) = auxEigvec.cols(bands);
        }
        else{
            eigvecKQStack.slice(i) = eigvecKStack.slice(i);
            eigvalKQStack.col(i) = eigvalKStack.col(i);
        };
        // The FT is calculated for vec kpoints starting in zero ALWAYS
        ftStack(i) = fourierTrans(kpoints(i) - kpoints(0), Ncell);
    };
    std::complex<double> ftZero = 0.0;

    // !!!!!!!!!!! Routines have to be fixed
    //cx_cube atomicGCoefsKstack = atomicGCoefs(eigvecKStack, motif, kpoints, N); // Needs to be fixed 
    //cx_cube atomicGCoefsKQstack = atomicGCoefs(eigvecKQStack, motif, kpoints + Q, N);

    cx_cube atomicGCoefsKstack = eigvecKStack;
    cx_cube atomicGCoefsKQstack = eigvecKQStack;

    //int i0 = int(nk/2);k2_index
    for (int i = 0; i < basisDimBSE; i++){

        k = states(i, 2);
        int v = bandToIndex[states(i, 0)];
        int c = bandToIndex[states(i, 1)];

        k_index = determineKIndex(k, kpoints);
        kQ_index = k_index;

        for (int j = i; j < basisDimBSE; j++){

            k2 = states(j, 2);
            int v2 = bandToIndex[states(j, 0)];
            int c2 = bandToIndex[states(j, 1)];

            k2_index = determineKIndex(k2, kpoints);
            k2Q_index = k2_index;

            // Using the atomic gauge
            cx_vec coefsK = atomicGCoefsKstack.slice(k_index).col(v);
            cx_vec coefsKQ = atomicGCoefsKQstack.slice(kQ_index).col(c);
            cx_vec coefsK2 = atomicGCoefsKstack.slice(k2_index).col(v2);
            cx_vec coefsK2Q = atomicGCoefsKQstack.slice(k2Q_index).col(c2);

            int kk2_index = abs(k2_index - k_index); // Always positive
            ftD = ftStack(kk2_index);

            D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);
            X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);
            //X = 0.0;
            //cout << arma::cdot(coefsK, coefsKQ) << endl;
            //cout << arma::cdot(coefsK2, coefsK2Q) << endl;

            //cout << i << j << endl;
            //cout << "Direct: " << D << endl;
            //cout << "Exchange: " << X << endl;

            if (i == j){
                HBS(i, j) = (eigvalKQStack.col(kQ_index)(c) - 
                             eigvalKStack.col(k_index)(v))/2.0 - (D - X)/2.0;
                HK(i, j) = eigvalKQStack(c, kQ_index) - eigvalKStack(v, k_index);
            }
            else{
                HBS(i, j) =  - (D - X);
            };
            if(false){
                cout << "i = " << i << ", j = " << j << "\n" << endl;
                cout << k_index << "--" << kpoints(k_index) << "--" << k << endl;
                cout << k2_index << "--" << kpoints(k2_index) << "--" << k2 << endl;
                cout << kQ_index << "--" << kpoints(kQ_index) << "--" << k+Q << endl;
                cout << k2Q_index << "--" << kpoints(k2Q_index) << "--" << k2+Q << endl;
                cout << kk2_index << "--" << kpoints(kk2_index) << "--" << k-k2 << endl;
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
    //cout << HBS << endl;
    //cout << HBS.submat(0,0 ,nk/2-1, nk/2-1) << endl;;
    //cout << HBS.submat(nk/2, nk/2, nk-1, nk-1) << endl;
    //cout << "Offdiag:" << endl;
    //cout << HBS.submat(0, nk/2, nk/2-1, nk-1) << endl;
    //cout << HBS.submat(nk/2, 0, nk-1, nk/2-1) << endl;
};

/* Compute expected value of tight-binding energy and potential term for 
an exciton eigenstate.
Input: cx_vec eigvec, cx_mat HBS, cx_mat HK. 
Output: Vector 2x1 (<T>, <V>) */
vec computeEnergies(const cx_vec& eigvec, const cx_mat& HBS, const mat& HK){

    std::complex<double> kineticEnergy = arma::cdot(eigvec, HK*eigvec);

    cx_mat HV = HBS - HK;
    std::complex<double> potentialEnergy = arma::cdot(eigvec, HV*eigvec);

    vec energies = {kineticEnergy.real(), potentialEnergy.real()};
    return energies;
};

/* Routine to compute the expected Sz spin value of the electron
and hole that form a given exciton. */
cx_vec spinX(const cx_vec& coefs, 
          const mat& basis, const vec& kpoints, int N){

    double k;

    // Initialize Sz for both electron and hole to zero
    cx_double electronSpin = 0;
    cx_double holeSpin = 0;
    double totalSpin = 0;
    int dimTB = 2*(N+1)*8;
    int dimX = basis.n_rows;
    double Q = basis(0,3);
    cout << Q << endl;

	cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
	cx_vec spinVector = arma::kron(arma::ones(dimTB/8, 1), spinEigvalues);
	cx_vec eigvec, spinEigvec, eigvec2;
    cx_double coefEx;

    // Initialize eigenvectors from H0
    int nk = (int)kpoints.n_elem;
    cx_cube eigvecKStack(dimTB, dimTB, nk);
    cx_cube eigvecKQStack(dimTB, dimTB, nk);

    cx_mat auxEigvec(dimTB, dimTB);
    cx_mat h;
    vec auxEigVal; // We do not want to store the eigenvalues

    for (int i = 0; i < nk; i++){
		h = hamiltonian(kpoints(i), H0, Ha, Hsoc);
        eig_sym(auxEigVal, eigvecKStack.slice(i), h);
        
        if(Q != 0){
            h = hamiltonian(kpoints(i) + Q, H0, Ha, Hsoc);
            eig_sym(auxEigVal, eigvecKQStack.slice(i), h);
        }
        else{
            eigvecKQStack.slice(i) = eigvecKStack.slice(i);
        };
    };

    // Initialize hole spin and electron spin operators
    int nbands = 2;
    int nbandsSq = nbands*nbands;
    vec valence = arma::regspace(2*(N+1)*5 - 3 - nbands + 1, 2*(N+1)*5 - 3);
    vec conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2 + nbands - 1);

<<<<<<< Updated upstream:old/libexcitons.cpp
    arma::cx_mat spinHole = arma::zeros<arma::cx_mat>(dimX, dimX);
    arma::cx_mat spinElectron = arma::zeros<arma::cx_mat>(dimX, dimX);

    arma::cx_mat spinHoleReduced = arma::zeros<arma::cx_mat>(nbands, nbands);
    arma::cx_mat spinElectronReduced = arma::zeros<arma::cx_mat>(nbands, nbands);

    arma::cx_mat vMatrix = arma::eye<cx_mat>(nbands, nbands);
    arma::cx_mat cMatrix = arma::eye<cx_mat>(nbands, nbands);

    for(int k = 0; k < nk; k++){

        for(int i = 0; i < nbands; i++){
            for(int j = 0; j < nbands; j++){
                eigvec = eigvecKStack.slice(k).col(valence(i));
                spinEigvec = eigvec % spinVector;
                eigvec = eigvecKStack.slice(k).col(valence(j));
                spinHoleReduced(i,j) = arma::cdot(eigvec, spinEigvec);

                eigvec = eigvecKQStack.slice(k).col(conduction(i));
                spinEigvec = eigvec % spinVector;
                eigvec = eigvecKQStack.slice(k).col(conduction(j));
                spinElectronReduced(i,j) = arma::cdot(eigvec, spinEigvec);
            }
        }
        spinHole.submat(k*nbandsSq, k*nbandsSq, (k+1)*nbandsSq - 1, (k+1)*nbandsSq - 1) = arma::kron(cMatrix, spinHoleReduced);
        spinElectron.submat(k*nbandsSq, k*nbandsSq, (k+1)*nbandsSq - 1, (k+1)*nbandsSq - 1) = arma::kron(spinElectronReduced, vMatrix);
    }

    // Perform tensor products with the remaining quantum numbers
    // - Valence:

    holeSpin = -arma::cdot(coefs, spinHole*coefs);
    electronSpin = arma::cdot(coefs, spinElectron*coefs);
    totalSpin = real((holeSpin + electronSpin));
    
    cx_vec results = {holeSpin, electronSpin, totalSpin};
=======
        int v = basis(n,0);
        int c = basis(n,1);
        double k = basis(n,2);
        int kIndex = determineKIndex(k, kpoints);

        for(int m = 0; m < (int)basis.n_rows; m++){

            int v2 = basis(m,0);
            int c2 = basis(m,1);
            double k2 = basis(m,2);
            int kIndex2 = determineKIndex(k2, kpoints);

            coefEx = conj(coefs(m))*coefs(n); // BSE coef. squared

            // First valence band (hole) spin
            eigvec = eigvecKStack.slice(kIndex).col(v);
            eigvec2 = eigvecKStack.slice(kIndex2).col(v2);
            spinEigvec = eigvec % spinVector;

            if(k == k2 && c == c2){
                holeSpin -= coefEx*arma::cdot(eigvec2, spinEigvec);
            }
            else{
                holeSpin -= 0;
            };
            //cout << "Valence band: " << v << endl;
            //cout << arma::cdot(eigvec, spinEigvec) << "\n" << endl;

            // Repeat for conduction electron
            eigvec = eigvecKQStack.slice(kIndex).col(c);
            eigvec2 = eigvecKQStack.slice(kIndex2).col(c2);
            spinEigvec = eigvec % spinVector;
            //cout << "Conduction band: " << c << endl;
            //cout << arma::cdot(eigvec, spinEigvec) << "\n" << endl;

            if(k == k2 && v == v2){
                electronSpin += coefEx*arma::cdot(eigvec2, spinEigvec);
            }
            else{
                electronSpin += 0;
            };
        } 
    }
    totalSpin += (abs(holeSpin) + abs(electronSpin));
    
    vec results = {abs(holeSpin), abs(electronSpin), totalSpin};
>>>>>>> Stashed changes:lib/libexcitons.cpp
    return results;
};

// ------------- Routines to compute Fermi Golden Rule -------------

/* Routine to compute the non-interacting electron-hole edge pair 
associated to a given energy. To do so, we run a search algorithm to 
find which k value matches the given energy.
Input: double energy, vec kpoints (provides the search grid)
vec gapEnergy: vector of gap energies associated to each k
Output: vec of coefficients in e-h pair basis associated to desired
pair */
cx_vec ehPairCoefs(double energy, const vec& kpoints, const vec& gapEnergy, bool zone = true){

    cx_vec coefs = arma::zeros<cx_vec>(kpoints.n_elem);
    int closestKindex;
    double eDiff, prevDiff;
    for(int n = 1; n < (int)kpoints.n_elem; n++){

        eDiff = gapEnergy(n) - energy;
        prevDiff = gapEnergy(n-1) - energy;
        if(abs(eDiff) < abs(prevDiff) && eDiff > 0){
            closestKindex = n;
        };
    };
    cout << closestKindex << endl;
    cout << "Selected k: " << kpoints(closestKindex) << "\t" << closestKindex << endl;
    cout << "Closest gap energy: " << gapEnergy(closestKindex) << endl;
    // By virtue of band symmetry, we expect n < nk/2
    int nk = kpoints.n_elem;
    if(zone == true){
        coefs(closestKindex) = 1./sqrt(2);
    }
    else if(zone == false){
        coefs(nk - closestKindex) = 1./sqrt(2);
    }

    cout << "Energy gap (-k): " << gapEnergy(closestKindex) << endl;
    cout << "Energy gap (k): " << gapEnergy(nk - closestKindex) << endl;

    return coefs;
};

/* Routine to compute the Fermi Golden Rule between a bulk exciton 
(the ground state one, but it is prepared for any) and an edge excitation
without interaction matching the initial energy.
Input: cx_vec coefs of bulk exciton, double Ei energy of bulk exc,
int nbands bulk bands that participate (nbands valence and nband conduction)
Output: double transition rate from bulk exciton to edge one */
double fermiGoldenRule(const cx_vec& initialCoefs, double initialE, int nbands, double Q, 
                       int N, int Ncell, const vec& kpoints){

    double transitionRate = 0;
    int nk = kpoints.n_elem;
    int nEdgeStates = nk; 
    mat bulkBasis = createBasis(N, Q, kpoints, nbands, 0);
    mat edgeBasis = createBasis(N, Q, kpoints, 0, nEdgeStates);
    int basisDimTB = 2*(N+1)*8;
    int nedge = edgeBasis.n_rows;
    int nbulk = bulkBasis.n_rows;

    cx_mat W = arma::zeros<cx_mat>(nedge, nbulk);

    cx_cube eigvecKStack(basisDimTB, basisDimTB, nk);
    cx_cube eigvecKQStack(basisDimTB, basisDimTB, nk);
    mat eigvalKStack(basisDimTB, nk);
    mat eigvalKQStack(basisDimTB, nk);
    cx_vec ftStack = arma::zeros<cx_vec>(nk);

    double k, k2;
    int ke_index, keQ_index, kb_index, kbQ_index;
    std::complex<double> ft, ftD, D, X;
    std::complex<double> ftX = fourierTrans(Q, Ncell);
    bool doFixCrossing = false;

    // ---------------- Load data in the stack ----------------
    mat motif = createMotiv(N);
    vec auxEigVal(basisDimTB);
    cx_mat auxEigvec(basisDimTB, basisDimTB);
    cx_mat h;

    int kBulkMultplicity = nbulk/nk;
    int kEdgeMultplicity = nedge/nk;

    for (int i = 0; i < nk; i++){
		h = hamiltonian(kpoints(i), H0, Ha, Hsoc);
        arma::eig_sym(auxEigVal, eigvecKStack.slice(i), h);

        if(doFixCrossing && (kpoints(i) > PI/a) && (kpoints(i) < 2*PI/a)){
            fixBandCrossing(auxEigVal, eigvecKStack.slice(i), N);
        };
        
        eigvalKStack.col(i) = auxEigVal;

        if(Q != 0){
            h = hamiltonian(kpoints(i) + Q, H0, Ha, Hsoc);
            arma::eig_sym(auxEigVal, eigvecKQStack.slice(i), h);

            if(doFixCrossing && (kpoints(i) + Q) > PI/a && (kpoints(i) + Q) < 2*PI/a){
                fixBandCrossing(auxEigVal, auxEigvec, N);
            };

            eigvalKQStack.col(i) = auxEigVal;
        }
        else{
            eigvecKQStack.slice(i) = eigvecKStack.slice(i);
            eigvalKQStack.col(i) = eigvalKStack.col(i);
        };
        // The FT is calculated for vec kpoints starting in zero ALWAYS
        ftStack(i) = fourierTrans(kpoints(i) - kpoints(0), Ncell);
    };
    std::complex<double> ftZero = 0.0;

    // !!!!!!!!!!! Routines have to be fixed
    //cx_cube atomicGCoefsKstack = atomicGCoefs(eigvecKStack, motif, kpoints, N); // Needs to be fixed 
    //cx_cube atomicGCoefsKQstack = atomicGCoefs(eigvecKQStack, motif, kpoints + Q, N);

    cx_cube atomicGCoefsKstack = eigvecKStack;
    cx_cube atomicGCoefsKQstack = eigvecKQStack;

    //int i0 = int(nk/2);k2_index
    // -------- Main loop (W initialization) --------
    for (int i = 0; i < nedge; i++){

        k = edgeBasis(i, 2);
        int v = edgeBasis(i, 0);
        int c = edgeBasis(i, 1);

        ke_index = (int)i/kEdgeMultplicity;
        keQ_index = ke_index;

        for (int j = 0; j < nbulk; j++){

            k2 = bulkBasis(j, 2);
            int v2 = bulkBasis(j, 0);
            int c2 = bulkBasis(j, 1);

            kb_index = (int)j/kBulkMultplicity;
            kbQ_index = kb_index;

            // Using the atomic gauge
            cx_vec coefsK = atomicGCoefsKstack.slice(ke_index).col(v);
            cx_vec coefsKQ = atomicGCoefsKQstack.slice(keQ_index).col(c);
            cx_vec coefsK2 = atomicGCoefsKstack.slice(kb_index).col(v2);
            cx_vec coefsK2Q = atomicGCoefsKQstack.slice(kbQ_index).col(c2);

            int kbke_index; // Always positive
            if(ke_index > kb_index){
                kbke_index = -(kb_index - ke_index);
            }
            else{
                kbke_index = kb_index - ke_index;
            }
            ftD = ftStack(kbke_index);

            D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);
            X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q, Ncell);

            // Compute matrix elements
            W(i, j) = - (D - X);

            if(false){
                cout << "i = " << i << ", j = " << j << "\n" << endl;
                cout << ke_index << "--" << kpoints(ke_index) << "--" << k << endl;
                cout << kb_index << "--" << kpoints(kb_index) << "--" << k2 << endl;
                cout << keQ_index << "--" << kpoints(keQ_index) << "--" << k+Q << endl;
                cout << kbQ_index << "--" << kpoints(kbQ_index) << "--" << k2+Q << endl;
                cout << kbke_index << "--" << kpoints(kbke_index) << "--" << k-k2 << endl;
                //cout << Q_index << "--" << kpoints(Q_index) << "--" << Q << endl;
            };
        };
    };
    
    // Compute gap energy between conduction and valence edge bands
    uword vband = 2*(N+1)*5-2; // Unsigned integer 
    uword cband = 2*(N+1)*5+1;
    uvec bands = {vband, cband};

    vec gapEnergy = (eigvalKStack.row(cband) - eigvalKStack.row(vband)).t();

    cx_vec ehCoefs1 = ehPairCoefs(initialE, kpoints, gapEnergy, true);
    cx_vec ehCoefs2 = ehPairCoefs(initialE, kpoints, gapEnergy, false);

    mat edgeBands = eigvalKStack.rows(bands);

    double delta = 0.01;

    double rho = densityOfStates(initialE, delta, edgeBands);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E10-16; // Units are eV*s
    cout << arma::cdot(ehCoefs1, W*initialCoefs) << endl;
    cout << arma::cdot(ehCoefs2, W*initialCoefs) << endl;
    transitionRate = 2*PI*(pow(abs(arma::cdot(ehCoefs1, W*initialCoefs)),2) + pow(abs(arma::cdot(ehCoefs2, W*initialCoefs)),2))*rho/hbar;

    return transitionRate;
};
