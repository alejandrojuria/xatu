#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "Zigzag.hpp"
#include "Exciton.hpp"

using namespace arma;
using namespace std::chrono;


// Constructor
Exciton::Exciton(int N, int Ncell, double Q, int nBulkBands, int nEdgeBands, vec specifyEdges) 
        : Zigzag(N){
    
    // Initialize basic attributes
    this->Ncell = Ncell;
    this->nk = 2*Ncell - 2;
    this->Q = Q;
    this->nBulkBands = nBulkBands;
    this->nEdgeBands = nEdgeBands;
    this->basisDimTB = 2*(N+1)*8;
    this->specifyEdges = specifyEdges;

    // Initialize derived attributes
    std::cout << "Creating BZ mesh... " << std::flush;
    createMesh();
    std::cout << "Initializing basis for BSE... " << std::flush;
    initializeBasis();
    generateBandDictionary();
    std::cout << "Diagonalizing H0 for all k points... " << std::flush;
    initializeResultsH0();
    std::cout << "Correctly initialized Exciton object" << std::endl;
};

// Destructor
Exciton::~Exciton(){
    std::cout << "Deleting exciton object... " << std::endl;
}

/*      =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x ò 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= 
Source: http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mstvh0_cpp.txt */
void Exciton::STVH0(double X, double *SH0) {
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


/* Calculate value of interaction potential (Keldysh). Units are eV.
   Input: double k. Output: complex double */
double Exciton::potential(double r){
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


/* Calculate lattice Fourier transform of Keldsyh potential
   Input: double k, int Ncell. Output:  complex double. Vk */
std::complex<double> Exciton::fourierTrans(double k){
    std::complex<double> i(0,1);
    std::complex<double> Vk = potential(0);
    for(int n = 0; n < Ncell; n++){
        Vk += 2*potential((n+1) * a)*cos(k * (double)(n+1)*a);
    };
    //potential((n+1) * a)*std::exp(i*k*(double)(n+1)*a)
    //2*potential((n+1) * a)*cos(k * (double)(n+1)*a)
    return Vk;
};

std::complex<double> Exciton::tDirect(std::complex<double> Vk,
                             const arma::cx_vec& coefsK, 
                             const arma::cx_vec& coefsKQ,
                             const arma::cx_vec& coefsK2, 
                             const arma::cx_vec& coefsK2Q)
                             {
    
    std::complex<double> D = 1./(2*Ncell + 1);
    cx_double I_first_pair = arma::cdot(coefsKQ, coefsK2Q);
    cx_double I_second_pair = arma::cdot(coefsK2, coefsK);

    if (abs(I_first_pair) < 1E-15){
        I_first_pair = 0.0;
    }
    if (abs(I_second_pair) < 1E-15){
        I_second_pair = 0.0;
    }
    D  *= Vk*I_first_pair*I_second_pair;

    return D;
};

std::complex<double> Exciton::tExchange(std::complex<double> VQ, 
                               const arma::cx_vec& coefsK, 
                               const arma::cx_vec& coefsKQ,
                               const arma::cx_vec& coefsK2, 
                               const arma::cx_vec& coefsK2Q)
                               {
    
    std::complex<double> X = 1./(2*Ncell + 1);
    cx_double I_first_pair = arma::cdot(coefsKQ, coefsK);
    cx_double I_second_pair = arma::cdot(coefsK2, coefsK2Q);

    if (abs(I_first_pair) < 1E-15){
        I_first_pair = 0.0;
    }
    if (abs(I_second_pair) < 1E-15){
        I_second_pair = 0.0;
    }
    X *= VQ*I_first_pair*I_second_pair;

    return X;
};


/* Method to create mesh of brillouin zone. Note that the last point will be removed, so the vector is 
initiallized with nk + 1 points */
void Exciton::createMesh(){
    vec kpoints = arma::linspace(0.0, 2*PI/a, nk + 2);
    kpoints = kpoints(arma::span(1, nk));

    this->kpoints = kpoints;
};


/* Overloading createBasis method to work with class attributes instead of given ones
*/
void Exciton::initializeBasis(){
    this->basisStates = createBasis(nBulkBands, nEdgeBands);
};

/* Initilise basis to be used in the construction of the BSE matrix.
   We consider only the latest valence band, and the first conduction band.
*/
arma::mat Exciton::createBasis(int nBulkBands, int nEdgeBands){

    int nk = kpoints.n_elem;
    vec valence, conduction;

    // Sanity check
    if (nBulkBands % 2 != 0 ){
        throw std::invalid_argument("Number of bulk bands has to be even");
    };

    if (nBulkBands == 0){
        valence = {};
        conduction = {};
    }
    else{
        valence = arma::regspace(2*(N+1)*5 - 3 - nBulkBands + 1, 2*(N+1)*5 - 3);
        conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2 + nBulkBands - 1);
    };

    if (nEdgeBands == 0){
        edgeC = {};
        edgeV = {};
    }
    else if(!specifyEdges.is_empty() && nEdgeBands == 1){
        edgeC = arma::regspace(2*(N+1)*5 + 0 + specifyEdges(1), 2*(N+1)*5 + 0 + specifyEdges(1));
        edgeV = arma::regspace(2*(N+1)*5 - 2 + specifyEdges(0), 2*(N+1)*5 - 2 + specifyEdges(0));
    }
    else{
        edgeC = arma::regspace(2*(N+1)*5 + 0, 2*(N+1)*5 + nEdgeBands - 1);
        edgeV = arma::regspace(2*(N+1)*5 - 2, 2*(N+1)*5 - 2 + nEdgeBands - 1); 
    };

    vec tValence;
    vec tConduction;
    vec auxBands;

    int nv = valence.n_elem;
    int nc = conduction.n_elem;

    int basisDimension = nv*nc*nk + nk*edgeC.n_elem*edgeV.n_elem + 
                         nk*edgeC.n_elem*nv + 
                         nk*edgeV.n_elem*nc;

    mat states = arma::zeros(basisDimension, 4);

    int it = 0;
    for (int i = 0; i < nk; i++){

        if (nEdgeBands > 0){
            tValence = arma::join_cols(valence, edgeV);
            tConduction = arma::join_cols(edgeC, conduction);
        }
        else{
            tValence = valence;
            tConduction = conduction;
        };
        
        for (int k = 0; k < (int)tConduction.n_elem; k++){
            for (int j = 0; j < (int)tValence.n_elem; j++){
            
                arma::mat state = { tValence(j), tConduction(k), kpoints(i), Q };
                states.row(it) = state;

                it++;
            };
        };
    };

    bandList = arma::sort(arma::conv_to<uvec>::from(arma::join_cols(tValence, tConduction)));
    basisStates = states;
    valenceBands = tValence;
    conductionBands = tConduction;

    return states;
};


/* Compute the basis elements for the spinful exciton problem. Reorders basis
in blocks of defined spin (so that they are diagonal for later calculation of eigenstates of BSE)
Generates only basis for bulk excitons */
void Exciton::createSOCBasis(){

    vec valence = arma::regspace(2*(N+1)*5 - 3 - nBulkBands + 1, 2*(N+1)*5 - 3);
    vec conduction = arma::regspace(2*(N+1)*5 + 2, 2*(N+1)*5 + 2 + nBulkBands - 1);

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

    this->bandList = arma::conv_to<uvec>::from(arma::join_cols(valence, conduction));
    this->basisStates = states;
};

/* Routine to fix the band crossing at the K=PI/a point, intented
to work ONLY for edge states. Does not work for band crossings happening
in other in other k points for bulk bands.
Input: vec eigenval, cx_mat eigenvec
Output: void (updates input references) */
void Exciton::fixBandCrossing(vec& eigenval, cx_mat& eigenvec){

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
int Exciton::determineKIndex(double k){
    int ndiv = kpoints.n_elem - 1;
    if(k > kpoints[ndiv]){
        k -= 2*PI/a;
    };
    return round((k - kpoints(0))*ndiv/(kpoints(ndiv) - kpoints(0)));
};

/* Routine to calculate the coefficients corresponding to wavefunctions
in the atomic gauge.
Input:  */
cx_cube Exciton::atomicGCoefs(const cx_cube& coefs){

    arma::cx_mat formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
    std::complex<double> i(0,1);

    cx_cube atomicCoefsStack = arma::zeros<cx_cube>(basisDimTB, basisDimTB, nk);
    for (int j = 0; j < nk; j++){
        
        // Initialise vector with form factors
        formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
        for (int k = 0; k < basisDimTB; k++){
            int kAtom = (int)k/8;
            rowvec atomPos = motif.row(kAtom);
            formFactorArray(k) += std::exp(-i * kpoints(j) * atomPos(1));
        };
        cx_mat formFactorMatrix = arma::kron(arma::ones(1, basisDimTB), formFactorArray);
        // % stands for element-wise multiplication
        atomicCoefsStack.slice(j) = coefs.slice(j) % formFactorMatrix;
    };

    return atomicCoefsStack;
}

void Exciton::generateBandDictionary(){

    // Determine bands that participate
    // Depends on basis ordering
    int nbands = nBulkBands + nEdgeBands;

    // Create dictionary that maps bands to indices for storage
    std::map<int, int> bandToIndex;
    for(int i = 0; i < nbands; i++){
        bandToIndex[valenceBands(i)] = i;
        bandToIndex[conductionBands(i)] = i + nbands;
    };

    this->bandToIndex = bandToIndex;
};

// Routine to save the relevant data in the stack for later computations
void Exciton::initializeResultsH0(){

    int nTotalBands = 2*(nBulkBands + nEdgeBands);

    cx_cube eigvecKStack(basisDimTB, nTotalBands, nk);
    cx_cube eigvecKQStack(basisDimTB, nTotalBands, nk);
    mat eigvalKStack(nTotalBands, nk);
    mat eigvalKQStack(nTotalBands, nk);
    arma::cx_vec ftStack = arma::zeros<arma::cx_vec>(nk);

    vec auxEigVal(basisDimTB);
    cx_mat auxEigvec(basisDimTB, basisDimTB);
    cx_mat h;

    for (int i = 0; i < nk; i++){
		h = hamiltonian(kpoints(i));
        arma::eig_sym(auxEigVal, auxEigvec, h);
    
        eigvalKStack.col(i) = auxEigVal(bandList);
        eigvecKStack.slice(i) = auxEigvec.cols(bandList);

        if(Q != 0){
            h = hamiltonian(kpoints(i) + Q);
            arma::eig_sym(auxEigVal, auxEigvec, h);

            eigvalKQStack.col(i) = auxEigVal(bandList);
            eigvecKQStack.slice(i) = auxEigvec.cols(bandList);
        }
        else{
            eigvecKQStack.slice(i) = eigvecKStack.slice(i);
            eigvalKQStack.col(i) = eigvalKStack.col(i);
        };
        // The FT is calculated for vec kpoints starting in zero ALWAYS
        ftStack(i) = fourierTrans(kpoints(i) - kpoints(0));
    };

    // !!!!!!!!!!! Routines have to be fixed
    //cx_cube atomicGCoefsKstack = atomicGCoefs(eigvecKStack, motif, kpoints, N); // Needs to be fixed 
    //cx_cube atomicGCoefsKQstack = atomicGCoefs(eigvecKQStack, motif, kpoints + Q, N);

    cx_cube atomicGCoefsKstack = eigvecKStack;
    cx_cube atomicGCoefsKQstack = eigvecKQStack;

    this->eigvalKStack = eigvalKStack;
    this->eigvalKQStack = eigvalKQStack;
    this->eigvecKStack = eigvecKStack;
    this->eigvecKQStack = eigvecKQStack;
    this->ftStack = ftStack;
};


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
void Exciton::BShamiltonian(){

    std::cout << "Initializing Bethe-Salpeter matrix... " << std::flush;

    int basisDimBSE = basisStates.n_rows;

    HBS = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK  = arma::zeros<mat>(basisDimBSE, basisDimBSE);

    std::complex<double> ft;
    std::complex<double> ftX = fourierTrans(Q);
    double threshold = 1E-10;

    #pragma omp parallel for schedule(static, 1) collapse(2)
    for (int i = 0; i < basisDimBSE; i++){
        for (int j = 0; j < basisDimBSE; j++){
            
            double k = basisStates(i, 2);
            int v = bandToIndex[basisStates(i, 0)];
            int c = bandToIndex[basisStates(i, 1)];

            int k_index = determineKIndex(k);
            int kQ_index = k_index;

            double k2 = basisStates(j, 2);
            int v2 = bandToIndex[basisStates(j, 0)];
            int c2 = bandToIndex[basisStates(j, 1)];

            int k2_index = determineKIndex(k2);
            int k2Q_index = k2_index;

            // Using the atomic gauge
            cx_vec coefsK = eigvecKStack.slice(k_index).col(v);
            cx_vec coefsKQ = eigvecKQStack.slice(kQ_index).col(c);
            cx_vec coefsK2 = eigvecKStack.slice(k2_index).col(v2);
            cx_vec coefsK2Q = eigvecKQStack.slice(k2Q_index).col(c2);

            int kk2_index = abs(k2_index - k_index); // Always positive
            std::complex<double> ftD = ftStack(kk2_index);

            std::complex<double> D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q);
            std::complex<double> X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q);

            if(abs(D) < threshold && abs(X) < threshold && i != j){
                continue;
            }

            if (i == j){
                HBS(i, j) = (eigvalKQStack.col(kQ_index)(c) - 
                             eigvalKStack.col(k_index)(v)) - (D - X);
                HK(i, j) = eigvalKQStack(c, kQ_index) - eigvalKStack(v, k_index);
            }
            else{
                HBS(i, j) =  - (D - X);
            };
        };
    };
    //HBS = HBS + HBS.t();
    std::cout << "Done" << std::endl;

    this->HBS = HBS;
    this->HK = HK;
};

/* Compute expected value of tight-binding energy and potential term for 
an exciton eigenstate.
Input: cx_vec eigvec, cx_mat HBS, cx_mat HK. 
Output: Vector 2x1 (<T>, <V>) */
vec Exciton::computeEnergies(const cx_vec& eigvec){

    std::complex<double> kineticEnergy = arma::cdot(eigvec, HK*eigvec);

    cx_mat HV = HBS - HK;
    std::complex<double> potentialEnergy = arma::cdot(eigvec, HV*eigvec);

    vec energies = {kineticEnergy.real(), potentialEnergy.real()};
    return energies;
};

/* Routine to compute the expected Sz spin value of the electron
and hole that form a given exciton. */
cx_vec Exciton::spinX(const cx_vec& coefs){

    // Initialize Sz for both electron and hole to zero
    cx_double electronSpin = 0;
    cx_double holeSpin = 0;
    double totalSpin = 0;
    int dimX = basisStates.n_rows;

	cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
	cx_vec spinVector = arma::kron(arma::ones(basisDimTB/8, 1), spinEigvalues);
	cx_vec eigvec, spinEigvec;

    // Initialize hole spin and electron spin operators
    int nbands = nBulkBands + nEdgeBands;
    int nbandsSq = nbands*nbands;

    arma::cx_mat spinHole = arma::zeros<arma::cx_mat>(dimX, dimX);
    arma::cx_mat spinElectron = arma::zeros<arma::cx_mat>(dimX, dimX);

    arma::cx_mat spinHoleReduced = arma::zeros<arma::cx_mat>(nbands, nbands);
    arma::cx_mat spinElectronReduced = arma::zeros<arma::cx_mat>(nbands, nbands);

    arma::cx_mat vMatrix = arma::eye<cx_mat>(nbands, nbands);
    arma::cx_mat cMatrix = arma::eye<cx_mat>(nbands, nbands);

    for(int k = 0; k < nk; k++){

        for(int i = 0; i < nbands; i++){
            int vIndex = bandToIndex[valenceBands(i)];
            int cIndex = bandToIndex[conductionBands(i)];
            for(int j = 0; j < nbands; j++){
                int vIndex2 = bandToIndex[valenceBands(j)];
                int cIndex2 = bandToIndex[conductionBands(j)];
                eigvec = eigvecKStack.slice(k).col(vIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = eigvecKStack.slice(k).col(vIndex2);
                spinHoleReduced(i,j) = arma::cdot(eigvec, spinEigvec);

                eigvec = eigvecKQStack.slice(k).col(cIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = eigvecKQStack.slice(k).col(cIndex2);
                spinElectronReduced(i,j) = arma::cdot(eigvec, spinEigvec);
            }
        }
        spinHole.submat(k*nbandsSq, k*nbandsSq, (k+1)*nbandsSq - 1, (k+1)*nbandsSq - 1) = arma::kron(cMatrix, spinHoleReduced);
        spinElectron.submat(k*nbandsSq, k*nbandsSq, (k+1)*nbandsSq - 1, (k+1)*nbandsSq - 1) = arma::kron(spinElectronReduced, vMatrix);
    }

    // Perform tensor products with the remaining quantum numbers
    holeSpin = -arma::cdot(coefs, spinHole*coefs);
    electronSpin = arma::cdot(coefs, spinElectron*coefs);
    totalSpin = real((holeSpin + electronSpin));
    
    cx_vec results = {holeSpin, electronSpin, totalSpin};
    return results;
};

// ------------- Routines to compute Fermi Golden Rule -------------

/* Method to compute density of states associated to electron-hole
pairs, particularized to edge bands only */
double Exciton::pairDensityOfStates(double energy, double delta){
    
    double dos = 0;
    for(int n = 0; n < (int)edgeC.n_elem; n++){
        for(int m = 0; m < (int)edgeV.n_elem; m++){
            for(int i = 0; i < (int)kpoints.n_elem; i++){

                uword vband = bandToIndex[edgeV(m)]; // Unsigned integer 
                uword cband = bandToIndex[edgeC(n)];

                double stateEnergy = eigvalKStack.col(i)(cband) - eigvalKStack.col(i)(vband);
                dos += -PI*imag(rGreenF(energy, delta, stateEnergy));
            };
        }
    }
    dos /= kpoints.n_elem*a;

    return dos;
}


/* Method to fix coefficients of degenerate eigenstates according
to expected properties of exciton eigenstates (k, -k symmetric) */
cx_mat Exciton::fixDegeneracyIteration(const cx_vec& eigvec, const cx_vec& eigvec_deg){

        arma::cx_vec rev_eigvec = arma::reverse(eigvec);
        arma::cx_vec rev_eigvec_deg = arma::reverse(eigvec_deg);

        double lowestVal = abs(eigvec(0));
        int index = 0;
        for(int i = 0; i < eigvec.n_elem; i++){
            if (abs(eigvec(i)) < lowestVal){
                lowestVal = abs(eigvec(i));
                index = i;
            };
        };
        cx_double r1 = eigvec(0) - rev_eigvec(nBulkBands*nBulkBands - 1);
        cx_double r2 = eigvec_deg(0) - rev_eigvec_deg(nBulkBands*nBulkBands - 1);
        double r1r2 = abs(r2/r1);
        double alpha = sqrt(r1r2*r1r2/(1 + r1r2*r1r2));
        double beta = sqrt(1 - alpha*alpha);
        cx_vec state = alpha*eigvec + beta*eigvec_deg;
        cx_vec state_deg = beta*eigvec - alpha*eigvec_deg;

        cx_mat states = arma::zeros<cx_mat>(eigvec.n_elem, 2);
        states.col(0) = state;
        states.col(1) = state_deg;

        return states;
}

/* Method to apply several times the degeneracy fixing algorithm
to converge wavefunction */
cx_mat Exciton::fixDegeneracy(const cx_vec& eigvec, const cx_vec& eigvec_deg, int iterations){

    cx_mat states;
    cx_vec state = eigvec;
    cx_vec state_deg = eigvec_deg;
    for(int i = 0; i < iterations; i++){
        states = fixDegeneracyIteration(state, state_deg);
        state = states.col(0);
        state_deg = states.col(1);
        double difference_norm = 0;
        for(int n = 0; n < state.n_elem; n++){
            int bandNumber = n%(nBulkBands*nBulkBands);
            int kindex = n/(nBulkBands*nBulkBands);
            difference_norm += abs(state(n) - state(state.n_elem - 1 - (nBulkBands*nBulkBands - 1 - bandNumber) - kindex*nBulkBands*nBulkBands));
            //cout << state_deg(n) << "--" << state_deg(state.n_elem - 1 - (nBulkBands*nBulkBands - 1 - bandNumber) - kindex*nBulkBands*nBulkBands) << endl;
        };
        cout << "Difference norm: " << sqrt(difference_norm) << endl;
    }
    return states;
};

/* Private method to create an e-h edge state corresponding to a 
wave packets centered in a given kpoint and with a given dispersion */
cx_vec Exciton::wavePacket(double kcenter, double dispersion){

    double coef = 1./(dispersion*sqrt(2*PI));
    cx_vec state = arma::zeros<cx_vec>(nk);
    for(int n = 0; n < nk; n++){
        state(n) = coef*exp(-(kpoints(n) - kcenter)*(kpoints(n) - kcenter)/(dispersion*dispersion));
    }

    return state;
}


/* Routine to compute the non-interacting electron-hole edge pair 
associated to a given energy. To do so, we run a search algorithm to 
find which k value matches the given energy.
Input: double energy, vec kpoints (provides the search grid)
vec gapEnergy: vector of gap energies associated to each k
Output: vec of coefficients in e-h pair basis associated to desired
pair */
cx_vec Exciton::ehPairCoefs(double energy, const vec& gapEnergy, bool zone){

    cx_vec coefs = arma::zeros<cx_vec>(nk);
    int closestKindex = -1;
    double eDiff, prevDiff;
    for(int n = 1; n < (int)kpoints.n_elem/2; n++){

        eDiff = gapEnergy(n) - energy;
        prevDiff = gapEnergy(n-1) - energy;
        if(abs(eDiff) < abs(prevDiff)){
            closestKindex = n;
        };
    };
    cout << closestKindex << endl;
    cout << "Selected k: " << kpoints(closestKindex) << "\t" << closestKindex << endl;
    cout << "Closest gap energy: " << gapEnergy(closestKindex) << endl;
    // By virtue of band symmetry, we expect n < nk/2
    int nk = kpoints.n_elem;
    double dispersion = PI/(16*a);
    if(zone == true){
        //coefs = wavePacket(kpoints(closestKindex), dispersion);
        coefs(closestKindex) = 1.;
    }
    else if(zone == false){
        //coefs = wavePacket(kpoints(nk - closestKindex), dispersion);
        coefs(nk - 1 - closestKindex) = 1.;
    }

    cout << "Energy gap (-k): " << gapEnergy(closestKindex) << endl;
    cout << "Energy gap (k): " << gapEnergy(nk - 1 - closestKindex) << endl;
    this->pairEnergy = gapEnergy(closestKindex);

    return coefs;
};

/* Routine to compute the Fermi Golden Rule between a bulk exciton 
(the ground state one, but it is prepared for any) and an edge excitation
without interaction matching the initial energy.
Input: cx_vec coefs of bulk exciton, double Ei energy of bulk exc,
int nbands bulk bands that participate (nbands valence and nband conduction)
Output: double transition rate from bulk exciton to edge one */
double Exciton::fermiGoldenRule(const cx_vec& initialCoefs, double initialE)
{

    double transitionRate = 0;
    mat bulkBasis = createBasis(nBulkBands, 0);
    mat edgeBasis = createBasis(0, nEdgeBands);
    int nedge = edgeBasis.n_rows;
    int nbulk = bulkBasis.n_rows;

    cx_mat W = arma::zeros<cx_mat>(nedge, nbulk);

    std::complex<double> ft;
    std::complex<double> ftX = fourierTrans(Q);

    // !!!!!!!!!!! Routines have to be fixed
    //cx_cube atomicGCoefsKstack = atomicGCoefs(eigvecKStack, motif, kpoints, N); // Needs to be fixed 
    //cx_cube atomicGCoefsKQstack = atomicGCoefs(eigvecKQStack, motif, kpoints + Q, N);

    cx_cube atomicGCoefsKstack = eigvecKStack;
    cx_cube atomicGCoefsKQstack = eigvecKQStack;

    //int i0 = int(nk/2);k2_index
    // -------- Main loop (W initialization) --------
    //#pragma omp parallel for schedule(static, 1) collapse(2)
    for (int i = 0; i < nedge; i++){
        for (int j = 0; j < nbulk; j++){

            double ke = edgeBasis(i, 2);
            int v = edgeBasis(i, 0);
            int c = edgeBasis(i, 1);

            int ke_index = determineKIndex(ke);
            int keQ_index = ke_index;

            int vIndexEdge = bandToIndex[v];
            int cIndexEdge = bandToIndex[c];

            double kb = bulkBasis(j, 2);
            int v2 = bulkBasis(j, 0);
            int c2 = bulkBasis(j, 1);

            int knumber = j/(nBulkBands*nBulkBands);
            int bandnumber = j%(nBulkBands*nBulkBands);

            int vIndexBulk = bandToIndex[v2];
            int cIndexBulk = bandToIndex[c2];

            int kb_index = determineKIndex(kb);
            int kbQ_index = kb_index;

            /*if(ke > PI/a){
                ke_index = nk - 1 - ke_index;
                keQ_index = ke_index;
            }
            if(kb > PI/a){
                kb_index = nk - 1 - kb_index;
                keQ_index = ke_index;
            }*/

            // Using the atomic gauge
            cx_vec coefsK = atomicGCoefsKstack.slice(ke_index).col(vIndexEdge);
            cx_vec coefsKQ = atomicGCoefsKQstack.slice(keQ_index).col(cIndexEdge);
            cx_vec coefsK2 = atomicGCoefsKstack.slice(kb_index).col(vIndexBulk);
            cx_vec coefsK2Q = atomicGCoefsKQstack.slice(kbQ_index).col(cIndexBulk);

            int kbke_index; // Always positive
            if(ke_index > kb_index){
                kbke_index = -(kb_index - ke_index);
            }
            else{
                kbke_index = kb_index - ke_index;
            }
            std::complex<double> ftD = ftStack(kbke_index);

            std::complex<double> D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q);
            std::complex<double> X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q);

            // Compute matrix elements
            //cout << bandnumber << knumber << endl;
            //cout << nedge - i - 1 << "--" << nbulk - 1 - (nBulkBands*nBulkBands - 1 - bandnumber) - knumber*nBulkBands*nBulkBands << endl;
            W(i, j) = - (D - X);
            //W(nedge - i - 1, nbulk - 1 - (nBulkBands*nBulkBands - 1 - bandnumber) - knumber*nBulkBands*nBulkBands) = conj(W(i,j));
        };
    };

    //cout << W << endl;
    //cout << W.col(3) << endl;
    //cout << W.col(nbulk - 1) << endl;
    //cout << W.row(1) << endl;
    //cout << W.row(nedge - 2) << endl;
    //cout << atomicGCoefsKstack.slice(3).col(bandToIndex[edgeV(0)]) << endl;
    //cout << atomicGCoefsKstack.slice(4).col(bandToIndex[edgeV(0)]) << endl;

    //cx_vec ones = arma::ones<cx_vec>(nbulk, 1);
    //cout << W*ones << endl;

    //cx_rowvec coefss = {0,0,0,1,1,0,0,0};
    //cout << coefss*W << endl;

    //coefss = {0,0,0,0,1,0,0,0};
    //cout << coefss*W << endl;

    //cout << initialCoefs << endl;
    //vec out_vec = abs(W*initialCoefs);
    //cout << out_vec << endl;
    //cout << out_vec(32) << "--" << out_vec(nk - 1 - 32);
    
    // Compute gap energy between conduction and valence edge bands
    uword vband = bandToIndex[edgeV(0)]; // Unsigned integer 
    uword cband = bandToIndex[edgeC(0)];
    uvec bands = {vband, cband};

    vec gapEnergy = (eigvalKStack.row(cband) - eigvalKStack.row(vband)).t();

    cx_vec ehCoefs1 = ehPairCoefs(initialE, gapEnergy, true);
    cx_vec ehCoefs2 = ehPairCoefs(initialE, gapEnergy, false);

    mat edgeBands = eigvalKStack.rows(bands);

    double delta = 0.01;
    double rho = pairDensityOfStates(pairEnergy, delta);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E-16; // Units are eV*s
    //cout << initialCoefs << endl;
    cout << "First overlap: " << pow(abs(arma::cdot(ehCoefs1, W*initialCoefs)),2) << endl;
    cout << "Second overlap: " << pow(abs(arma::cdot(ehCoefs2, W*initialCoefs)),2) << endl;
    transitionRate = 2*PI*(pow(abs(arma::cdot(ehCoefs1, W*initialCoefs)),2) + pow(abs(arma::cdot(ehCoefs2, W*initialCoefs)),2))*rho/hbar;

    return transitionRate;
};
