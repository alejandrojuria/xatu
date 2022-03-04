#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "System.hpp"
#include "GExciton.hpp"

using namespace arma;
using namespace std::chrono;

GExciton::GExciton(){
    throw std::invalid_argument("Error: GExciton must be invoked with at least one parameter (systemfile)");
};

GExciton::GExciton(std::string filename, int ncell, const arma::ivec& bands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          System(filename){

    // Constructor for GExciton class. If bands vector is given, it is used in all
    // calculations instead of parameters nbands, nrmbands.
    // int ncell: Number of unit cells along ONE direction.
    // Initialize basic attributes
    this->ncell_      = ncell;
    this->totalCells_ = pow(ncell, ndim);
    this->Q_          = Q;
    this->bands_      = bands;
    this->eps_m_      = parameters(0);
    this->eps_s_      = parameters(1);
    this->r0_         = parameters(2);
    this->cutoff_     = ncell/2.5; // Default value, seems to preserve crystal point group

    if(r0 == 0){
        throw std::invalid_argument("Error: r0 must be non-zero");
    }

    if (bands.n_elem > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }

    // arma::ivec is implemented with typedef s64
    std::vector<arma::s64> valence, conduction;
    for(int i = 0; i < bands.n_elem; i++){
        if (bands(i) <= 0){
            valence.push_back(bands(i) + fermiLevel);
        }
        else{
            conduction.push_back(bands(i) + fermiLevel);
        }
    }
    this->valenceBands_ = arma::ivec(valence);
    this->conductionBands_ = arma::ivec(conduction);

    if(!bands.empty()){
        std::cout << "Correctly initialized Exciton object" << std::endl;
    }
};

GExciton::GExciton(std::string filename, int ncell, int nbands, int nrmbands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          GExciton(filename, ncell, {}, parameters, Q) {
    
    if (2*nbands > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }
    this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands - nrmbands + 1, 
                                                     fermiLevel - nrmbands);
    this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, 
                                                        fermiLevel + nbands + nrmbands);
    this->bands_ = arma::join_cols(valenceBands, conductionBands) - fermiLevel;
    this->excitonbasisdim_ = nk*valenceBands.n_elem*conductionBands.n_elem;

    std::cout << "Correctly initialized Exciton object" << std::endl;
};


// This constructor is intented to run directly 
GExciton::GExciton(std::string modelfile, std::string excitonfile) : System(modelfile){
    ExcitonConfiguration excitonconfig = ExcitonConfiguration(excitonfile);
    initializeExcitonAttributes(excitonconfig);
    bool useApproximation = excitonconfig.excitonInfo.useApproximation;
    initializeHamiltonian(useApproximation);
}

// Destructor
GExciton::~GExciton(){
    std::cout << "Deleting exciton object... " << std::endl;
}


/* ------------------------------ Setters ------------------------------ */
void GExciton::setUnitCells(int ncell){
    if(ncell > 0){
        ncell_ = ncell;
    }
    else{
        std::cout << "ncell must be a positive number" << std::endl;
    }
}

void GExciton::setBands(const arma::ivec& bands){
    bands_ = bands;
    std::vector<arma::s64> valence, conduction;
    for(int i = 0; i < bands.n_elem; i++){
        if (bands(i) <= 0){
            valence.push_back(bands(i) + fermiLevel);
        }
        else{
            conduction.push_back(bands(i) + fermiLevel);
        }
    }
    this->valenceBands_ = arma::ivec(valence);
    this->conductionBands_ = arma::ivec(conduction);
}

void GExciton::setBands(int nbands, int nrmbands){
    if(nbands > 0 && nrmbands > 0){
        this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands + 1, fermiLevel - nrmbands);
        this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, fermiLevel + nbands);
        this->bands_ = arma::join_rows(valenceBands, conductionBands);
    }
    else{
        std::cout << "Included bands and removed bands must be positive numbers" << std::endl;
    }
}

void GExciton::setQ(const arma::rowvec& Q){
    if(Q.n_elem == 3){
        Q_ = Q;
    }
    else{
        std::cout << "Q vector must be 3d" << std::endl;
    }
    
}

void GExciton::setParameters(const arma::rowvec& parameters){
    if(parameters.n_elem == 3){
        eps_m_ = parameters(0);
        eps_s_ = parameters(1);
        r0_    = parameters(2);
    }
    else{
        std::cout << "parameters array must be 3d (eps_m, eps_s, r0)" << std::endl;
    }
}

void GExciton::setCutoff(double cutoff){
    if(cutoff > 0){
        cutoff_ = cutoff;
        if(cutoff > ncell){
            std::cout << "Warning: cutoff is higher than number of unit cells" << std::endl;
        }
    }
    else{
        std::cout << "cutoff must be a positive number" << std::endl;
    }
}

void GExciton::setParameters(double eps_m, double eps_s, double r0){
    // TODO: Introduce additional comprobations regarding value of parameters (positive)
    eps_m_ = eps_m;
    eps_s_ = eps_s;
    r0_    = r0;
}

/*      =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x ò 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= 
Source: http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mstvh0_cpp.txt */
void GExciton::STVH0(double X, double *SH0) {
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
double GExciton::potential(double r){
    double eps_bar = (eps_m + eps_s)/2;
    double SH0;
    double cutoff = arma::norm(bravaisLattice.row(0)) * cutoff_ + 1E-5;
    double R = abs(r)/r0;
    double potential_value;
    if(r == 0){
        STVH0(a/r0, &SH0);
        potential_value = ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(a/r0));
    }
    else if (r > cutoff){
        potential_value = 0.0;
    }
    else{
        STVH0(R, &SH0);
        potential_value = ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(R));
    };

    return potential_value;
    
};


/* Calculate lattice Fourier transform of Keldsyh potential
   Input: double k, int ncell. Output:  complex double. Vk */
std::complex<double> GExciton::fourierTransform(arma::rowvec k, const arma::mat& cells, bool useApproximation){
    std::complex<double> imag(0,1);
    //std::complex<double> Vk = potential(0);
    std::complex<double> Vk = 0.0;

    if (useApproximation){
        for(int n = 0; n < cells.n_rows; n++){
            arma::rowvec cell = cells.row(n);
            double module = arma::norm(cell);
            Vk += potential(module)*std::exp(imag*arma::dot(k, cell));
	    }
        Vk *= totalCells;
    }

    else{ // This might be wrong if already using truncated cells
        for (int n = 0; n < cells.n_rows; n++){
            arma::rowvec cell = cells.row(n);
            for (int m = 0; m < cells.n_rows; m++){
                arma::rowvec cell2 = cells.row(m);
                double module = arma::norm(cell - cell2);
                Vk += potential(module)*std::exp(imag*arma::dot(k, cell - cell2));
            }
        }
    };
    return Vk;
};

std::complex<double> GExciton::tDirect(std::complex<double> Vk,
                             const arma::cx_vec& coefsK, 
                             const arma::cx_vec& coefsKQ,
                             const arma::cx_vec& coefsK2, 
                             const arma::cx_vec& coefsK2Q)
                             {
    
    std::complex<double> D = 1./pow(totalCells, 2);
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

std::complex<double> GExciton::tExchange(std::complex<double> VQ, 
                               const arma::cx_vec& coefsK, 
                               const arma::cx_vec& coefsKQ,
                               const arma::cx_vec& coefsK2, 
                               const arma::cx_vec& coefsK2Q)
                               {
    
    std::complex<double> X = 1./pow(totalCells, 2);
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

void GExciton::initializeExcitonAttributes(const ExcitonConfiguration& config){
    ncell_ = config.excitonInfo.ncell;
    bands_ = config.excitonInfo.bands;
    Q_ = config.excitonInfo.Q;

};

void GExciton::initializePotentialMatrix(){

    arma::mat cells_coefs = generateCombinations(ncell, ndim);
    int dimRows = natoms*natoms*norbitals*norbitals;
    int dimCols = cells_coefs.n_rows*cells_coefs.n_rows;

    arma::mat potentialMat = arma::zeros<arma::mat>(dimRows, dimCols);
    //arma::vec potentialVector = arma::zeros<vec>(dimRows, 1);
    
    int ncells = cells_coefs.n_rows;
    arma::mat cells = arma::zeros(cells_coefs.n_rows, 3);
    for (int n = 0; n < cells.n_rows; n++){
        arma::rowvec cell_vector = arma::zeros(1, 3);
        for (int i = 0; i < ndim; i++){
            cell_vector += cells_coefs.row(n)(i) * bravaisLattice.row(i);
        }
        cells.row(n) = cell_vector;
    };

    vec ones = arma::ones(natoms*norbitals, 1);

    // With spin
    arma::mat motif = arma::kron(arma::ones(2, 1), arma::kron(this->motif, arma::ones(norbitals/2, 1)));
    
    // Without spin
    //arma::mat motif = arma::kron(this->motif, arma::ones(norbitals, 1));

    arma::mat motif_combinations = arma::kron(motif, ones) - arma::kron(ones, motif);

    #pragma omp parallel for schedule(static, 1) collapse(2)
    for(int i = 0; i < ncells; i++){
        for(int j = 0; j < ncells; j++){
            arma::mat position = arma::kron(cells.row(i) - cells.row(j), arma::ones(dimRows, 1)) - motif_combinations;
            for(int n = 0; n < dimRows; n++){ // Aqui habia puesto solo un natoms -> MAL (?)
                potentialMat.col(i*ncells + j)(n) = potential(arma::norm(position.row(n)));
            };
        };
    };

    std::cout << "Potential matrix computed" << std::endl;
    this->potentialMat = potentialMat;
};

/* Exact implementation of interaction term, valid for both direct and exchange */
std::complex<double> GExciton::exactInteractionTerm(const arma::cx_vec& coefsK1, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsK3, 
                                     const arma::cx_vec& coefsK4, 
                                     const arma::rowvec& k){
    
    cx_vec firstCoefArray = arma::conj(coefsK1) % coefsK3;
    cx_vec secondCoefArray = arma::conj(coefsK2) % coefsK4;
    cx_vec coefVector = arma::kron(secondCoefArray, firstCoefArray);

    std::complex<double> i(0,1);

    arma::mat cells_coefs = generateCombinations(ncell, ndim);
    arma::mat cells = arma::zeros(cells_coefs.n_rows, 3);
    for (int n = 0; n < cells.n_rows; n++){
        arma::rowvec cell_vector = arma::zeros(1, 3);
        for (int j = 0; j < ndim; j++){
            cell_vector += cells_coefs.row(n)(j) * bravaisLattice.row(j);
        }
        cells.row(n) = cell_vector;
    };
    cx_vec expArray = arma::zeros<cx_vec>(cells.n_rows*cells.n_rows);
    for (int n = 0; n < cells.n_rows; n++){
        for(int m = 0; m < cells.n_rows; m++){
            expArray(n*cells.n_rows + m) = std::exp(i*arma::dot(k, cells.row(n) - cells.row(m)));
        }
    }

    cx_vec result = potentialMat.st()*coefVector;
    result = result % expArray;

    double ncells = pow(ncell, ndim);
    
    std::complex<double> term = arma::sum(result)/(ncells*ncells);
    return term;
};

/* Overloading createBasis method to work with class attributes instead of given ones
*/
void GExciton::initializeBasis(){
    this->basisStates_ = createBasis(conductionBands, valenceBands);
    //createSOCBasis();
};

/* Initilise basis to be used in the construction of the BSE matrix.
   We consider only the latest valence band, and the first conduction band.
*/
arma::imat GExciton::createBasis(const arma::ivec& conductionBands, 
                                const arma::ivec& valenceBands){

    arma::imat states = arma::zeros<arma::imat>(excitonbasisdim, 3);
    int it = 0;
    for (int i = 0; i < nk; i++){
        for (int k = 0; k < (int)conductionBands.n_elem; k++){
            for (int j = 0; j < (int)valenceBands.n_elem; j++){

                arma::irowvec state = { valenceBands(j), conductionBands(k), i };
                states.row(it) = state;
                it++;
            };
        };
    };

    basisStates_ = states;
    bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));

    return states;
};

/* Method to generate a basis which is a subset of the basis considered for the
exciton. Its main purpose is to allow computation of Fermi golden rule between
two specified subbasis. */
arma::imat GExciton::specifyBasisSubset(const arma::ivec& bands){

    // Check if given bands vector corresponds to subset of bands
    try{
        for (const auto& band : bands){
            for (const auto& reference_band : bandList){
                if ((band + fermiLevel - reference_band) == 0) {
                    continue;
                }
            }
            throw "Error: Given band list must be a subset of the exciton one"; 
        };
    }
    catch (std::string e){
        std::cerr << e;
    };

    int reducedBasisDim = nk*bands.n_elem;
    std::vector<arma::s64> valence, conduction;
    for(int i = 0; i < bands.n_elem; i++){
        if (bands(i) <= 0){
            valence.push_back(bands(i) + fermiLevel);
        }
        else{
            conduction.push_back(bands(i) + fermiLevel);
        }
    }
    arma::ivec valenceBands = arma::ivec(valence);
    arma::ivec conductionBands = arma::ivec(conduction);

    arma::imat states = createBasis(conductionBands, valenceBands);

    return states;
}


/* Compute the basis elements for the spinful exciton problem. Reorders basis
in blocks of defined spin (so that they are diagonal for later calculation of eigenstates of BSE)
Generates only basis for bulk excitons.
This routine does not work with vec bands (for now) -> Now it works? (overloading constructors) */
void GExciton::useSpinfulBasis(){

    arma::imat states = arma::zeros<arma::imat>(excitonbasisdim, 3);

    int counter = 0;
    for (int vIndex = 0; vIndex < valenceBands.n_elem; vIndex++){
        for (int cIndex = 0; cIndex < conductionBands.n_elem; cIndex++){
            for(int i = 0; i < nk; i++){
                
                arma::irowvec state = {valenceBands(vIndex), conductionBands(cIndex), i};
                states.row(counter) = state;

                counter++;
            };
        };
    };

    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));
    this->basisStates_ = states;
};


// /* Routine to fix the band crossing at the K=PI/a point, intented
// to work ONLY for edge states. Does not work for band crossings happening
// in other in other k points for bulk bands.
// Input: vec eigenval, cx_mat eigenvec
// Output: void (updates input references) */
// void GExciton::fixBandCrossing(vec& eigenval, cx_mat& eigenvec){

//     double auxEigenval;
//     cx_vec auxEigenvec;
//     int vband = 2*(N+1)*5 - 2;
//     int cband = 2*(N+1)*5 + 1;

//     // Valence band
//     auxEigenval = eigenval(vband);
// 	eigenval(vband) = eigenval(vband + 1);
// 	eigenval(vband + 1) = auxEigenval;

// 	auxEigenvec = eigenvec.col(vband);
// 	eigenvec.col(vband) = eigenvec.col(vband + 1);
// 	eigenvec.col(vband + 1) = auxEigenvec;

// 	// Conduction band
// 	auxEigenval = eigenval(cband);
// 	eigenval(cband) = eigenval(cband - 1);
// 	eigenval(cband - 1) = auxEigenval;

// 	auxEigenvec = eigenvec.col(cband);
// 	eigenvec.col(cband) = eigenvec.col(cband - 1);
// 	eigenvec.col(cband - 1) = auxEigenvec;

//     return;
// };


// /* Routine to calculate the index i associated to k within the 
// kpoints vector.
// Input: double k, vec kpoints. Output: index i */
// int GExciton::determineKIndex(double k){
//     int ndiv = kpoints.n_elem - 1;
//     if(k > kpoints[ndiv]){
//         k -= 2*PI/a;
//     };
//     return round((k - kpoints(0))*ndiv/(kpoints(ndiv) - kpoints(0)));
// };


// NOTE: This one may be important (or not, apparently not so much)
// /* Routine to calculate the coefficients corresponding to wavefunctions
// in the atomic gauge.
// Input:  */
// cx_cube GExciton::atomicGCoefs(const cx_cube& coefs){

//     arma::cx_mat formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
//     std::complex<double> i(0,1);

//     cx_cube atomicCoefsStack = arma::zeros<cx_cube>(basisDimTB, basisDimTB, nk);
//     for (int j = 0; j < nk; j++){
        
//         // Initialise vector with form factors
//         formFactorArray = arma::zeros<arma::cx_vec>(basisDimTB, 1);
//         for (int k = 0; k < basisDimTB; k++){
//             int kAtom = (int)k/8;
//             rowvec atomPos = motif.row(kAtom);
//             formFactorArray(k) += std::exp(-i * kpoints(j) * atomPos(1));
//         };
//         cx_mat formFactorMatrix = arma::kron(arma::ones(1, basisDimTB), formFactorArray);
//         // % stands for element-wise multiplication
//         atomicCoefsStack.slice(j) = coefs.slice(j) % formFactorMatrix;
//     };

//     return atomicCoefsStack;
// }

void GExciton::generateBandDictionary(){


    // Create dictionary that maps bands to indices for storage
    std::map<int, int> bandToIndex;
    for(int i = 0; i < bandList.n_elem; i++){
        bandToIndex[bandList(i)] = i;
    };

    this->bandToIndex = bandToIndex;
};

// Routine to save the relevant data in the stack for later computations
void GExciton::initializeResultsH0(){

    int nTotalBands = bandList.n_elem;
    double radius = arma::norm(bravaisLattice.row(0)) * cutoff_;
    // By default
    arma::mat cells = truncateSupercell(ncell, radius);

    cx_cube eigvecKStack(basisdim, nTotalBands, nk);
    cx_cube eigvecKQStack(basisdim, nTotalBands, nk);
    arma::mat eigvalKStack(nTotalBands, nk);
    arma::mat eigvalKQStack(nTotalBands, nk);
    arma::cx_mat ftStack(nk, nk);

    vec auxEigVal(basisdim);
    cx_mat auxEigvec(basisdim, basisdim);
    cx_mat h;

    // Progress bar variables
    int step = 1;
	int displayNext = step;
	int percent = 0;

    for (int i = 0; i < nk; i++){
		h = hamiltonian(kpoints.row(i));
        arma::eig_sym(auxEigVal, auxEigvec, h);
    
        eigvalKStack.col(i) = auxEigVal(bandList);
        eigvecKStack.slice(i) = auxEigvec.cols(bandList);

        if(arma::norm(Q) != 0){
            h = hamiltonian(kpoints.row(i) + Q);
            arma::eig_sym(auxEigVal, auxEigvec, h);

            eigvalKQStack.col(i) = auxEigVal(bandList);
            eigvecKQStack.slice(i) = auxEigvec.cols(bandList);
        }
        else{
            eigvecKQStack.slice(i) = eigvecKStack.slice(i);
            eigvalKQStack.col(i) = eigvalKStack.col(i);
        };
        // The FT is calculated for vec kpoints starting in zero ALWAYS
        #pragma omp parallel for
        for (int j = 0; j < nk; j++){
            ftStack(i, j) = fourierTransform(kpoints.row(i) - kpoints.row(j), cells, true);
        }

        // Formatted progress indicator
		percent = (100 * (i + 1)) / nk ;
		if (percent >= displayNext){
            cout << "\r" << "[" << std::string(percent / 5, '|') << std::string(100 / 5 - percent / 5, ' ') << "]";
            cout << percent << "%";
            std::cout.flush();
            displayNext += step;
        }
    };
    std::cout << std::endl;
    // ftStack(0) = 0.0;

    // !!!!!!!!!!! Routines have to be fixed
    //cx_cube atomicGCoefsKstack = atomicGCoefs(eigvecKStack, motif, kpoints, N); // Needs to be fixed 
    //cx_cube atomicGCoefsKQstack = atomicGCoefs(eigvecKQStack, motif, kpoints + Q, N);

    cx_cube atomicGCoefsKstack = eigvecKStack;
    cx_cube atomicGCoefsKQstack = eigvecKQStack;

    this->eigvalKStack_ = eigvalKStack;
    this->eigvalKQStack_ = eigvalKQStack;
    this->eigvecKStack_ = eigvecKStack;
    this->eigvecKQStack_ = eigvecKQStack;
    this->ftStack = ftStack;
    this->ftX = fourierTransform(Q, cells, false);
};


/* Routine to initialize the required variables to construct the Bethe-Salpeter Hamiltonian */
void GExciton:: initializeHamiltonian(bool useApproximation){

    if(bands.empty()){
        throw std::invalid_argument("Error: Exciton object must have some bands");
    }
    if(nk == 0){
        throw std::invalid_argument("Error: BZ mesh must be initialized first");
    }

    this->excitonbasisdim_ = nk*valenceBands.n_elem*conductionBands.n_elem;

    std::cout << "Initializing basis for BSE... " << std::flush;
    initializeBasis();
    generateBandDictionary();

    std::cout << "Diagonalizing H0 for all k points... " << std::endl;
    initializeResultsH0();

    if (!useApproximation){
        std::cout << "Calculating potential for all lattice positions... " << std::flush;
        initializePotentialMatrix();
    };
}


/* Initialize BSE hamiltonian matrix and kinetic matrix. Recursive approach:
Instead of calculating the energies and coeficients dinamically, which
is too expensive, instead we first calculate those for each k, save them
in the stack, and then call them consecutively as we build the matrix.
Analogously, we calculate the Fourier transform of the potential beforehand,
saving it in the stack so that it can be later called in the matrix element
calculation.
Input: int N (cells finite direction), vec states, int ncells (periodic 
direction), int nEdgeStates. Output: None (updates previously declared matrices) 
BEWARE: Does not work for Q < 0 (Expected to use reflection symmetry)*/
void GExciton::BShamiltonian(const arma::imat& basis, bool useApproximation){

    std::cout << "Initializing Bethe-Salpeter matrix... " << std::flush;

    arma::imat basisStates = this->basisStates;
    if (!basis.is_empty()){
        basisStates = basis;
    };

    int basisDimBSE = basisStates.n_rows;
    std::cout << "BSE dimension : " << basisDimBSE << std::endl;

    HBS_ = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK_   = arma::zeros<arma::mat>(basisDimBSE, basisDimBSE);

    std::complex<double> ft;
    // std::complex<double> ftX = fourierTransform(Q);
    double threshold = 1E-10;

    // To be able to parallelize over the triangular matrix, we build
    int loopLength = basisDimBSE*(basisDimBSE + 1)/2.;

    // https://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients
    #pragma omp parallel for
    for (int n = 0; n < loopLength; n++){
        int ii = loopLength - 1 - n;
        int k  = floor((sqrt(8*ii + 1) - 1)/2);
        int i = basisDimBSE - 1 - k;
        int j = basisDimBSE - 1 - ii + k*(k+1)/2;
    
        double k_index = basisStates(i, 2);
        int v = bandToIndex[basisStates(i, 0)];
        int c = bandToIndex[basisStates(i, 1)];

        int kQ_index = k_index;

        double k2_index = basisStates(j, 2);
        int v2 = bandToIndex[basisStates(j, 0)];
        int c2 = bandToIndex[basisStates(j, 1)];

        int k2Q_index = k2_index;

        // Using the atomic gauge
        cx_vec coefsK = eigvecKStack.slice(k_index).col(v);
        cx_vec coefsKQ = eigvecKQStack.slice(kQ_index).col(c);
        cx_vec coefsK2 = eigvecKStack.slice(k2_index).col(v2);
        cx_vec coefsK2Q = eigvecKQStack.slice(k2Q_index).col(c2);

        std::complex<double> D, X;
        if (useApproximation){
            // int kk2_index = k2_index - k_index; // Always positive
            std::complex<double> ftD;
            ftD = ftStack(k2_index, k_index);
            D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q);
            X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q);
        }
        else{
            arma::rowvec kD = kpoints.row(k2_index) - kpoints.row(k_index);
            D = exactInteractionTerm(coefsKQ, coefsK2, coefsK2Q, coefsK, kD);
            X = exactInteractionTerm(coefsKQ, coefsK2, coefsK, coefsK2Q, Q);
        };

        if (i == j){
            HBS_(i, j) = (eigvalKQStack.col(kQ_index)(c) - 
                            eigvalKStack.col(k_index)(v))/2. - (D - X)/2.;
            HK_(i, j) = eigvalKQStack(c, kQ_index) - eigvalKStack(v, k_index);
            
        }
        else{
            HBS_(i, j) =  - (D - X);
        };
        
    }
        
    HBS_ = HBS + HBS.t();
    std::cout << "Done" << std::endl;
};

// Routine to diagonalize the BSE and return a Result object
Result GExciton::diagonalize(){
    std::cout << "Diagonalizing BSE... " << std::flush;
    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::eig_sym(eigval, eigvec, HBS);
    std::cout << "Done" << std::endl;

    Result results = Result(*this, eigval, eigvec);

    return results;
}

// ------------- Routines to compute Fermi Golden Rule -------------

/* Method to compute density of states associated to electron-hole
pairs, particularized to edge bands only */
double GExciton::pairDensityOfStates(const arma::ivec& valence, const arma::ivec& conduction, double energy, double delta){
    
    double dos = 0;
    for(int n = 0; n < (int)valence.n_elem; n++){
        for(int m = 0; m < (int)conduction.n_elem; m++){
            for(int i = 0; i < (int)kpoints.n_elem; i++){

                uword vband = bandToIndex[valence(n)]; // Unsigned integer 
                uword cband = bandToIndex[conduction(m)];
                //cx_vec coefK = eigvecKStack.slice(i).col(vband);
                //cx_vec coefKQ = eigvecKQStack.slice(i).col(cband);
                //double ftD = abs(ftStack(0));
                //double D = abs(tDirect(ftD, coefK, coefKQ, coefK, coefKQ));
                //cout << D << endl;

                double stateEnergy = eigvalKStack.col(i)(cband) - eigvalKStack.col(i)(vband);
                dos += -PI*imag(rGreenF(energy, delta, stateEnergy));
            };
        }
    }
    //if(!specifyEdges.is_empty() && specifyEdges(0) != specifyEdges(1)){ // Add degenerate bands
    //    dos *= 4;
    //}
    dos /= (kpoints.n_elem*a);

    return dos;
}


// /* Method to fix coefficients of degenerate eigenstates according
// to expected properties of exciton eigenstates (k, -k symmetric) */
// cx_mat GExciton::fixDegeneracyIteration(const cx_vec& eigvec, const cx_vec& eigvec_deg){

//         arma::cx_vec rev_eigvec = arma::reverse(eigvec);
//         arma::cx_vec rev_eigvec_deg = arma::reverse(eigvec_deg);

//         double lowestVal = abs(eigvec(0));
//         int index = 0;
//         for(int i = 0; i < eigvec.n_elem; i++){
//             if (abs(eigvec(i)) < lowestVal){
//                 lowestVal = abs(eigvec(i));
//                 index = i;
//             };
//         };
//         cx_double r1 = eigvec(0) - rev_eigvec(eigvec.n_elem- 1);
//         cx_double r2 = eigvec_deg(0) - rev_eigvec_deg(eigvec.n_elem - 1);
//         double r1r2 = abs(r2/r1);
//         double alpha = sqrt(r1r2*r1r2/(1 + r1r2*r1r2));
//         double beta = sqrt(1 - alpha*alpha);
//         cx_vec state = alpha*eigvec + beta*eigvec_deg;
//         cx_vec state_deg = beta*eigvec - alpha*eigvec_deg;

//         cx_mat states = arma::zeros<cx_mat>(eigvec.n_elem, 2);
//         states.col(0) = state;
//         states.col(1) = state_deg;

//         return states;
// }

// /* Method to apply several times the degeneracy fixing algorithm
// to converge wavefunction */
// cx_mat GExciton::fixDegeneracy(const cx_vec& eigvec, const cx_vec& eigvec_deg, int iterations){

//     cx_mat states;
//     cx_vec state = eigvec;
//     cx_vec state_deg = eigvec_deg;
//     for(int i = 0; i < iterations; i++){
//         states = fixDegeneracyIteration(state, state_deg);
//         state = states.col(0);
//         state_deg = states.col(1);
//         double difference_norm = 0;
//         for(int n = 0; n < state.n_elem; n++){
//             int bandNumber = n%(nBulkBands*nBulkBands);
//             int kindex = n/(nBulkBands*nBulkBands);
//             difference_norm += abs(state(n) - state(state.n_elem - 1 - (nBulkBands*nBulkBands - 1 - bandNumber) - kindex*nBulkBands*nBulkBands));
//             //cout << state_deg(n) << "--" << state_deg(state.n_elem - 1 - (nBulkBands*nBulkBands - 1 - bandNumber) - kindex*nBulkBands*nBulkBands) << endl;
//         };
//         cout << "Difference norm: " << sqrt(difference_norm) << endl;
//     }
//     return states;
// };

/* Private method to create an e-h edge state corresponding to a 
wave packets centered in a given kpoint and with a given dispersion */
cx_vec GExciton::wavePacket(double kcenter, double dispersion){

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
cx_vec GExciton::ehPairCoefs(double energy, const vec& gapEnergy, bool zone){

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
double GExciton::fermiGoldenRule(const cx_vec& initialCoefs, double initialE)
{

    double transitionRate = 0;
    // Previously here I had bands={} in the routines below
    arma::imat bulkBasis = specifyBasisSubset({});
    arma::imat edgeBasis = specifyBasisSubset({});
    int nedge = edgeBasis.n_rows;
    int nbulk = bulkBasis.n_rows;
    int nBulkBands = nbands - nrmbands;
    int nEdgeBands = nrmbands;

    cx_mat W = arma::zeros<cx_mat>(nedge, nbulk);

    std::complex<double> ft;
    //std::complex<double> ftX = fourierTransform(Q);

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

            double ke_index = edgeBasis(i, 2);
            int v = edgeBasis(i, 0);
            int c = edgeBasis(i, 1);

            int keQ_index = ke_index;

            int vIndexEdge = bandToIndex[v];
            int cIndexEdge = bandToIndex[c];

            double kb_index = bulkBasis(j, 2);
            int v2 = bulkBasis(j, 0);
            int c2 = bulkBasis(j, 1);

            int knumber = j/(nBulkBands*nBulkBands);
            int bandnumber = j%(nBulkBands*nBulkBands);

            int vIndexBulk = bandToIndex[v2];
            int cIndexBulk = bandToIndex[c2];

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

            // Aproximated interaction
            std::complex<double> D = tDirect(ftD, coefsK, coefsKQ, coefsK2, coefsK2Q);
            std::complex<double> X = tExchange(ftX, coefsK, coefsKQ, coefsK2, coefsK2Q);

            // Exact interaction (TO BE IMPLEMENTED YET)
            //std::complex<double> D = exactInteractionTerm(coefsK, coefsK2, coefsKQ, coefsK2Q, kArray, motif, potentialMatrix, ncell)
            //std::complex<double> X = exactInteractionTerm(coefsK, coefsK2, coefsKQ, coefsK2Q, )

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
    uword vband = bandToIndex[fermiLevel]; // Unsigned integer 
    uword cband = bandToIndex[fermiLevel + nrmbands];
    uvec bands = {vband, cband};

    vec gapEnergy = (eigvalKStack.row(cband) - eigvalKStack.row(vband)).t();

    cx_vec ehCoefs1 = ehPairCoefs(initialE, gapEnergy, true);
    cx_vec ehCoefs2 = ehPairCoefs(initialE, gapEnergy, false);

    arma::mat edgeBands = eigvalKStack.rows(bands);

    double delta = 2.4/(2*ncell); // Adjust delta depending on number of k points
    double rho = pairDensityOfStates(valenceBands, conductionBands, pairEnergy, delta);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E-16; // Units are eV*s
    //cout << initialCoefs << endl;
    cout << "First t.r. (-k): " << 2*PI*pow(abs(arma::cdot(ehCoefs1, W*initialCoefs)),2)*rho/hbar << endl;
    cout << "Second t.r. (k): " << 2*PI*pow(abs(arma::cdot(ehCoefs2, W*initialCoefs)),2)*rho/hbar << endl;
    transitionRate = 2*PI*(pow(abs(arma::cdot(ehCoefs1, W*initialCoefs)),2) + pow(abs(arma::cdot(ehCoefs2, W*initialCoefs)),2))*rho/hbar;

    return transitionRate;
};
