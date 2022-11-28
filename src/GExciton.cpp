#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iomanip>

#include "System.hpp"
#include "GExciton.hpp"
#include "utils.hpp"
#include "davidson.hpp"

using namespace arma;
using namespace std::chrono;

void GExciton::initializeExcitonAttributes(int ncell, const arma::ivec& bands, 
                                      const arma::rowvec& parameters, const arma::rowvec& Q){
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
}

void GExciton::initializeExcitonAttributes(const ExcitonConfiguration& cfg){

    int ncell        = cfg.excitonInfo.ncell;
    int nbands       = cfg.excitonInfo.nbands;
    arma::ivec bands = cfg.excitonInfo.bands;
    arma::rowvec parameters = {cfg.excitonInfo.eps(0), cfg.excitonInfo.eps(1), cfg.excitonInfo.eps(2)};
    arma::rowvec Q   = cfg.excitonInfo.Q;

    if (bands.empty()){
        bands = arma::regspace<arma::ivec>(- nbands + 1, nbands);
    }

    initializeExcitonAttributes(ncell, bands, parameters, Q);

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

GExciton::GExciton(const SystemConfiguration& config, int ncell, const arma::ivec& bands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
                  System(config) {

    initializeExcitonAttributes(ncell, bands, parameters, Q);

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

GExciton::GExciton(const SystemConfiguration& config, int ncell, int nbands, int nrmbands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          GExciton(config, ncell, {}, parameters, Q) {
    
    if (2*nbands > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }
    this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands - nrmbands + 1, 
                                                     fermiLevel - nrmbands);
    this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, 
                                                        fermiLevel + nbands + nrmbands);
    this->bands_ = arma::join_cols(valenceBands, conductionBands) - fermiLevel;

    std::cout << "Correctly initialized Exciton object" << std::endl;
};


// This constructor is intented to run directly 
GExciton::GExciton(const SystemConfiguration& config, const ExcitonConfiguration& excitonConfig) : System(config){
    initializeExcitonAttributes(excitonConfig);
    std::cout << "Correctly initialized Exciton object" << std::endl;
}

GExciton::GExciton(const System& system, int ncell, const arma::ivec& bands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
                  System(system) {

    initializeExcitonAttributes(ncell, bands, parameters, Q);

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

GExciton::GExciton(const System& system, int ncell, int nbands, int nrmbands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          GExciton(system, ncell, {}, parameters, Q) {
    
    if (2*nbands > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }

    this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands - nrmbands + 1, 
                                                     fermiLevel - nrmbands);
    this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, 
                                                        fermiLevel + nbands + nrmbands);
    this->bands_ = arma::join_cols(valenceBands, conductionBands) - fermiLevel;

    std::cout << "Correctly initialized Exciton object" << std::endl;
};

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

void GExciton::setParameters(double eps_m, double eps_s, double r0){
    // TODO: Introduce additional comprobations regarding value of parameters (positive)
    eps_m_ = eps_m;
    eps_s_ = eps_s;
    r0_    = r0;
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

void GExciton::setGauge(std::string gauge){
    if(gauge != "lattice" && gauge != "atomic"){
        throw std::invalid_argument("setGauge(): gauge must be either lattice or atomic");
    }
    this->gauge_ = gauge;
}

void GExciton::setMode(std::string mode){
    if(mode != "realspace" && mode != "reciprocalspace"){
        throw std::invalid_argument("setMode(): mode must be either realspace or reciprocalspace");
    }
    this->mode_ = mode;
}

void GExciton::setReciprocalVectors(int nReciprocalVectors){
    if(nReciprocalVectors < 0){
        throw std::invalid_argument("setReciprocalVectors(): given number must be positive");
    }
    this->nReciprocalVectors_ = nReciprocalVectors;
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

/*---------------------------------------- Fourier transforms ----------------------------------------*/

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

    Vk /= pow(totalCells, 2);
    return Vk;
};

double GExciton::analyticFourierTransform(arma::rowvec q){
    double radius = cutoff*arma::norm(reciprocalLattice.row(0));
    double potential = 0;
    double eps_bar = (eps_m + eps_s)/2;
    double eps = 1E-5;
    double unitCellArea = 5.41;

    double qnorm = arma::norm(q);
    if (qnorm < eps){
        potential = 0;
    }
    else{
        potential = 1/(qnorm*(1 + r0*qnorm));
    }
    
    potential = potential*ec*1E10/(2*eps0*eps_bar*unitCellArea*totalCells);
    return potential;
}

double GExciton::fourierTransformFromCoefs(const arma::vec& Acoefs, const arma::vec& Dcoefs, 
                                           const arma::rowvec& kpoint_crystal, int Ncut){
    double VF = 0.0;
    double trigm_cos,trigm_sin, trign_cos,trign_sin, dA, dD;
    double kx_temp = kpoint_crystal(0);
    double ky_temp = kpoint_crystal(1);
    for (int m = 0; m < Ncut; m++ ){ //column
        for (int n = 0; n < Ncut; n++ ){ // row
            trigm_cos = cos(2* M_PI * m * kx_temp); // cos(2pi m qx)
            trigm_sin = sin(2* M_PI * m * kx_temp); // sin(2pi m qx)
            trign_cos = cos(2* M_PI * n * ky_temp); // cos(2pi m qy)
            trign_sin = sin(2* M_PI * n * ky_temp); // sin(2pi m qy)
            dA = Acoefs[n + Ncut*m]*trigm_cos*trign_cos;
            dD = Dcoefs[n + Ncut*m]*trigm_sin*trign_sin;
            VF= VF + dA + dD;
        }
    }

    VF *= (27.2); // Conversion to eV
    return VF;
}

/* Routine to compute the lattice Fourier transform with the potential displaced by some
vectors of the motif */
std::complex<double> GExciton::motifFourierTransform(int fAtomIndex, int sAtomIndex, const arma::rowvec& k, const arma::mat& cells){

    std::complex<double> imag(0,1);
    std::complex<double> Vk = 0.0;
    arma::rowvec firstAtom = motif.row(fAtomIndex).subvec(0, 2);
    arma::rowvec secondAtom = motif.row(sAtomIndex).subvec(0, 2);

    for(int n = 0; n < cells.n_rows; n++){
        arma::rowvec cell = cells.row(n);
        double module = arma::norm(cell + firstAtom - secondAtom);
        Vk += potential(module)*std::exp(imag*arma::dot(k, cell));
    }
    Vk /= pow(totalCells, 1);

    return Vk;
}

arma::cx_mat GExciton::extendMotifFT(const arma::cx_mat& motifFT){
    arma::cx_mat extendedMFT = arma::zeros<arma::cx_mat>(this->basisdim, this->basisdim);
    int rowIterator = 0;
    int colIterator = 0;
    for(unsigned int atom_index_r = 0; atom_index_r < this->motif.n_rows; atom_index_r++){
        int species_r = motif.row(atom_index_r)(3);
        int norbitals_r = orbitals(species_r);
        colIterator = 0;
        for(unsigned int atom_index_c = 0; atom_index_c < this->motif.n_rows; atom_index_c++){
            int species_c = motif.row(atom_index_c)(3);
            int norbitals_c = orbitals(species_c);
            extendedMFT.submat(rowIterator, colIterator, 
                               rowIterator + norbitals_r - 1, colIterator + norbitals_c - 1) = 
                          motifFT(atom_index_r, atom_index_c) * arma::ones(norbitals_r, norbitals_c);
            colIterator += norbitals_c;
        }
        rowIterator += norbitals_r;
    }

    return extendedMFT;
}


/*------------------------------------ Interaction matrix elements ------------------------------------*/

/* Exact implementation of interaction term, valid for both direct and exchange */
std::complex<double> GExciton::exactInteractionTermMFT(const arma::cx_vec& coefsK1, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsK3, 
                                     const arma::cx_vec& coefsK4,
                                     const arma::cx_mat& motifFT){
    
    arma::cx_vec firstCoefArray = arma::conj(coefsK1) % coefsK3;
    arma::cx_vec secondCoefArray = arma::conj(coefsK2) % coefsK4;
    std::complex<double> term = arma::dot(firstCoefArray, extendMotifFT(motifFT) * secondCoefArray);
    
    //arma::cx_vec coefVector = arma::kron(firstCoefArray, secondCoefArray);
    //std::complex<double> term = arma::dot(coefVector, extendMotifFT(motifFT));
    
    return term;
};

/* Exact implementation of interaction term, valid for both direct and exchange */
std::complex<double> GExciton::interactionTermFT(const arma::cx_vec& coefsK, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsKQ, 
                                     const arma::cx_vec& coefsK2Q,
                                     const arma::rowvec& k, 
                                     const arma::rowvec& k2,
                                     const arma::rowvec& kQ, 
                                     const arma::rowvec& k2Q,
                                     int nrcells){
    
    std::complex<double> Ic, Iv;
    std::complex<double> term = 0;
    double radius = cutoff * arma::norm(reciprocalLattice.row(0));
    arma::mat reciprocalVectors = truncateReciprocalSupercell(nrcells, radius);

    for(int i = 0; i < reciprocalVectors.n_rows; i++){
        auto G = reciprocalVectors.row(i);

        Ic = blochCoherenceFactor(coefsK2Q, coefsKQ, k2Q, kQ, G);
        Iv = blochCoherenceFactor(coefsK2, coefsK, k2, k, G);

        term += conj(Ic)*Iv*analyticFourierTransform(k-k2+G);
    }

    return term;
};

std::complex<double> GExciton::blochCoherenceFactor(const arma::cx_vec& coefs1, const arma::cx_vec& coefs2,
                                                    const arma::rowvec& k1, const arma::rowvec& k2,
                                                    const arma::rowvec& G){

    std::complex<double> imag(0, 1);
    arma::cx_vec coefs = arma::conj(coefs1) % coefs2;
    arma::cx_vec phases(basisdim);
    for(int i = 0; i < natoms; i++){
        int species = motif.row(i)(3);
        arma::rowvec atomPosition = motif.row(i).subvec(0, 2);
        phases.subvec(i*orbitals(species), (i+1)*orbitals(species) - 1) = 
        exp(imag*arma::dot(k1 - k2 - G, atomPosition));
    }

    std::complex<double> factor = arma::dot(coefs, phases);

    return factor;
}

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

// NOTE: This one may be important (or not, apparently not so much)
// /* Routine to calculate the coefficients corresponding to wavefunctions
// in the atomic gauge.
// Input:  */
arma::cx_vec GExciton::latticeToAtomicGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

    arma::cx_vec phases(basisdim);
    std::complex<double> imag(0, 1);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.row(atomIndex)(3);
        for(int orbIndex = 0; orbIndex < orbitals(species); orbIndex++){
            arma::rowvec atomPosition = motif.row(atomIndex).subvec(0, 2);
            phases(it) = exp(-imag*arma::dot(k, atomPosition));
            it++;
        }
    }

    arma::cx_vec atomicCoefs = coefs % phases;
    return atomicCoefs;
}

arma::cx_vec GExciton::atomicToLatticeGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

    arma::cx_vec phases(basisdim);
    std::complex<double> imag(0, 1);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.row(atomIndex)(3);
        for(int orbIndex = 0; orbIndex < orbitals(species); orbIndex++){
            arma::rowvec atomPosition = motif.row(atomIndex).subvec(0, 2);
            phases(it) = exp(-imag*arma::dot(k, atomPosition));
            it++;
        }
    }

    arma::cx_vec atomicCoefs = coefs % phases;
    return atomicCoefs;
}

arma::cx_mat GExciton::fixGlobalPhase(arma::cx_mat& coefs){
    arma::cx_rowvec sums = arma::sum(coefs);
    //arma::cout << "Before fix: " << sums << arma::endl;
    std::complex<double> imag(0, 1);
    for(int j = 0; j < sums.n_elem; j++){
        double phase = arg(sums(j));
        coefs.col(j) *= exp(-imag*phase);
    }

    return coefs;
}

void GExciton::generateBandDictionary(){


    // Create dictionary that maps bands to indices for storage
    std::map<int, int> bandToIndex;
    for(int i = 0; i < bandList.n_elem; i++){
        bandToIndex[bandList(i)] = i;
    };

    this->bandToIndex = bandToIndex;
};

void GExciton::initializeMotifFT(int i, const arma::mat& cells){
    // Uses hermiticity of V
    for(int fAtomIndex = 0; fAtomIndex < natoms; fAtomIndex++){
        for(int sAtomIndex = fAtomIndex; sAtomIndex < natoms; sAtomIndex++){
            ftMotifStack(fAtomIndex, sAtomIndex, i) = 
            motifFourierTransform(fAtomIndex, sAtomIndex, kpoints.row(i), cells);
            ftMotifStack(sAtomIndex, fAtomIndex, i) = conj(ftMotifStack(fAtomIndex, sAtomIndex, i));
        }   
    }
}

// Routine to save the relevant data in the stack for later computations
void GExciton::initializeResultsH0(bool triangular){

    int nTotalBands = bandList.n_elem;
    double radius = arma::norm(bravaisLattice.row(0)) * cutoff_;
    arma::mat cells = truncateSupercell(ncell, radius);

    this->eigvecKStack_  = arma::cx_cube(basisdim, nTotalBands, nk);
    this->eigvecKQStack_ = arma::cx_cube(basisdim, nTotalBands, nk);
    this->eigvalKStack_  = arma::mat(nTotalBands, nk);
    this->eigvalKQStack_ = arma::mat(nTotalBands, nk);
    this->ftMotifStack   = arma::cx_cube(natoms, natoms, nk);

    vec auxEigVal(basisdim);
    cx_mat auxEigvec(basisdim, basisdim);
    cx_mat h;

    // Progress bar variables
    int step = 1;
	int displayNext = step;
	int percent = 0;

    calculateInverseReciprocalMatrix();
    std::complex<double> imag(0, 1);

    for (int i = 0; i < nk; i++){
        arma::rowvec k = kpoints.row(i);
        solveBands(k, auxEigVal, auxEigvec, triangular);

        auxEigvec = fixGlobalPhase(auxEigvec);
        eigvalKStack_.col(i) = auxEigVal(bandList);
        eigvecKStack_.slice(i) = auxEigvec.cols(bandList);

        if(arma::norm(Q) != 0){
            arma::rowvec kQ = kpoints.row(i) + Q;
            solveBands(kQ, auxEigVal, auxEigvec, triangular);

            auxEigvec = fixGlobalPhase(auxEigvec);
            eigvalKQStack_.col(i) = auxEigVal(bandList);
            eigvecKQStack_.slice(i) = auxEigvec.cols(bandList);
        }
        else{
            eigvecKQStack_.slice(i) = eigvecKStack.slice(i);
            eigvalKQStack_.col(i) = eigvalKStack.col(i);
        };

        // BIGGEST BOTTLENECK OF THE CODE
        if(this->mode == "realspace"){
            initializeMotifFT(i, cells);     
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

    std::cout << "\nDone" << std::endl;
};

/* Routine to initialize the required variables to construct the Bethe-Salpeter Hamiltonian */
void GExciton::initializeHamiltonian(bool triangular){

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
    initializeResultsH0(triangular);
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
void GExciton::BShamiltonian(const arma::imat& basis){

    arma::imat basisStates = this->basisStates;
    if (!basis.is_empty()){
        basisStates = basis;
    };

    int basisDimBSE = basisStates.n_rows;
    std::cout << "BSE dimension: " << basisDimBSE << std::endl;
    std::cout << "Initializing Bethe-Salpeter matrix... " << std::flush;

    HBS_ = arma::zeros<cx_mat>(basisDimBSE, basisDimBSE);
    HK_   = arma::zeros<arma::mat>(basisDimBSE, basisDimBSE);
    
    // To be able to parallelize over the triangular matrix, we build
    int loopLength = basisDimBSE*(basisDimBSE + 1)/2.;

    // https://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients
    #pragma omp parallel for
    for(int n = 0; n < loopLength; n++){

        arma::cx_vec coefsK, coefsK2, coefsKQ, coefsK2Q;

        int ii = loopLength - 1 - n;
        int m  = floor((sqrt(8*ii + 1) - 1)/2);
        int i = basisDimBSE - 1 - m;
        int j = basisDimBSE - 1 - ii + m*(m+1)/2;
    
        double k_index = basisStates(i, 2);
        int v = bandToIndex[basisStates(i, 0)];
        int c = bandToIndex[basisStates(i, 1)];
        int kQ_index = k_index;

        double k2_index = basisStates(j, 2);
        int v2 = bandToIndex[basisStates(j, 0)];
        int c2 = bandToIndex[basisStates(j, 1)];
        int k2Q_index = k2_index;


        // Using the atomic gauge
        if(gauge == "atomic"){
            coefsK = latticeToAtomicGauge(eigvecKStack.slice(k_index).col(v), kpoints.row(k_index));
            coefsKQ = latticeToAtomicGauge(eigvecKQStack.slice(kQ_index).col(c), kpoints.row(kQ_index));
            coefsK2 = latticeToAtomicGauge(eigvecKStack.slice(k2_index).col(v2), kpoints.row(k2_index));
            coefsK2Q = latticeToAtomicGauge(eigvecKQStack.slice(k2Q_index).col(c2), kpoints.row(k2Q_index));
        }
        else{
            coefsK = eigvecKStack.slice(k_index).col(v);
            coefsKQ = eigvecKQStack.slice(kQ_index).col(c);
            coefsK2 = eigvecKStack.slice(k2_index).col(v2);
            coefsK2Q = eigvecKQStack.slice(k2Q_index).col(c2);
        }

        std::complex<double> D, X;
        if (mode == "realspace"){
            int effective_k_index = findEquivalentPointBZ(kpoints.row(k2_index) - kpoints.row(k_index), ncell);
            arma::cx_mat motifFT = ftMotifStack.slice(effective_k_index);
            D = exactInteractionTermMFT(coefsKQ, coefsK2, coefsK2Q, coefsK, motifFT);
            X = 0;

        }
        else if (mode == "reciprocalspace"){
            arma::rowvec k = kpoints.row(k_index);
            arma::rowvec k2 = kpoints.row(k2_index);
            D = interactionTermFT(coefsK, coefsK2, coefsKQ, coefsK2Q, k, k2, k, k2, this->nReciprocalVectors);
            X = 0;
        }
        
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
Result GExciton::diagonalize(std::string method, int nstates){
    std::cout << "Solving BSE with ";
    arma::vec eigval;
    arma::cx_mat eigvec;

    if (method == "diag"){
        std::cout << "exact diagonalization... " << std::flush;
        arma::eig_sym(eigval, eigvec, HBS);
    }
    else if (method == "davidson"){
        std::cout << "Davidson method... " << std::flush;
        davidson_method(eigval, eigvec, HBS, nstates);
    }
    
    std::cout << "Done" << std::endl;
    Result results = Result(*this, eigval, eigvec);

    return results;
}

// ------------- Symmetries -------------
arma::mat GExciton::C3ExcitonBasisRep(){
    calculateInverseReciprocalMatrix();
    arma::mat C3 = arma::zeros(excitonbasisdim, excitonbasisdim);
    int nbandCombinations = valenceBands.n_elem * conductionBands.n_elem;
    if(kpoints.empty()){
        throw std::logic_error("kpoints must be initialized first");
    }
    for(unsigned int i = 0; i < kpoints.n_rows; i++){
        arma::rowvec rotatedK = rotateC3(kpoints.row(i));
        int kIndex = findEquivalentPointBZ(rotatedK, ncell);

        for(int j = 0; j < nbandCombinations; j++){
            C3(kIndex*nbandCombinations + j, i*nbandCombinations + j) = 1;
        }
    }

    return C3;
}

// ------------- Routines to compute Fermi Golden Rule -------------

/* Method to compute density of states associated to non-interacting electron-hole pairs.
Considers only the bands defined as the basis for excitons. */
double GExciton::pairDensityOfStates(double energy, double delta) const{
    
    double dos = 0;
    for(int v = 0; v < (int)valenceBands.n_elem; v++){
        for(int c = 0; c < (int)conductionBands.n_elem; c++){
            for(int i = 0; i < nk; i++){

                arma::uword vband = bandToIndex.at(valenceBands(v)); // Unsigned integer 
                arma::uword cband = bandToIndex.at(conductionBands(c));

                double stateEnergy = eigvalKStack.col(i)(cband) - eigvalKStack.col(i)(vband);
                dos += -PI*imag(rGreenF(energy, delta, stateEnergy));
            };
        }
    }

    return dos;
}


/* Routine to compute the non-interacting electron-hole edge pair 
associated to a given energy. To do so, we run a search algorithm to 
find which k value matches the given energy.
Input: double energy, vec kpoints (provides the search grid)
vec gapEnergy: vector of gap energies associated to each k
Output: vec of coefficients in e-h pair basis associated to desired
pair */
cx_vec GExciton::ehPairCoefs(double energy, const vec& gapEnergy, std::string side){

    cx_vec coefs = arma::zeros<cx_vec>(nk);
    int closestKindex = -1;
    double eDiff, prevDiff;
    for(int n = 1; n < nk/2; n++){

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
    double dispersion = PI/(16*a);
    if(side == "left"){
        coefs(closestKindex) = 1.;
    }
    else if(side == "right"){
        coefs(nk - 1 - closestKindex) = 1.;
    }

    cout << "Energy gap (-k): " << gapEnergy(closestKindex) << endl;
    cout << "Energy gap (k): " << gapEnergy(nk - 1 - closestKindex) << endl;
    this->pairEnergy = gapEnergy(closestKindex);

    return coefs;
};


/* Method to compute the transition rate from one exciton to a general non-interacting electron-hole pair.
It checks that the energy of the final state is close to the specified energy. Otherwise, it throws an error. 
It assumes that both exciton systems have been initialized already. */
double GExciton::fermiGoldenRule(const GExciton& targetExciton, const arma::cx_vec& initialState, const arma::cx_vec& finalState, double energy){

    double transitionRate = 0;
    arma::imat initialBasis = basisStates;
    arma::imat finalBasis = targetExciton.basisStates;
    cx_mat W = arma::zeros<cx_mat>(finalBasis.n_rows, initialBasis.n_rows);

    // -------- Main loop (W initialization) --------
    #pragma omp parallel for schedule(static, 1) collapse(2)
    for (int i = 0; i < finalBasis.n_rows; i++){
        for (int j = 0; j < initialBasis.n_rows; j++){

            arma::cx_vec coefsK, coefsK2, coefsKQ, coefsK2Q;

            int vf = targetExciton.bandToIndex.at(finalBasis(i, 0));
            int cf = targetExciton.bandToIndex.at(finalBasis(i, 1));
            double kf_index = finalBasis(i, 2);
            
            int vi = bandToIndex[initialBasis(j, 0)];
            int ci = bandToIndex[initialBasis(j, 1)];
            double ki_index = initialBasis(j, 2);

            // Using the atomic gauge
            if(gauge == "atomic"){
                coefsK = latticeToAtomicGauge(targetExciton.eigvecKStack.slice(kf_index).col(vf), kpoints.row(kf_index));
                coefsKQ = latticeToAtomicGauge(targetExciton.eigvecKQStack.slice(kf_index).col(cf), kpoints.row(kf_index));
                coefsK2 = latticeToAtomicGauge(eigvecKStack.slice(ki_index).col(vi), kpoints.row(ki_index));
                coefsK2Q = latticeToAtomicGauge(eigvecKQStack.slice(ki_index).col(ci), kpoints.row(ki_index));
            }
            else{
                coefsK = targetExciton.eigvecKStack.slice(kf_index).col(vf);
                coefsKQ = targetExciton.eigvecKQStack.slice(kf_index).col(cf);
                coefsK2 = eigvecKStack.slice(ki_index).col(vi);
                coefsK2Q = eigvecKQStack.slice(ki_index).col(ci);
            }

            std::complex<double> D, X;
            if (mode == "realspace"){
                int effective_k_index = findEquivalentPointBZ(kpoints.row(ki_index) - kpoints.row(kf_index), ncell);
                arma::cx_mat motifFT = ftMotifStack.slice(effective_k_index);
                D = exactInteractionTermMFT(coefsKQ, coefsK2, coefsK2Q, coefsK, motifFT);
                X = 0;

            }
            else if (mode == "reciprocalspace"){
                arma::rowvec k = kpoints.row(ki_index);
                arma::rowvec k2 = kpoints.row(kf_index);
                D = interactionTermFT(coefsK, coefsK2, coefsKQ, coefsK2Q, k, k2, k, k2, this->nReciprocalVectors);
                X = 0;
            }
            
            W(i, j) = - (D - X);                
        };
    };

    double delta = 2.4/(2*ncell); // Adjust delta depending on number of k points
    double rho = targetExciton.pairDensityOfStates(pairEnergy, delta);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E-16; // Units are eV*s

    transitionRate = 2*PI*pow(abs(arma::cdot(finalState, W*initialState)),2)*rho/hbar;

    return transitionRate;

}


void GExciton::printInformation(){
    cout << std::left << std::setw(30) << "Number of cells: " << ncell << endl;
    cout << std::left << std::setw(30) << "Valence bands:";
    for (int i = 0; i < valenceBands.n_elem; i++){
        cout << valenceBands(i) << "\t";
    }
    cout << endl;

    cout << std::left << std::setw(30) << "Conduction bands: ";
    for (int i = 0; i < conductionBands.n_elem; i++){
        cout << conductionBands(i) << "\t";
    }
    cout << endl;

    cout << std::left << std::setw(30) << "Gauge used: " << gauge << "\n" << endl;
}