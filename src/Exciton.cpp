#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iomanip>

#include "xatu/System.hpp"
#include "xatu/Exciton.hpp"
#include "xatu/utils.hpp"
#include "xatu/davidson.hpp"

using namespace arma;
using namespace std::chrono;

namespace xatu {

/**
 * Method to set the values of the attributes of an exciton object.
 * @param ncell Number of unit cells per axis.
 * @param bands Vector with the indices of the bands that form the exciton.
 * @param parameters Dielectric constants and screening length.
 * @param Q Center-of-mass momentum of the exciton.
 * @return void 
 */
void Exciton::initializeExcitonAttributes(int ncell, const arma::ivec& bands, 
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

/**
 * Method to set the attributes of an exciton object from a ExcitonConfiguration object.
 * @details Overload of the method to use a configuration object. Based on the parametric method.
 * @param cfg ExcitonConfiguration object from parsed file.
 * @return void
 */
void Exciton::initializeExcitonAttributes(const ExcitonConfiguration& cfg){

    int ncell        = cfg.excitonInfo.ncell;
    int nbands       = cfg.excitonInfo.nbands;
    arma::ivec bands = cfg.excitonInfo.bands;
    arma::rowvec parameters = {cfg.excitonInfo.eps(0), cfg.excitonInfo.eps(1), cfg.excitonInfo.eps(2)};
    arma::rowvec Q   = cfg.excitonInfo.Q;

    if (bands.empty()){
        bands = arma::regspace<arma::ivec>(- nbands + 1, nbands);
    }

    initializeExcitonAttributes(ncell, bands, parameters, Q);

    if(cfg.excitonInfo.submeshFactor != 1){
        this->totalCells_ = pow(ncell * cfg.excitonInfo.submeshFactor, ndim);
    }
    this->exchange = cfg.excitonInfo.exchange;

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
    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));

    // Set flags
    this->exchange = cfg.excitonInfo.exchange;
    this->scissor_ = cfg.excitonInfo.scissor;
    this->mode_    = cfg.excitonInfo.mode;
    this->nReciprocalVectors_ = cfg.excitonInfo.nReciprocalVectors;
}

/**
 * Exciton constructor from a SystemConfiguration object and a vector with the bands that form
 * the exciton, as well as the other parameters.
 * @param config SystemConfiguration object from config file.
 * @param ncell Number of unit cells along each axis.
 * @param bands Vector with the indices of the bands that form the exciton.
 * @param parameters Vector with dielectric constants and screening length.
 * @param Q Center-of-mass momentum.
 */
Exciton::Exciton(const SystemConfiguration& config, int ncell, const arma::ivec& bands, 
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
    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));
};

/**
 * Exciton constructor from a SystemConfiguration object. One specifies the number of valence and conduction
 * bands, as well as the other parameters.
 * @param config SystemConfiguration object from config file.
 * @param ncell Number of unit cells along each axis.
 * @param nbands Number of bands (same number for both valence and conduction) that form the exciton.
 * @param nrmbands Number of bands to be removed with respect to the Fermi level. 
 * @param parameters Vector with dielectric constants and screening length.
 * @param Q Center-of-mass momentum.
 */
Exciton::Exciton(const SystemConfiguration& config, int ncell, int nbands, int nrmbands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          Exciton(config, ncell, {}, parameters, Q) {
    
    if (2*nbands > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }
    this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands - nrmbands + 1, 
                                                     fermiLevel - nrmbands);
    this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, 
                                                        fermiLevel + nbands + nrmbands);
    this->bands_ = arma::join_cols(valenceBands, conductionBands) - fermiLevel;
    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));
};


/**
 * Exciton constructor from SystemConfiguration and ExcitonConfiguration.
 * @param config SystemConfiguration object.
 * @param excitonConfig ExcitonConfiguration object.
 */ 
Exciton::Exciton(const SystemConfiguration& config, const ExcitonConfiguration& excitonConfig) : System(config){
    initializeExcitonAttributes(excitonConfig);
}

Exciton::Exciton(const System& system, int ncell, const arma::ivec& bands, 
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
    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));
};

/**
 * Exciton constructor from an already initialized System object, and all the exciton parameters.
 * @param system System object where excitons are computed.
 * @param ncell Number of unit cells along each axis.
 * @param nbands Number of bands (same for valence and conduction) that form the exciton.
 * @param nrmbands Number of bands to be removed with respect to the Fermi level. 
 * @param parameters Dielectric constant and screening length.
 * @param Q Center-of-mass momentum of the exciton.
 */
Exciton::Exciton(const System& system, int ncell, int nbands, int nrmbands, 
                  const arma::rowvec& parameters, const arma::rowvec& Q) : 
          Exciton(system, ncell, {}, parameters, Q) {
    
    if (2*nbands > basisdim){
        cout << "Error: Number of bands cannot be higher than actual material bands" << endl;
        exit(1);
    }

    this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands - nrmbands + 1, 
                                                     fermiLevel - nrmbands);
    this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, 
                                                        fermiLevel + nbands + nrmbands);
    this->bands_ = arma::join_cols(valenceBands, conductionBands) - fermiLevel;
    this->bandList_ = arma::conv_to<arma::uvec>::from(arma::join_cols(valenceBands, conductionBands));
};

/** 
 * Exciton destructor.
 * @details Used mainly for debugging; the message should be removed at some point.
 */
Exciton::~Exciton(){};


/* ------------------------------ Setters ------------------------------ */

/**
 * Sets the number of unit cells along each axis.
 * @param ncell Number of unit cells per axis.
 * @return void
 */
void Exciton::setUnitCells(int ncell){
    if(ncell > 0){
        ncell_ = ncell;
    }
    else{
        std::cout << "ncell must be a positive number" << std::endl;
    }
}

/**
 * Sets the bands involved in the exciton calculation from a vector.
 * @param bands Vector of integers corresponding to the indices of the bands.
 * @return void 
 */
void Exciton::setBands(const arma::ivec& bands){
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

/**
 * Sets the bands involved in the exciton calculation specifying the number of bands
 * above and below the Fermi level.
 * @param nbands Number of valence (conduction) bands used.
 * @param nrmbands Number of valence (conduction) bands removed from calculation.
 */
void Exciton::setBands(int nbands, int nrmbands){
    if(nbands > 0 && nrmbands > 0){
        this->valenceBands_ = arma::regspace<arma::ivec>(fermiLevel - nbands + 1, fermiLevel - nrmbands);
        this->conductionBands_ = arma::regspace<arma::ivec>(fermiLevel + 1 + nrmbands, fermiLevel + nbands);
        this->bands_ = arma::join_rows(valenceBands, conductionBands);
    }
    else{
        std::cout << "Included bands and removed bands must be positive numbers" << std::endl;
    }
}

/**
 * Sets the center-of-mass momentum of the exciton.
 * @param Q Momentum vector.
 * @return void 
 */
void Exciton::setQ(const arma::rowvec& Q){
    if(Q.n_elem == 3){
        Q_ = Q;
    }
    else{
        std::cout << "Q vector must be 3d" << std::endl;
    }
    
}

/**
 * Method to set the parameters of the Keldysh potential, namely the environmental
 * dielectric constants and the effective screening length.
 * @param parameters Vector with the three parameters, '{eps_m, eps_s, r0}'. 
 * @return void
 */
void Exciton::setParameters(const arma::rowvec& parameters){
    if(parameters.n_elem == 3){
        eps_m_ = parameters(0);
        eps_s_ = parameters(1);
        r0_    = parameters(2);
    }
    else{
        std::cout << "parameters array must be 3d (eps_m, eps_s, r0)" << std::endl;
    }
}

/**
 * Sets the parameters of the Keldysh potential.
 * @param eps_m Dielectric constant of embedding medium.
 * @param eps_s Dielectric constant of substrate.
 * @param r0 Effective screeening length.
 * @return void 
 */
void Exciton::setParameters(double eps_m, double eps_s, double r0){
    // TODO: Introduce additional comprobations regarding value of parameters (positive)
    eps_m_ = eps_m;
    eps_s_ = eps_s;
    r0_    = r0;
}

/**
 * Sets the cutoff over unit cells used in the calculation of the lattice Fourier transform
 * for the interactions.
 * @param cutoff Number of unit cells to consider.
 * @return void 
 */
void Exciton::setCutoff(double cutoff){
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

/**
 * Sets the gauge used for the Bloch basis, either 'lattice' or 'atomic'.
 * @param gauge Gause to be used, default to 'lattice'.
 * @return void
 */
void Exciton::setGauge(std::string gauge){
    if(gauge != "lattice" && gauge != "atomic"){
        throw std::invalid_argument("setGauge(): gauge must be either lattice or atomic");
    }
    this->gauge_ = gauge;
}

/**
 * Sets the type of calculation used to obtain the exciton spectrum. It can be 'realspace' (default) 
 * or 'reciprocalspace'.
 * @param mode Calculation model.
 * @return void 
 */
void Exciton::setMode(std::string mode){
    if(mode != "realspace" && mode != "reciprocalspace"){
        throw std::invalid_argument("setMode(): mode must be either realspace or reciprocalspace");
    }
    this->mode_ = mode;
}

/**
 * Sets the number of reciprocal vectors to use if the exciton calculation is set to 'reciprocalspace'.
 * @param nReciprocalVector Number of reciprocal vectors to sum over.
 * @return void 
 */
void Exciton::setReciprocalVectors(int nReciprocalVectors){
    if(nReciprocalVectors < 0){
        throw std::invalid_argument("setReciprocalVectors(): given number must be positive");
    }
    this->nReciprocalVectors_ = nReciprocalVectors;
}

/**
 * Sets the value of the scissor cut of change the gap of the system.
 * @param shift Value of scissor cut (in eV). Can be positive or negative.
 * @return void 
 */
void Exciton::setScissor(double shift){
    this->scissor_ = shift;
}

/**
 * To toggle on or off the exchange term in the interaction matrix elements.
 * @param exchange Either true of false
 * @return void
*/
void Exciton::setExchange(bool exchange){
    this->exchange = exchange;
}


/*---------------------------------------- Potentials ----------------------------------------*/

/** 
 * Purpose: Compute Struve function H0(x).
 * Source: http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mstvh0_cpp.txt 
 * @param X x --- Argument of H0(x) ( x ò 0 )
 * @param SH0 SH0 --- H0(x). The return value is written to the direction of the pointer.
*/
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


/** 
 * Calculate value of interaction potential (Keldysh). Units are eV.
 * @details If the distance is zero, then the interaction is renormalized to be V(a) since
 * V(0) is infinite, where a is the lattice parameter. Also, for r > cutoff the interaction is taken to be zero.
 * @param r Distance at which we evaluate the potential.
 * @return Value of Keldysh potential, V(r).
 */
double Exciton::potential(double r){
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

/**
 * Calculate lattice Fourier transform of Keldsyh potential.
 * @param k kpoint where we compute the FT.
 * @param cells Matrix with position of unit cells over which we sum to obtain the lattice FT.
 * @param useApproximation Boolean to toggle simplified FT calculation (one summation instead of two over unit cells).
 * @return Lattice Fourier transform of the potential, lFT[V](k).
 */
std::complex<double> Exciton::fourierTransform(arma::rowvec k, const arma::mat& cells, bool useApproximation){
    std::complex<double> imag(0, 1);
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

/**
 * Evaluates the Fourier transform of the Keldysh potential, which is an analytical expression.
 * @param q kpoint where we evaluate the FT.
 * @return Fourier transform of the potential at q, FT[V](q).
 */
double Exciton::analyticFourierTransform(arma::rowvec q){
    double radius = cutoff*arma::norm(reciprocalLattice.row(0));
    double potential = 0;
    double eps_bar = (eps_m + eps_s)/2;
    double eps = arma::norm(reciprocalLattice.row(0))/totalCells;

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

/**
 * Routine to compute the lattice Fourier transform with the potential displaced by some
 * vectors of the motif.
 * @param fAtomIndex Index of first atom of the motif.
 * @param sAtomIndex Index of second atom of the motif.
 * @param k kpoint where we evaluate the FT.
 * @param cells Matrix with the unit cells over which we sum to compute the lattice FT.
 * @return Motif lattice Fourier transform of the Keldysh potential at k.
 */
std::complex<double> Exciton::motifFourierTransform(int fAtomIndex, int sAtomIndex, const arma::rowvec& k, const arma::mat& cells){

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

/**
 * Method to extend the motif Fourier transform matrix to match the dimension of the
 * one-particle basis. 
 * @param motifFT Matrix storing the motif Fourier transform to be extended.
 * @return Extended matrix.
 */
arma::cx_mat Exciton::extendMotifFT(const arma::cx_mat& motifFT){
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

/** 
 * Real space implementation of interaction term, valid for both direct and exchange.
 * To compute the direct term, the expected order is (ck,v'k',c'k',vk).
 * For the exchange term, the order is (ck,v'k',vk,c'k').
 * @param coefsK1 First eigenstate vector.
 * @param coefsK2 Second eigenstate vector.
 * @param coefsK3 Third eigenstate vector.
 * @param coefsK4 Fourth eigenstate vector.
 * @param motifFT Motif Fourier transform.
 * @return Interaction term.
 */
std::complex<double> Exciton::exactInteractionTermMFT(const arma::cx_vec& coefsK1, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsK3, 
                                     const arma::cx_vec& coefsK4,
                                     const arma::cx_mat& motifFT){
    
    arma::cx_vec firstCoefArray = arma::conj(coefsK1) % coefsK3;
    arma::cx_vec secondCoefArray = arma::conj(coefsK2) % coefsK4;
    std::complex<double> term = arma::dot(firstCoefArray, extendMotifFT(motifFT) * secondCoefArray);
        
    return term;
};

/**
 * Reciprocal space implementation of interaction term, valid for both direct and exchange.
 * @param coefsK Vector of eigenstate |v,k>.
 * @param coefsK2 Vector of eigenstate |v',k'>.
 * @param coefsKQ Vector of eigenstate |c,k+Q>.
 * @param coefsK2Q Vector of eigenstate |c',k'+Q>.
 * @param k kpoint corresponding to k.
 * @param k2 kpoint corresponding to k'.
 * @param kQ kpoint corresponding to k + Q.
 * @param k2Q kpoint corresponding to k' + Q.
 * @return Interaction term.
 */
std::complex<double> Exciton::interactionTermFT(const arma::cx_vec& coefsK, 
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

        Ic = blochCoherenceFactor(coefsKQ, coefsK2Q, kQ, k2Q, G);
        Iv = blochCoherenceFactor(coefsK, coefsK2, k, k2, G);

        term += Ic*conj(Iv)*analyticFourierTransform(k - k2 + G);
    }

    return term;
};

/**
 * Calculation of Bloch coherence factors, required to compute the interaction terms in reciprocal space.
 * @param coefs1 Vector of eigenstate |n,k>.
 * @param coefs2 Vector of eigenstate |n',k'>.
 * @param k1 kpoint k.
 * @param k2 kpoint k'.
 * @param G Reciprocal lattice vector used to compute the coherence factor.
 * @return Bloch coherence factor I evaluated at G for states |nk>, |n'k'>.
 */
std::complex<double> Exciton::blochCoherenceFactor(const arma::cx_vec& coefs1, const arma::cx_vec& coefs2,
                                                    const arma::rowvec& k1, const arma::rowvec& k2,
                                                    const arma::rowvec& G){

    std::complex<double> imag(0, 1);
    arma::cx_vec coefs = arma::conj(coefs1) % coefs2;
    arma::cx_vec phases = arma::ones<arma::cx_vec>(basisdim);
    for(int i = 0; i < natoms; i++){
        int species = motif.row(i)(3);
        arma::rowvec atomPosition = motif.row(i).subvec(0, 2);
        phases.subvec(i*orbitals(species), (i+1)*orbitals(species) - 1) *= 
        exp(imag*arma::dot(k1 - k2 + G, atomPosition));
    }

    std::complex<double> factor = arma::dot(coefs, phases);

    return factor;
}


/*------------------------------------ Electron-hole pair basis ------------------------------------*/

/**
 * Initialise basis to be used in the construction of the BSE matrix.
 * @param conductionBands Conduction bands that will populate the electrons of the exciton.
 * @param valenceBands Valence bands to be populated by holes of the exciton.
 * @return Matrix where each row denotes an electron-hole pair, '{v, c, k}'.
 */
arma::imat Exciton::createBasis(const arma::ivec& conductionBands, 
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

    return states;
};

/**
 * Overload of createBasis method to work with class attributes instead of given ones.
 * @return void.
 */
void Exciton::initializeBasis(){
    this->basisStates_ = createBasis(conductionBands, valenceBands);
    //createSOCBasis();
};

/**
 * Method to generate a basis which is a subset of the basis considered for the
 * exciton. Its main purpose is to allow computation of Fermi golden rule between
 * two specified subbasis. 
 * @param bands Band subset of the originally specified for the exciton.
 * @return Matrix with the states corresponding to the specified subset.
 */
arma::imat Exciton::specifyBasisSubset(const arma::ivec& bands){

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


/**
 * Compute the basis elements for the spinful exciton problem. Reorders basis
 * in blocks of defined spin (so that they are diagonal for later calculation of eigenstates of BSE).
 * @return void
 */
void Exciton::useSpinfulBasis(){

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

/** 
 * Routine to calculate the coefficients corresponding to wavefunctions in the atomic gauge.
 * @param coefs Lattice gauge state coefficients on which we perform the gauge transformation.
 * @param k kpoint required to perform the transformation.
 * @return Atomic gauge state coefficients.
 */
arma::cx_vec Exciton::latticeToAtomicGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

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

/**
 * Method to transform the single-particle state coefficients from atomic to lattice gauge.
 * @param coefs Atomic gauge coefficients.
 * @param k kpoint used in the transformation.
 * @return Lattice gauge coefficients.
 */
arma::cx_vec Exciton::atomicToLatticeGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

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

/**
 * Criterium to fix the phase of the single-particle eigenstates after diagonalization.
 * @details The prescription we take here is to impose that the sum of all the coefficients is real.
 * @return Fixed coefficients. 
 */
arma::cx_mat Exciton::fixGlobalPhase(arma::cx_mat& coefs){

    arma::cx_rowvec sums = arma::sum(coefs);
    std::complex<double> imag(0, 1);
    for(int j = 0; j < sums.n_elem; j++){
        double phase = arg(sums(j));
        coefs.col(j) *= exp(-imag*phase);
    }

    return coefs;
}

/**
 * Creates a dictionary that maps bands to indices for storage.
 * @return void
 */
void Exciton::generateBandDictionary(){

    std::map<int, int> bandToIndex;
    for(int i = 0; i < bandList.n_elem; i++){
        bandToIndex[bandList(i)] = i;
    };

    this->bandToIndex = bandToIndex;
};

/**
 * Method to compute the motif FT matrix at a given k vector.
 * @param k k vector where we compute the motif FT.
 * @param cells Matrix of unit cells over which the motif FT is computed.
 * @return void
 */
arma::cx_mat Exciton::motifFTMatrix(const arma::rowvec& k, const arma::mat& cells){
    // Uses hermiticity of V
    arma::cx_mat motifFT = arma::zeros<arma::cx_mat>(natoms, natoms);

    for(int fAtomIndex = 0; fAtomIndex < natoms; fAtomIndex++){
        for(int sAtomIndex = fAtomIndex; sAtomIndex < natoms; sAtomIndex++){
            motifFT(fAtomIndex, sAtomIndex) = motifFourierTransform(fAtomIndex, sAtomIndex, k, cells);
            motifFT(sAtomIndex, fAtomIndex) = conj(motifFT(fAtomIndex, sAtomIndex));
        }   
    }

    return motifFT;
}

/**
 * Method to initialize the motif Fourier transform for all possible motif combination 
 * at a given kpoint.
 * @param i Index of kpoint.
 * @param cells Matrix of unit cells over which the motif FT is computed.
 * @return void
 */
void Exciton::initializeMotifFT(int i, const arma::mat& cells){
    ftMotifStack.slice(i) = motifFTMatrix(meshBZ_.row(i), cells);
}


/**
 * Main method to compute all the relevant single-particle quantities (bands, eigenstates and fourier transforms),
 * to compute the Bethe-Salpeter equation.
 * @details It precomputes and saves the relevant data in the heap for later computations.
 * @param triangular Boolean to specify whether the Hamiltonian matrices are triangular (default = false).
 * @return void
 */ 
void Exciton::initializeResultsH0(bool triangular){

    int nTotalBands = bandList.n_elem;
    double radius = arma::norm(bravaisLattice.row(0)) * cutoff_;
    arma::mat cells = truncateSupercell(ncell, radius);

    this->eigvecKStack_  = arma::cx_cube(basisdim, nTotalBands, nk);
    this->eigvecKQStack_ = arma::cx_cube(basisdim, nTotalBands, nk);
    this->eigvalKStack_  = arma::mat(nTotalBands, nk);
    this->eigvalKQStack_ = arma::mat(nTotalBands, nk);
    this->ftMotifStack   = arma::cx_cube(natoms, natoms, meshBZ_.n_rows);
    this->ftMotifQ       = arma::cx_mat(natoms, natoms);

    vec auxEigVal(basisdim);
    cx_mat auxEigvec(basisdim, basisdim);
    cx_mat h;

    // Progress bar variables
    int step = 1;
	int displayNext = step;
	int percent = 0;

    calculateInverseReciprocalMatrix();
    std::complex<double> imag(0, 1);

    std::cout << "Diagonalizing H0 for all k points... " << std::flush;
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
        
    };
    std::cout << "Done" << std::endl;

    if(this->mode == "realspace"){
        std::cout << "Computing lattice Fourier transform..." << std::endl;
        for (unsigned int i = 0; i < meshBZ_.n_rows; i++){
            // BIGGEST BOTTLENECK OF THE CODE
            initializeMotifFT(i, cells);     

		    percent = (100 * (i + 1)) / meshBZ_.n_rows ;
		    if (percent >= displayNext){
                cout << "\r" << "[" << std::string(percent / 5, '|') << std::string(100 / 5 - percent / 5, ' ') << "]";
                cout << percent << "%";
                std::cout.flush();
                displayNext += step;
            }
        }
        std::cout << "\nDone" << std::endl;
    }

    if(this->exchange){
        this->ftMotifQ = motifFTMatrix(this->Q, cells);
    }
};

/**
 * Routine to initialize the required variables to construct the Bethe-Salpeter Hamiltonian.
 * @param triangular Boolean to specify whether the single-particle Hamiltonian matrices are triangular.
 * @return void.
 */
void Exciton::initializeHamiltonian(bool triangular){

    if(bands.empty()){
        throw std::invalid_argument("Error: Exciton object must have some bands");
    }
    if(nk == 0){
        throw std::invalid_argument("Error: BZ mesh must be initialized first");
    }

    this->excitonbasisdim_ = nk*valenceBands.n_elem*conductionBands.n_elem;

    std::cout << "Initializing basis for BSE... " << std::flush;
    initializeBasis();
    //useSpinfulBasis();
    generateBandDictionary();

    initializeResultsH0(triangular);
}


/**
 * Initialize BSE hamiltonian matrix and kinetic matrix.
 * @details Instead of calculating the energies and coeficients dinamically, which
 * is too expensive, instead we first calculate those for each k, save them
 * in the heap, and then call them consecutively as we build the matrix.
 * Analogously, we calculate the Fourier transform of the potential beforehand,
 * saving it in the stack so that it can be later called in the matrix element
 * calculation.
 * Also note that this routine involves a omp parallelization when building the matrix.
 * @param basis Subset of the exciton basis to build the BSE. If none, defaults to
 * the complete or original basis.
 * @return void
 */
void Exciton::BShamiltonian(const arma::imat& basis){

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
    long int loopLength = basisDimBSE*(basisDimBSE + 1)/2.;

    // https://stackoverflow.com/questions/242711/algorithm-for-index-numbers-of-triangular-matrix-coefficients
    #pragma omp parallel for
    for(long int n = 0; n < loopLength; n++){

        arma::cx_vec coefsK, coefsK2, coefsKQ, coefsK2Q;

        long int ii = loopLength - 1 - n;
        long int m  = floor((sqrt(8*ii + 1) - 1)/2);
        long int i = basisDimBSE - 1 - m;
        long int j = basisDimBSE - 1 - ii + m*(m+1)/2;
    
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

        std::complex<double> D, X = 0.0;
        if (mode == "realspace"){
            int effective_k_index = findEquivalentPointBZ(kpoints.row(k2_index) - kpoints.row(k_index), ncell);
            arma::cx_mat motifFT = ftMotifStack.slice(effective_k_index);
            D = exactInteractionTermMFT(coefsKQ, coefsK2, coefsK2Q, coefsK, motifFT);
            if(this->exchange){
                X = exactInteractionTermMFT(coefsKQ, coefsK2, coefsK, coefsK2Q, this->ftMotifQ);
            }            
        }
        else if (mode == "reciprocalspace"){
            arma::rowvec k = kpoints.row(k_index);
            arma::rowvec k2 = kpoints.row(k2_index);
            D = interactionTermFT(coefsK, coefsK2, coefsKQ, coefsK2Q, k, k2, k, k2, this->nReciprocalVectors);
            if(this->exchange){
                X = interactionTermFT(coefsK2Q, coefsK2, coefsKQ, coefsK, k2 + Q, k2, k + Q, k, this->nReciprocalVectors);
            }
        }
        
        if (i == j){
            HBS_(i, j) = (this->scissor + 
                          eigvalKQStack.col(kQ_index)(c) - eigvalKStack.col(k_index)(v))/2. 
                          - (D - X)/2.;
            HK_(i, j) = eigvalKQStack(c, kQ_index) - eigvalKStack(v, k_index);
            
        }
        else{
            HBS_(i, j) =  - (D - X);
        };
    }
       
    HBS_ = HBS + HBS.t();
    std::cout << "Done" << std::endl;
};

/**
 * Routine to diagonalize the BSE and return a Result object.
 * @param method Method to diagonalize the BSE, either 'diag' (standard diagonalization) 
 * 'davidson' (iterative diagonalization) or 'sparse' (Lanczos).
 * @param nstates Number of states to be stored from the diagonalization.
 * @return Result object storing the exciton energies and states.
 */ 
Result Exciton::diagonalize(std::string method, int nstates){
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
    else if (method == "sparse"){
        std::cout << "Lanczos method..." << std::flush;

        arma::cx_vec cx_eigval;
        arma::eigs_gen(cx_eigval, eigvec, arma::sp_cx_mat(HBS), nstates, "sr");
        eigval = arma::sort(real(cx_eigval));
    }
    
    std::cout << "Done" << std::endl;
    Result results = Result(*this, eigval, eigvec);

    return results;
}

// ------------- Symmetries -------------
/**
 * Method to obtain the C3 rotation operator in the basis of electron-hole pairs of the exciton.
 * @return Matrix representation of C3. 
 */
arma::mat Exciton::C3ExcitonBasisRep(){
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

/**
 * Method to compute density of states associated to non-interacting electron-hole pairs.
 * Considers only the bands defined as the basis for excitons.
 * @param energy Energy at which we evaluate the pair DoS.
 * @param delta Broadening used to smooth the DoS.
 * @return DoS at E.
 */
double Exciton::pairDensityOfStates(double energy, double delta) const {
    
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
    dos /= (a*nk);

    return dos;
}


/** 
 * Routine to compute and write to a file the density of states of non-interacting e-h pairs.
 * @param file Pointer to file.
 * @param delta DoS broadening.
 * @param n Number of points on which we evaluate the DoS.
 * @return void
 */
void Exciton::writePairDOS(FILE* file, double delta, int n){

    double eMin = eigvalKStack.min();
    double eMax = eigvalKQStack.max();
    arma::vec energies = arma::linspace(0, (eMax - eMin)*1.1, n);
    for (double energy : energies){
        double dos = pairDensityOfStates(energy, delta);
        fprintf(file, "%f\t%f\n", energy, dos);
    }
}


/**
 * Routine to compute the non-interacting electron-hole edge pair associated to a given energy.
 * @details We run a search algorithm to find which k value matches the given energy.
 * @param energy Energy at which we want the non-interacting e-h pair.
 * @param gapEnergy Values of the gap at all kpoints.
 * @param side Whether to obtain the pair at +k or -k.
 * @return Vec of coefficients in e-h pair basis associated to desired pair.
 */
cx_vec Exciton::ehPairCoefs(double energy, const vec& gapEnergy, std::string side){

    cx_vec coefs = arma::zeros<cx_vec>(nk);
    int closestKindex = -1;
    double eDiff;
    double currentEnergy = gapEnergy(0) - energy;

    for(int n = 1; n < nk/2; n++){
        
        eDiff = gapEnergy(n) - energy;
        if(abs(eDiff) < abs(currentEnergy)){
            closestKindex = n;
            currentEnergy = eDiff;
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

    return coefs;
};

/**
 * Method to compute the transition rate from one exciton to a general non-interacting electron-hole pair.
 * @param targetExciton Exciton object representing the final system.
 * @param initialState Exciton state (coefficients) from which the transition happens.
 * @param finalState Final exciton state in the transition.
 * @param energy Energy of the initial exciton state.
 * @return Transition rate from initialState to finalState.
 */ 
double Exciton::fermiGoldenRule(const Exciton& targetExciton, const arma::cx_vec& initialState, const arma::cx_vec& finalState, double energy){

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
                arma::rowvec k = kpoints.row(kf_index);
                arma::rowvec k2 = kpoints.row(ki_index);
                D = interactionTermFT(coefsK, coefsK2, coefsKQ, coefsK2Q, k, k2, k, k2, this->nReciprocalVectors);
                X = 0;
            }
            
            W(i, j) = - (D - X);                
        };
    };

    double delta = 2.4/(2*ncell); // Adjust delta depending on number of k points
    double rho = targetExciton.pairDensityOfStates(energy, delta);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E-16; // Units are eV*s

    transitionRate = 2*PI*std::norm(arma::cdot(finalState, W*initialState))*rho/hbar;

    return transitionRate;
}


/**
 * Method to identify a k point corresponding to a non-interacting electron-hole pair in the defined system
 * with the energy specified.
 * @param targetExciton Exciton object representing the general final states in the transition.
 * @param energy Energy of the initial state.
 * @param side Whether the transition takes place to an state with +k or -k.
 * @param increasing Used to specify whether the gap increases or decreases with k.
 * @return k vector of the equivalent electron-hole pair.
*/
arma::rowvec Exciton::findElectronHolePair(const Exciton& targetExciton, double energy, std::string side, bool increasing){

// First identify k edge e-h pair with same energy as exciton
    double n = 10; // Submeshing
    arma::rowvec min_k, max_k, kmin, kmax;
    if (side == "right"){
        min_k = kpoints.row(nk/2);
        max_k = -kpoints.row(0);
    }
    else if(side == "left"){
        max_k = kpoints.row(0);
        min_k = kpoints.row(nk/2 - 1);
    }
    
    arma::rowvec k;
    arma::cx_vec coefsK, coefsKQ, auxCoefsK, auxCoefsKQ;
    double threshold = 1E-8;
    arma::vec eigval;
    arma::cx_mat eigvec;
    int currentIndex;
    double currentEnergy = 0, vEnergy, cEnergy, gap, prevGap;
    double prevEnergy = currentEnergy;
    prevGap = 0;
    
    while(abs(currentEnergy - energy) > threshold){
        for(double i = 0; i <= n; i++){
            
            k = min_k * (1 - i/n) + max_k * i/n;
            targetExciton.solveBands(k, eigval, eigvec);

            eigval = eigval(targetExciton.bandList);
            vEnergy = eigval(0);

            eigvec = fixGlobalPhase(eigvec);
            eigvec = eigvec.cols(targetExciton.bandList);
            auxCoefsK = eigvec.col(0);

            if(arma::norm(Q) != 0){
                arma::rowvec kQ = k + Q;
                targetExciton.solveBands(kQ, eigval, eigvec);

                eigval = eigval(targetExciton.bandList);
                eigvec = fixGlobalPhase(eigvec);
                eigvec = eigvec.cols(targetExciton.bandList);
            }
            auxCoefsKQ = eigvec.col(1);
            cEnergy = eigval(1);

            gap = cEnergy - vEnergy;
            if (!increasing && (gap <= energy) && (prevGap > energy)){
                currentIndex = i;
                currentEnergy = gap;

                coefsK = auxCoefsK;
                coefsKQ = auxCoefsKQ;

                kmin = min_k * (1 - (currentIndex - 1)/n) + max_k * (currentIndex - 1)/n;
                kmax = min_k * (1 - (currentIndex + 1)/n) + max_k * (currentIndex + 1)/n;
            }
            if (increasing && (gap > energy) && (prevGap <= energy)){
                currentIndex = i;
                currentEnergy = gap;

                coefsK = auxCoefsK;
                coefsKQ = auxCoefsKQ;

                kmin = min_k * (1 - (currentIndex - 1)/n) + max_k * (currentIndex - 1)/n;
                kmax = min_k * (1 - (currentIndex + 1)/n) + max_k * (currentIndex + 1)/n;
            }
            prevGap = gap;

        }
        k = min_k * (1 - currentIndex/n) + max_k * currentIndex/n;
        min_k = kmin;
        max_k = kmax;
        arma::cout << "Current edge pair energy: " << currentEnergy << arma::endl;
        arma::cout << "Target energy: " << energy << "\n" << arma::endl;

        if (currentEnergy == prevEnergy){
            n += 1;
        }
        prevEnergy = currentEnergy;
    }

    arma::cout << "k: " << k << arma::endl;

    bool computeOccupations = false;
    if (computeOccupations){
        //////// Specific for Bi ribbon; must be deleted afterwards.
        int N = targetExciton.basisdim;
        double l_e_edge_occ = arma::norm(coefsKQ.subvec(0, 15));
        double r_e_edge_occ = arma::norm(coefsKQ.subvec(N - 16, N - 1));
        double l_h_edge_occ = arma::norm(coefsK.subvec(0, 15));
        double r_h_edge_occ = arma::norm(coefsK.subvec(N - 16, N - 1));

        std::cout << "left e occ.: " << l_e_edge_occ << "\nright e occ: " << r_e_edge_occ << std::endl;
        std::cout << "Total e occ.: " << std::sqrt(l_e_edge_occ*l_e_edge_occ + r_e_edge_occ*r_e_edge_occ) << arma::endl;
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "left h occ.: " << l_h_edge_occ << "\nright h occ: " << r_h_edge_occ << std::endl;
        std::cout << "Total h occ.: " << std::sqrt(l_h_edge_occ*l_h_edge_occ + r_h_edge_occ*r_h_edge_occ) << arma::endl;
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "Total e-h pair edge occu.: " << std::sqrt(l_e_edge_occ*l_e_edge_occ + r_e_edge_occ*r_e_edge_occ) + 
                    std::sqrt(l_h_edge_occ*l_h_edge_occ + r_h_edge_occ*r_h_edge_occ) << std::endl;
    }

    return k;
};

/**
 * Method to compute the transition to an edge e-h pair with the same energy (up to some error) as the bulk exciton.
 * @param targetExciton Exciton object representing the general final states in the transition.
 * @param initialState Exciton state from which the transition happens.
 * @param energy Energy of the initial state.
 * @param side Whether the transition takes place to an state with +k or -k.
 * @param increasing Used to specify whether the gap increases or decreases with k.
 * @return Transition rate from the initial exciton state to a non-interacting e-h pair.
 */
double Exciton::edgeFermiGoldenRule(const Exciton& targetExciton, const arma::cx_vec& initialState, double energy, std::string side, bool increasing){

    double transitionRate = 0;
    arma::imat initialBasis = basisStates;

    arma::rowvec k = findElectronHolePair(targetExciton, energy, side, increasing);

    arma::vec eigval;
    arma::cx_mat eigvec;
    arma::cx_vec coefsK, coefsKQ;

    targetExciton.solveBands(k, eigval, eigvec);

    eigvec = fixGlobalPhase(eigvec);
    eigvec = eigvec.cols(targetExciton.bandList);
    coefsK = eigvec.col(0);

    if(arma::norm(Q) != 0){
        arma::rowvec kQ = k + Q;
        targetExciton.solveBands(kQ, eigval, eigvec);

        eigvec = fixGlobalPhase(eigvec);
        eigvec = eigvec.cols(targetExciton.bandList);
    }
    coefsKQ = eigvec.col(1);
    
    // Now compute motif FT using k of edge pair
    
    double radius = arma::norm(bravaisLattice.row(0)) * cutoff_;
    arma::mat cells = truncateSupercell(ncell, radius);

    arma::cx_cube ftMotifStack = arma::cx_cube(natoms, natoms, meshBZ_.n_rows);
    
    for(int i = 0; i < nk; i++){
        for(int fAtomIndex = 0; fAtomIndex < natoms; fAtomIndex++){
            for(int sAtomIndex = fAtomIndex; sAtomIndex < natoms; sAtomIndex++){
                ftMotifStack(fAtomIndex, sAtomIndex, i) = 
                motifFourierTransform(fAtomIndex, sAtomIndex, meshBZ_.row(i) - k, cells);
                ftMotifStack(sAtomIndex, fAtomIndex, i) = conj(ftMotifStack(fAtomIndex, sAtomIndex, i));
            }   
        }
    }
    
    arma::cx_vec W = arma::zeros<arma::cx_vec>(initialBasis.n_rows);

    // -------- Main loop (W initialization) --------
    #pragma omp parallel for
    for (int i = 0; i < initialBasis.n_rows; i++){

        arma::cx_vec coefsK2, coefsK2Q;
        
        int vi = bandToIndex[initialBasis(i, 0)];
        int ci = bandToIndex[initialBasis(i, 1)];
        double ki_index = initialBasis(i, 2);

        // Using the atomic gauge
        if(gauge == "atomic"){
            coefsK2 = latticeToAtomicGauge(eigvecKStack.slice(ki_index).col(vi), kpoints.row(ki_index));
            coefsK2Q = latticeToAtomicGauge(eigvecKQStack.slice(ki_index).col(ci), kpoints.row(ki_index));
        }
        else{
            coefsK2 = eigvecKStack.slice(ki_index).col(vi);
            coefsK2Q = eigvecKQStack.slice(ki_index).col(ci);
        }

        std::complex<double> D, X;
        if (mode == "realspace"){
            arma::cx_mat motifFT = ftMotifStack.slice(ki_index);
            D = exactInteractionTermMFT(coefsKQ, coefsK2, coefsK2Q, coefsK, motifFT);
            X = 0;

        }
        else if (mode == "reciprocalspace"){
            arma::rowvec k2 = kpoints.row(ki_index);
            D = interactionTermFT(coefsK, coefsK2, coefsKQ, coefsK2Q, k, k2, k, k2, this->nReciprocalVectors);
            X = 0;
        }
        
        W(i) = - (D - X);                
    };

    double delta = 2.0/targetExciton.nk; // Adjust delta depending on number of k points
    double rho = targetExciton.pairDensityOfStates(energy, delta);
    cout << "DoS value: " << rho << endl;
    double hbar = 6.582119624E-16; // Units are eV*s

    transitionRate = (ncell*a)*2*PI*std::norm(arma::dot(W, initialState))*rho/hbar;

    return transitionRate;
}

/**
 * Method to print information about the exciton.
 * @return void 
 */
void Exciton::printInformation(){
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
    cout << "\n" << endl;

    cout << std::left << std::setw(30) << "Gauge used: " << gauge << endl;
    cout << std::left << std::setw(30) << "Calculation mode: " << mode << endl;
    if(mode == "reciprocalspace"){
        cout << std::left << std::setw(30) << "nG: " << nReciprocalVectors << endl;
    }
    if(exchange){
        cout << std::left << std::setw(30) << "Exchange: " << (exchange ? "True" : "False") << endl;
    }
    if(arma::norm(Q) > 1E-7){
        cout << std::left << std::setw(30) << "Q: "; 
        for (auto qi : Q){
            cout << qi << "  ";
        }
        cout << endl;
    }
    cout << std::left << std::setw(30) << "Scissor cut: " << scissor_ << endl;
}

}