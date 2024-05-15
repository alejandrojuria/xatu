#pragma once
#include <math.h>
#include <iomanip>
#include <armadillo>
#include "xatu/Exciton.hpp"

namespace xatu {

/* ------------------------------ Setters ------------------------------ */
/**
 * Sets the number of unit cells along each axis.
 * @param ncell Number of unit cells per axis. 
 * @return void
 */
template <typename T>
void Exciton<T>::setUnitCells(int ncell){
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
template <typename T>
void Exciton<T>::setBands(const arma::ivec& bands){
    bands_ = bands;
    std::vector<arma::s64> valence, conduction;
    for(int i = 0; i < bands.n_elem; i++){
        if (bands(i) <= 0){
            valence.push_back(bands(i) + system->highestValenceBand);
        }
        else{
            conduction.push_back(bands(i) + system->highestValenceBand);
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
template <typename T>
void Exciton<T>::setBands(int nbands, int nrmbands){
    int highestValenceBand = system->highestValenceBand;
    if(nbands > 0 && nrmbands > 0){
        this->valenceBands_ = arma::regspace<arma::ivec>(highestValenceBand - nbands + 1, highestValenceBand - nrmbands);
        this->conductionBands_ = arma::regspace<arma::ivec>(highestValenceBand + 1 + nrmbands, highestValenceBand + nbands);
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
template <typename T>
void Exciton<T>::setQ(const arma::rowvec& Q){
    if(Q.n_elem == 3){
        this->Q_ = Q;
    }
    else{
        std::cout << "Q vector must be 3d" << std::endl;
    }
}

/**
 * Sets the cutoff over unit cells used in the calculation of the lattice Fourier transform
 * for the interactions.
 * @param cutoff Number of unit cells to consider.
 * @return void 
 */
template <typename T>
void Exciton<T>::setCutoff(double cutoff){
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
 * Sets the value of the scissor cut of change the gap of the system->
 * @param shift Value of scissor cut (in eV). Can be positive or negative.
 * @return void 
 */
template <typename T>
void Exciton<T>::setScissor(double shift){
    this->scissor_ = shift;
}

/**
 * To toggle on or off the exchange term in the interaction matrix elements.
 * @param exchange Either true of false
 * @return void
*/
template <typename T>
void Exciton<T>::setExchange(bool exchange){
    this->exchange = exchange;
}

/*------------------------------------ Electron-hole pair basis ------------------------------------*/

/**
 * Initialise basis to be used in the construction of the BSE matrix.
 * @param conductionBands Conduction bands that will populate the electrons of the exciton.
 * @param valenceBands Valence bands to be populated by holes of the exciton.
 * @return Matrix where each row denotes an electron-hole pair, '{v, c, k}'.
 */
template <typename T>
arma::imat Exciton<T>::createBasis(const arma::ivec& conductionBands, 
                                const arma::ivec& valenceBands){

    arma::imat states = arma::zeros<arma::imat>(dimBSE, 3);
    int it = 0;
    for (uint32_t i = 0; i < system->nk; i++){
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
template <typename T>
void Exciton<T>::initializeBasis(){
    this->basisStates_ = createBasis(conductionBands, valenceBands);
};

/**
 * Wrapper for the Brillouin zone mesh method in the System class,
 * to avoid using the pointer to the System object.
 * @param ncell Number of cells in the mesh.
*/
template <typename T>
void Exciton<T>::brillouinZoneMesh(int ncell){
    system->brillouinZoneMesh(ncell);
}

/**
 * Criterium to fix the phase of the single-particle eigenstates after diagonalization.
 * @details The prescription we take here is to impose that the sum of all the coefficients is real.
 * @return Fixed coefficients. 
 */
template <typename T>
arma::cx_mat Exciton<T>::fixGlobalPhase(arma::cx_mat& coefs){

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
template <typename T>
void Exciton<T>::generateBandDictionary(){

    std::map<int, int> bandToIndex;
    for(int i = 0; i < bandList.n_elem; i++){
        bandToIndex[bandList(i)] = i;
    };

    this->bandToIndex = bandToIndex;
};

/**
 * Method to print information about the exciton.
 * @return void 
 */
template <typename T>
void Exciton<T>::printInformation(){
    std::cout << std::left << std::setw(30) << "Number of cells: " << ncell << std::endl;
    std::cout << std::left << std::setw(30) << "Valence bands:";
    for (int i = 0; i < valenceBands.n_elem; i++){
        std::cout << valenceBands(i) << "\t";
    }
    std::cout << std::endl;

    std::cout << std::left << std::setw(30) << "Conduction bands: ";
    for (int i = 0; i < conductionBands.n_elem; i++){
        std::cout << conductionBands(i) << "\t";
    }
    std::cout << "\n" << std::endl;

    if(exchange){
        std::cout << std::left << std::setw(30) << "Exchange: " << (exchange ? "True" : "False") << std::endl;
    }
    if(arma::norm(Q) > 1E-7){
        std::cout << std::left << std::setw(30) << "Q: "; 
        for (auto qi : Q){
            std::cout << qi << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::left << std::setw(30) << "Scissor cut: " << scissor_ << std::endl;
}

}