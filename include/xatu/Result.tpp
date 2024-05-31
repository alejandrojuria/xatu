#ifndef RESULT_TPP
#define RESULT_TPP

#include "xatu/Result.hpp"

namespace xatu {

template <typename T>
Result<T>::Result(Exciton<T>* exciton, arma::vec& eigval, arma::cx_mat& eigvec) : 
    exciton_(exciton), system_(exciton->system) {
        m_eigval = eigval;
        m_eigvec = eigvec;
};

/* -------------------- Observables -------------------- */

/**
 * Computes the kinetic energy of the specified exciton.
 * @details The kinetic energy is defined as the part of the exciton energy that is coming from
 * the bands. It is computed using the diagonal part of the BSE matrix without the interactions, HK.
 * @param stateindex Index of exciton.
 * @return Kinetic energy of the exciton.
 */
template <typename T>
double Result<T>::kineticEnergy(int stateindex){
    arma::cx_vec coefs = eigvec.col(stateindex);
    arma::vec HK = arma::zeros(exciton->dimBSE);

    for (uint32_t n = 0; n < exciton->dimBSE; n++){
        int v = exciton_->bandToIndex[exciton->basisStates(n, 0)];
        int c = exciton_->bandToIndex[exciton->basisStates(n, 1)];
        uint32_t k_index = exciton->basisStates(n, 2);
        uint32_t kQ_index = k_index;

        HK(n) = exciton->eigvalKQStack.col(kQ_index)(c) - exciton->eigvalKStack.col(k_index)(v);
    }

    std::complex<double> energy = arma::cdot(coefs, HK % coefs);
    return energy.real();
}

/**
 * Computes the potential energy of the specified exciton.
 * @details The potential energy is the part of the exciton energy that comes from the
 * electrostatic interaction. It is calculated from V = E - K.
 * @param stateindex Index of exciton.
 * @return Potential energy of the exciton.
 */
template <typename T>
double Result<T>::potentialEnergy(int stateindex){
    std::complex<double> energy = eigval(stateindex) - kineticEnergy(stateindex);
    return energy.real();
}

/**
 * Computes the binding energy of the specified exciton. Requires specifying the gap
 * of the system.
 * @details The binding energy of an exciton is defined as the difference between the excitation energy 
 * of the exciton and the gap, i.e. EB = E - G
 * @param stateindex Index of the exciton.
 * @param gap Gap of the system
 * @return Binding energy of the exciton.
 */
template <typename T>
double Result<T>::bindingEnergy(int stateindex, double gap){
    double energy;
    if (gap == -1){
        gap = determineGap();
    }
    else if (gap < 0){
        std::cout << "Provided gap value must be positive" << std::endl;
    }

    energy = eigval(stateindex) - gap;
    return energy;
}

/**
 * Routine to compute the gap from the bands based on the position of the centre
 * of the exciton and the bands used. Beware: The gap is computed using only the 
 * bands that are used in the exciton formation.
 * @return Gap of the system.
 */
template <typename T>
double Result<T>::determineGap(){
    int stateindex = 0; // Ground state
    int kIndex = findExcitonPeak(stateindex);
    int valence = exciton->bandToIndex[exciton->valenceBands.max()];
    int conduction = exciton->bandToIndex[exciton->conductionBands.min()];

    double gap = exciton->eigvalKStack.col(kIndex)(conduction) - 
                 exciton->eigvalKStack.col(kIndex)(valence);
    return gap;
}


/**
 * Method to compute the Fourier transform of the exciton envelope function A(k).
 * @param stateindex Index of exciton eigenstate.
 * @param electron_position Position of the electron where we evaluate the Fourier transform.
 * @param hole_position Position of the hole.
 * @return Fourier transform of the amplitude, evaluated at Re - Rh. 
 */
template <typename T>
double Result<T>::ftExcitonEnvelope(int stateindex, const arma::rowvec& electron_position, 
                                       const arma::rowvec& hole_position){

    arma::cx_vec coefs = eigvec.col(stateindex);

    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = system->basisdim;
    int nk = system->nk;
    arma::cx_double ft;
    std::complex<double> imag(0,1);
    // Matrix dimension is: #orbitals*#motif*#orbitals

    for (int i = 0; i < nk; i++){

        kBlockEnd += (exciton->dimBSE/nk);
        
        arma::cx_vec A = coefs(arma::span(kBlock, kBlockEnd - 1));
        //arma::cx_vec v = exciton->eigvecKStack.slice(i).col();
        //arma::cx_vec c = exciton->eigvecKStack.slice(i).col();
        double summed_coefs = pow(arma::norm(A), 2);

        ft += summed_coefs * std::exp(-imag*arma::dot(system->kpoints.col(i), electron_position - hole_position));

        kBlock = kBlockEnd;
    };
    double result = std::real(ft);

    return result;
}

/** 
 * Routine to compute the expected Sz spin value of the electron
 * and hole that form a given exciton.
 * @param stateIndex Index of the exciton state.
 * @return Vector with the total spin of the exciton, the spin of the hole and that of the electron
 */
template <typename T>
arma::cx_vec Result<T>::spinX(int stateIndex){
    arma::cx_vec coefs = eigvec.col(stateIndex);
    arma::cx_vec spin = spinX(coefs);

    return spin;
}

/* -------------------- Output -------------------- */

/**
 * Method to write to file the exciton reciprocal amplitude (squared k-wavefunction).
 * @details The squared reciprocal wavefunction at each kpoint is given by the sum of the square of the
 * electron-hole ampltiudes for different pairs of bands at one same k.
 * @param statecoefs State corresponding to array of electron-hole pair coefficients, not
 * necessarily an exciton eigenstate.
 * @param textfile Pointer to file to write the reciprocal amplitude.
 */
template <typename T>
void Result<T>::writeReciprocalAmplitude(const arma::cx_vec& statecoefs, FILE* textfile){
    fprintf(textfile, "kx\tky\tkz\tProb.\n");
    int nbandsCombinations = exciton->conductionBands.n_elem * exciton->valenceBands.n_elem;
    for (int32_t i = 0; i < system->kpoints.n_cols; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(system->kpoints.col(1) - system->kpoints.col(0)); // L2 norm instead of l2
        fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    system->kpoints.col(i)(0), system->kpoints.col(i)(1), system->kpoints.col(i)(2), coef);
    };
    fprintf(textfile, "#\n");
}

/**
 * Method to write the reciprocal amplitude of one exciton eigenstate.
 * @details Overload of the method; uses the previous one (more general), to write the
 * reciprocal amplitude for an eigenstate.
 * @param stateindex Index of exciton.
 * @param textfile Pointer to file to write.
 */
template <typename T>
void Result<T>::writeReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeReciprocalAmplitude(statecoefs, textfile);
}

/**
 * Method to write the reciprocal amplitude of a given state on a extended Brillouin zone.
 * @details This method determines the minimum box that bounds the reciprocal unit cell, and
 * then writes the reciprocal amplitude on each point of the box using the periodicity of the BZ.
 * @param statecoefs Coefficients of state.
 * @param textfile Pointer to file. 
 */
template <typename T>
void Result<T>::writeExtendedReciprocalAmplitude(const arma::cx_vec& statecoefs, FILE* textfile){
    int nbandsCombinations = exciton->conductionBands.n_elem * exciton->valenceBands.n_elem;
    double boxLimit = boundingBoxBZ();

    for (int32_t i = 0; i < system->kpoints.n_cols; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(system->kpoints.col(1) - system->kpoints.col(0)); // L2 norm instead of l2
        arma::mat cells = system->generateCombinations(3, true);
        for(unsigned int n = 0; n < cells.n_cols; n++){
            arma::colvec cell = arma::colvec(3);
            for(int j = 0; j < system->ndim; j++){
                cell += cells.col(n)(j)*system->Gbasis.col(j);
            }
            arma::colvec displaced_k = system->kpoints.col(i) + cell;
            if(abs(displaced_k(0)) < boxLimit && abs(displaced_k(1)) < boxLimit){
                fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    displaced_k(0), displaced_k(1), displaced_k(2), coef);
            }
        }
        
    }
    fprintf(textfile, "#\n");
}

/**
 * Writes the extended reciprocal amplitude of a state, given its index, to a text file, 
 * using k-points within a bounding box of the Brillouin zone.
 * @param stateindex Index of the state in the eigenvectors matrix
 * @param textfile Pointer to a file where the extended reciprocal amplitude will be written
 * @return void
 */
template <typename T>
void Result<T>::writeExtendedReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeExtendedReciprocalAmplitude(statecoefs, textfile);
}

/**
 * Method to write the phase and module of each exciton coefficient. 
 * Intended to be used with excitons formed only with one electron-hole pair for each k. 
 * @details Throws an error if used with excitons formed by multiple bands.
 * @param statecoefs Coefficients of state (not necessarily an exciton eigenstate).
 * @param textfile Pointer to file.
 */
template <typename T>
void Result<T>::writePhase(const arma::cx_vec& statecoefs, FILE* textfile){
    if(exciton->bandList.n_elem != 2){
        throw std::logic_error("writePhase requires only one valence and conduction bands");
    }
    fprintf(textfile, "kx\tky\tkz\tMod.\tArg.\n");
    int nbandsCombinations = exciton->conductionBands.n_elem * exciton->valenceBands.n_elem;
    double module, phase;
    for (int i = 0; i < system->kpoints.n_cols; i++){
        module = abs(statecoefs(i));
        phase = arg(statecoefs(i));
        fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    system->kpoints.row(i)(0), system->kpoints.row(i)(1), system->kpoints.row(i)(2), 
                    module, phase);
    }
    fprintf(textfile, "#\n");
}

/**
 * Method to write to a file the phase of an exciton eigenstate.
 * @details Overload of the general method; uses it to write the phase of the eigenstate
 * and as such will throw an error if used with excitons formed by more than one pair of bands.
 * @param stateindex Index of exciton.
 * @param textfile Pointer of file.
 */
template <typename T>
void Result<T>::writePhase(int stateindex, FILE* textfile){
    arma::cx_vec coefs = eigvec.col(stateindex);
    writePhase(coefs, textfile);
}

/**
 * Writes the extended phase of a given state to a text file, using k-points within a bounding box of the Brillouin zone.
 * Note that this method does not work with excitons formed by more than one pair of bands.
 * @param statecoefs Vector representing the coefficients of a given state (not necessarily an exciton eigenstate).
 * @param textfile Pointer to a file where the phases will be written
 * @throws std::logic_error if the number of valence and conduction bands is different from one (i.e. one pair of bands)
 * @return void
 */
template <typename T>
void Result<T>::writeExtendedPhase(const arma::cx_vec& statecoefs, FILE* textfile){
    if(exciton->bandList.n_elem != 2){
        throw std::logic_error("writeExtendedPhase requires only one valence and conduction bands");
    }
    fprintf(textfile, "kx\tky\tkz\tMod.\tArg.\n");
    int nbandsCombinations = exciton->conductionBands.n_elem * exciton->valenceBands.n_elem;
    double boxLimit = boundingBoxBZ();
    double module, phase;
    for (int32_t i = 0; i < system->kpoints.n_cols; i++){
        module = abs(statecoefs(i));
        module /= arma::norm(system->kpoints.col(1) - system->kpoints.col(0)); // L2 norm instead of l2

        phase = arg(statecoefs(i));
        
        arma::mat cells = system->generateCombinations(3, system->ndim, true);
        for(unsigned int n = 0; n < cells.n_cols; n++){
            arma::colvec cell = cells.col(n)(0)*system->reciprocalLattice.col(0) + 
                                cells.col(n)(1)*system->reciprocalLattice.col(1);
            arma::colvec displaced_k = system->kpoints.col(i) + cell;
            if(abs(displaced_k(0)) < boxLimit && abs(displaced_k(1)) < boxLimit){
                fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    displaced_k(0), displaced_k(1), displaced_k(2), module, phase);
            }
        }
        
    }
    fprintf(textfile, "#\n");
}

/**
 * Writes the phase of an exciton eigenstate over an extended BZ. Works only
 * with excitons formed with one pair of bands.
 * @details Overload of the method based on the most general one for any state.
 * @param stateindex Index of the state.
 * @param textfile Pointer to file.
 * @throws std::logic_error if the exciton have more than one pair of bands.
 * @return void
 */
template <typename T>
void Result<T>::writeExtendedPhase(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeExtendedPhase(statecoefs, textfile);
}

/**
 * Writes the probability density of finding the electron at a given position, having
 * fixed the position of the hole.
 * @param statecoefs State whose probability density we want to determine.
 * @param holeIndex Index of atom of the motif where we fix the hole.
 * @param holeCell Unit cell where we fix the hole.
 * @param textfile File to write the amplitudes.
 * @param ncells Number of unit cells where we compute the amplitudes.
 * @return void
 */
template <typename T>
void Result<T>::writeRealspaceAmplitude(const arma::cx_vec& statecoefs, int holeIndex,
                                     const arma::colvec& holeCell, FILE* textfile, int ncells){

    arma::colvec holePosition = system->motif.col(holeIndex).subvec(0, 2) + holeCell;
    fprintf(textfile, "%11.8lf\t%11.8lf\t%14.11lf\n", holePosition(0), holePosition(1), 0.0);

    double radius = arma::norm(system->Rbasis.col(0)) * ncells;
    arma::mat cellCombinations = system->truncateSupercell(exciton->ncell, radius);
    arma::vec coefs = arma::zeros(cellCombinations.n_cols*system->motif.n_cols);
    uint64_t it = 0;

    // Compute probabilities
    for(uint32_t cellIndex = 0; cellIndex < cellCombinations.n_cols; cellIndex++){
        arma::colvec cell = cellCombinations.col(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < system->motif.n_cols; atomIndex++){
            //coefs(it) = atomCoefficientSquared(atomIndex, cell, holeCell, RScoefs);
            coefs(it) = realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
    }

    // Write probabilities to file
    it = 0;
    for(uint32_t cellIndex = 0; cellIndex < cellCombinations.n_cols; cellIndex++){
        arma::colvec cell = cellCombinations.col(cellIndex);
        for(unsigned int atomIndex = 0; atomIndex < system->motif.n_cols; atomIndex++){
            arma::colvec position = system->motif.col(atomIndex).subvec(0, 2) + cell;
            fprintf(textfile, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            position(0), position(1), coefs(it));
            it++;
        }
    }
    fprintf(textfile, "#\n");                              
}

/**
 * Method to write the probability density of finding the electron of an exciton eigenstate,
 * having fixed the position of the hole.
 * @details Overload of the method for general states.
 * @param stateindex Index of exciton.
 * @param holeIndex Index of atom where we put the hole.
 * @param holeCell Unit cell where the hole is fixed.
 * @param textfile File to write the probability density.
 * @param ncells Number of unit cells where we compute the probability density.
 * @return void 
 */
template <typename T>
void Result<T>::writeRealspaceAmplitude(int stateindex, int holeIndex, 
                                     const arma::colvec& holeCell, FILE* textfile, int ncells){

    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeRealspaceAmplitude(statecoefs, holeIndex, holeCell, textfile, ncells);
}

/** 
 * Method to write the eigenvalues in ascending order into a file. 
 * @param textfile Pointer to file
 * @param n Number of eigenvalues to write. If not specified, all eigenvalues are written.
 * @return void 
 */
template <typename T>
void Result<T>::writeEigenvalues(FILE* textfile, int64_t n){

    if(n > exciton->dimBSE || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }

    fprintf(textfile, "%d\t", exciton->dimBSE);
    uint32_t maxEigval = (n == 0) ? exciton->dimBSE : n;  
    for(unsigned int i = 0; i < maxEigval; i++){
        fprintf(textfile, "%11.7lf\t", eigval(i));
    }
    fprintf(textfile, "\n");
}

/**
 * Method to write eigenstates in ascending order into a file. 
 * @param textfile Pointer to file.
 * @param n Optional argument to specify number of states to write to a file. 
 * @return void
 */
template <typename T>
void Result<T>::writeStates(FILE* textfile, int64_t n){
    if(n > exciton->dimBSE || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }
    // First write basis
    fprintf(textfile, "%d\n", exciton->dimBSE);
    for(unsigned int i = 0; i < exciton->dimBSE; i++){
        arma::irowvec state = exciton->basisStates.row(i);
        arma::colvec kpoint = system->kpoints.col(state(2));
        int v = state(0);
        int c = state(1);
        fprintf(textfile, "%11.7lf\t%11.7lf\t%11.7lf\t%d\t%d\n", 
                kpoint(0), kpoint(1), kpoint(2), v, c);
    }

    uint32_t nstates = (n == 0) ? exciton->dimBSE : n;  
    for(uint32_t i = 0; i < nstates; i++){
        for(uint32_t j = 0; j < exciton->dimBSE; j++){
            fprintf(textfile, "%11.7lf\t%11.7lf\t", 
                    std::real(eigvec.col(i)(j)), std::imag(eigvec.col(i)(j)));
        }
        fprintf(textfile, "\n");
    }
}

/**
 * Writes the total, electron and hole spin of the first n excitons.
 * @param n Number of excitons to compute and write spin.
 * @param textfile Textfile where the spins are written.
 */
template <typename T>
void Result<T>::writeSpin(int64_t n, FILE* textfile){

    if(n > exciton->dimBSE || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }

    int64_t maxState = (n == 0) ? exciton->dimBSE : n;  
    fprintf(textfile, "n\tSt\tSe\tSh\n");
    for(int64_t i = 0; i < maxState; i++){
        auto spin = spinX(i);
        fprintf(textfile, "%ld\t%11.7lf\t%11.7lf\t%11.7lf\n", i, real(spin(0)), real(spin(1)), real(spin(2)));
    }
}

/**
 * Routine to find the k index where the exciton has the maximum amplitude, which
 * usually corresponds with the band gap location.
 * @details This routine might fail in the multiband case since it does not use the
 * reciprocal wavefunction to determine the kpoint.
 * @param stateindex Index of exciton.
 * @return Index of kpoint where the exciton peaks.
 */ 
template <typename T>
int Result<T>::findExcitonPeak(int32_t stateindex){
    int32_t index = eigvec.col(stateindex).index_max();
    int32_t bandCombinations = exciton->valenceBands.n_elem*exciton->conductionBands.n_elem;
    index = (int32_t)index/bandCombinations;
    return index;
}

/**
 * Routine to determine the size of a box that contains the reciprocal unit cell as given
 * by the BZ mesh. Intended to use with full BZ meshes. 
 * @return Half the side of the box.
 */ 
template <typename T>
double Result<T>::boundingBoxBZ(){
    double max_x = arma::max(system->kpoints.row(0));
    double max_y = arma::max(system->kpoints.row(1));

    double value = (max_x > max_y) ? max_x : max_y;
    return value;
}

}

#endif