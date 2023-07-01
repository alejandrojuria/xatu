#include "xatu/Result.hpp"
#include <complex>

namespace xatu {

/**
 * Computes the kinetic energy of the specified exciton.
 * @details The kinetic energy is defined as the part of the exciton energy that is coming from
 * the bands. It is computed using the diagonal part of the BSE matrix without the interactions, HK.
 * @param stateindex Index of exciton.
 * @return Kinetic energy of the exciton.
 */
double Result::kineticEnergy(int stateindex){
    arma::cx_vec coefs = eigvec.col(stateindex);
    std::complex<double> energy = arma::cdot(coefs, exciton.HK*coefs);
    return energy.real();
}

/**
 * Computes the potential energy of the specified exciton.
 * @details The potential energy is the part of the exciton energy that comes from the
 * electrostatic interaction. It is calculated from V = E - K.
 * @param stateindex Index of exciton.
 * @return Potential energy of the exciton.
 */
double Result::potentialEnergy(int stateindex){
    arma::cx_vec coefs = eigvec.col(stateindex);
    arma::cx_mat HV = exciton.HBS - exciton.HK;
    std::complex<double> energy = arma::cdot(coefs, HV*coefs);
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
double Result::bindingEnergy(int stateindex, double gap){
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
double Result::determineGap(){
    int stateindex = 0; // Ground state
    int kIndex = findExcitonPeak(stateindex);
    int valence = exciton.bandToIndex[exciton.valenceBands.max()];
    int conduction = exciton.bandToIndex[exciton.conductionBands.min()];

    double gap = exciton.eigvalKStack.col(kIndex)(conduction) - 
                 exciton.eigvalKStack.col(kIndex)(valence);
    return gap;
}

/**
 * Routine to find the k index where the exciton has the maximum amplitude, which
 * usually corresponds with the band gap location.
 * @details This routine might fail in the multiband case since it does not use the
 * reciprocal wavefunction to determine the kpoint.
 * @param stateindex Index of exciton.
 * @return Index of kpoint where the exciton peaks.
 */ 
int Result::findExcitonPeak(int stateindex){
    int index = eigvec.col(stateindex).index_max();
    int bandCombinations = exciton.valenceBands.n_elem*exciton.conductionBands.n_elem;
    index = (int)index/bandCombinations;
    return index;
}

/** 
 * Routine to compute the expected Sz spin value of the electron
 * and hole that form a given exciton.
 * @param stateindex Index of the exciton.
 * @return Vector with the spin of the hole, of the electron and the total spin of the exciton.
 */
arma::cx_vec Result::spinX(int stateindex){

    arma::cx_vec coefs = eigvec.col(stateindex);
    
    // Initialize Sz for both electron and hole to zero
    arma::cx_double electronSpin = 0;
    arma::cx_double holeSpin = 0;
    double totalSpin = 0;
    int dimX = exciton.basisStates.n_rows;

    arma::cx_vec spinEigvalues = {1./2, -1./2};
    arma::cx_vec spinVector = arma::zeros<arma::cx_vec>(exciton.basisdim);
    int vecIterator = 0;
    for(int atomIndex = 0; atomIndex < exciton.natoms; atomIndex++){
        int species = exciton.motif.row(atomIndex)(3);
        int norbitals = exciton.orbitals(species);
        spinVector.subvec(vecIterator, vecIterator + norbitals - 1) = 
                          arma::kron(arma::ones(norbitals/2), spinEigvalues);
        vecIterator += exciton.orbitals(species);
    }
    
	arma::cx_vec eigvec, spinEigvec;

    // Initialize hole spin and electron spin operators
    int nvbands = exciton.valenceBands.n_elem;
    int ncbands = exciton.conductionBands.n_elem;
    int npairs = nvbands*ncbands;

    arma::cx_mat spinHole = arma::zeros<arma::cx_mat>(dimX, dimX);
    arma::cx_mat spinElectron = arma::zeros<arma::cx_mat>(dimX, dimX);

    arma::cx_mat vMatrix = arma::eye<arma::cx_mat>(nvbands, nvbands);
    arma::cx_mat cMatrix = arma::eye<arma::cx_mat>(ncbands, ncbands);

    // Initialize list of pairs of valence-conduction bands
    arma::mat bandPairs = arma::zeros(npairs, 2);
    int i = 0;
    for(double v : exciton.valenceBands){
        for(double c : exciton.conductionBands){
            bandPairs.row(i) = arma::rowvec{v, c};
            i++;
        }
    }

    for(unsigned int k = 0; k < exciton.kpoints.n_rows; k++){
        arma::cx_mat spinHoleReduced = arma::zeros<arma::cx_mat>(nvbands, nvbands);
        arma::cx_mat spinElectronReduced = arma::zeros<arma::cx_mat>(ncbands, ncbands);
        for(int i = 0; i < nvbands; i++){
            int vIndex = exciton.bandToIndex[exciton.valenceBands(i)];
            for(int j = 0; j < nvbands; j++){
                int vIndex2 = exciton.bandToIndex[exciton.valenceBands(j)];
                eigvec = exciton.eigvecKStack.slice(k).col(vIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = exciton.eigvecKStack.slice(k).col(vIndex2);
                spinHoleReduced(i,j) = arma::cdot(eigvec, spinEigvec);
            }
        }
        for(int i = 0; i < ncbands; i++){
            int cIndex = exciton.bandToIndex[exciton.conductionBands(i)];
            for(int j = 0; j < ncbands; j++){
                int cIndex2 = exciton.bandToIndex[exciton.conductionBands(j)];
                eigvec = exciton.eigvecKQStack.slice(k).col(cIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = exciton.eigvecKQStack.slice(k).col(cIndex2);
                spinElectronReduced(i,j) = arma::cdot(eigvec, spinEigvec);
            }
        }
            
        spinHole.submat(k*npairs, k*npairs, (k+1)*npairs - 1, (k+1)*npairs - 1) = arma::kron(cMatrix, spinHoleReduced);
        spinElectron.submat(k*npairs, k*npairs, (k+1)*npairs - 1, (k+1)*npairs - 1) = arma::kron(spinElectronReduced, vMatrix);
    }

    // Perform tensor products with the remaining quantum numbers
    holeSpin = -arma::cdot(coefs, spinHole*coefs);
    electronSpin = arma::cdot(coefs, spinElectron*coefs);
    totalSpin = real((holeSpin + electronSpin));
    
    arma::cx_vec results = {totalSpin, holeSpin, electronSpin};
    return results;
}

/**
 * Method to diagonalize the C3 rotation operator in a exciton degenerate subspace.
 * @details This method is intended to be used with systems with C3 rotational symmetry;
 * note that its current implementation is not correct (lacking additional phases).
 * @param states Vector storing the indices of the states of the degenerate subspace.
 * @return Eigenvectors of C3 in the degenerate exciton basis.
 */
arma::cx_mat Result::diagonalizeC3(const arma::vec& states){
    arma::mat C3 = exciton.C3ExcitonBasisRep();
    arma::cx_mat degenerateSubspaceC3 = arma::zeros<arma::cx_mat>(states.n_elem, states.n_elem);
    arma::cx_vec state = eigvec.col(states(0));
    arma::cout << "First C3: " << arma::cdot(state, C3*state) << arma::endl;
    arma::cout << "Second C3: " << arma::cdot(state, C3*C3*state) << arma::endl;
    arma::cout << "Third C3: " << arma::cdot(state, C3*C3*C3*state) << arma::endl;
    for(unsigned int i = 0; i < states.n_elem; i++){
        for(unsigned int j = 0; j < states.n_elem; j++){
            arma::cx_vec rowState = eigvec.col(states(i));
            arma::cx_vec columnState = eigvec.col(states(j));

            degenerateSubspaceC3(i, j) = arma::cdot(rowState, C3*columnState);
        }
    }

    std::cout << degenerateSubspaceC3 << std::endl;
    arma::cx_vec eigvalC3;
    arma::cx_mat eigvecC3;
    arma::eig_gen(eigvalC3, eigvecC3, degenerateSubspaceC3);
    std::cout << "C3 eigval:\n" << eigvalC3 << std::endl;
    return eigvecC3;
}

/** TODO: DELETE
 * Method to symmetrize exciton states of a degenerate subspace according to C3 
 * rotational symmetry. Beware: this method most likely returns incorrect results.
 * @param state
 */
arma::cx_mat Result::symmetrizeStates(const arma::cx_vec& state, const arma::cx_vec& degState){
    arma::mat C3 = exciton.C3ExcitonBasisRep();
    double alpha, phase;
    arma::vec values = arma::linspace(0, 1, 100);
    arma::vec phaseValues = arma::linspace(0, 2*PI, 100);
    std::complex<double> imag(0, 1);
    arma::cx_vec linearCombination, rotated;
    double target = -1;
    double scalar;
    for(double value : values){
        for(double theta : phaseValues){
            linearCombination = state*value*exp(imag*theta) + degState*sqrt(1 - value*value);
            rotated = C3*linearCombination;
            scalar = abs(arma::dot(linearCombination, rotated));
            if(target < scalar){
                target = scalar;
                alpha = value;
            }
        }
    }
    std::cout << "alpha value is: " << alpha << std::endl;
    arma::cx_vec symmetrizedState    = state*alpha*exp(imag*phase) + degState*sqrt(1 - alpha*alpha);
    arma::cx_vec degSymmetrizedState = state*sqrt(1 - alpha*alpha) - degState*alpha*exp(imag*phase);
    arma::cx_mat states(exciton.excitonbasisdim, 2);
    states.col(0) = symmetrizedState;
    states.col(1) = degSymmetrizedState;

    return states;
}


/**
 * Method to write to file the exciton reciprocal amplitude (squared k-wavefunction).
 * @details The squared reciprocal wavefunction at each kpoint is given by the sum of the square of the
 * electron-hole ampltiudes for different pairs of bands at one same k.
 * @param statecoefs State corresponding to array of electron-hole pair coefficients, not
 * necessarily an exciton eigenstate.
 * @param textfile Pointer to file to write the reciprocal amplitude.
 */
void Result::writeReciprocalAmplitude(const arma::cx_vec& statecoefs, FILE* textfile){
    fprintf(textfile, "kx\tky\tkz\tProb.\n");
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    for (int i = 0; i < exciton.kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.kpoints.row(1) - exciton.kpoints.row(0)); // L2 norm instead of l2
        fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    exciton.kpoints.row(i)(0), exciton.kpoints.row(i)(1), exciton.kpoints.row(i)(2), coef);
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
void Result::writeReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeReciprocalAmplitude(statecoefs, textfile);
};

/**
 * Method to write the phase and module of each exciton coefficient. 
 * Intended to be used with excitons formed only with one electron-hole pair for each k. 
 * @details Throws an error if used with excitons formed by multiple bands.
 * @param statecoefs Coefficients of state (not necessarily an exciton eigenstate).
 * @param textfile Pointer to file.
 */
void Result::writePhase(const arma::cx_vec& statecoefs, FILE* textfile){
    if(exciton.bandList.n_elem != 2){
        throw std::logic_error("writePhase requires only one valence and conduction bands");
    }
    fprintf(textfile, "kx\tky\tkz\tMod.\tArg.\n");
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    double module, phase;
    for (int i = 0; i < exciton.kpoints.n_rows; i++){
        module = abs(statecoefs(i));
        phase = arg(statecoefs(i));
        fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    exciton.kpoints.row(i)(0), exciton.kpoints.row(i)(1), exciton.kpoints.row(i)(2), 
                    module, phase);
    };
    fprintf(textfile, "#\n");
}

/**
 * Method to write to a file the phase of an exciton eigenstate.
 * @details Overload of the general method; uses it to write the phase of the eigenstate
 * and as such will throw an error if used with excitons formed by more than one pair of bands.
 * @param stateindex Index of exciton.
 * @param textfile Pointer of file.
 */
void Result::writePhase(int stateindex, FILE* textfile){
    arma::cx_vec coefs = eigvec.col(stateindex);
    writePhase(coefs, textfile);
}

/**
 * Method to write the reciprocal amplitude of a given state on a extended Brillouin zone.
 * @details This method determines the minimum box that bounds the reciprocal unit cell, and
 * then writes the reciprocal amplitude on each point of the box using the periodicity of the BZ.
 * @param statecoefs Coefficients of state.
 * @param textfile Pointer to file. 
 */
void Result::writeExtendedReciprocalAmplitude(const arma::cx_vec& statecoefs, FILE* textfile){
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    double boxLimit = boundingBoxBZ();

    for (int i = 0; i < exciton.kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.kpoints.row(1) - exciton.kpoints.row(0)); // L2 norm instead of l2
        arma::mat cells = exciton.generateCombinations(3, exciton.ndim, true);
        for(unsigned int n = 0; n < cells.n_rows; n++){
            arma::rowvec cell = arma::rowvec(3);
            for(int j = 0; j < exciton.ndim; j++){
                cell += cells.row(n)(j)*exciton.reciprocalLattice.row(j);
            }
            arma::rowvec displaced_k = exciton.kpoints.row(i) + cell;
            if(abs(displaced_k(0)) < boxLimit && abs(displaced_k(1)) < boxLimit){
                fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    displaced_k(0), displaced_k(1), displaced_k(2), coef);
            }
        }
        
    };
    fprintf(textfile, "#\n");
}

/**
 * Writes the extended reciprocal amplitude of a state, given its index, to a text file, 
 * using k-points within a bounding box of the Brillouin zone.
 * @param stateindex Index of the state in the eigenvectors matrix
 * @param textfile Pointer to a file where the extended reciprocal amplitude will be written
 * @return void
 */
void Result::writeExtendedReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeExtendedReciprocalAmplitude(statecoefs, textfile);
}

/**
 * Writes the extended phase of a given state to a text file, using k-points within a bounding box of the Brillouin zone.
 * Note that this method does not work with excitons formed by more than one pair of bands.
 * @param statecoefs Vector representing the coefficients of a given state (not necessarily an exciton eigenstate).
 * @param textfile Pointer to a file where the phases will be written
 * @throws std::logic_error if the number of valence and conduction bands is different from one (i.e. one pair of bands)
 * @return void
 */
void Result::writeExtendedPhase(const arma::cx_vec& statecoefs, FILE* textfile){
    if(exciton.bandList.n_elem != 2){
        throw std::logic_error("writeExtendedPhase requires only one valence and conduction bands");
    }
    fprintf(textfile, "kx\tky\tkz\tMod.\tArg.\n");
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    double boxLimit = boundingBoxBZ();
    double module, phase;
    for (int i = 0; i < exciton.kpoints.n_rows; i++){
        module = abs(statecoefs(i));
        module /= arma::norm(exciton.kpoints.row(1) - exciton.kpoints.row(0)); // L2 norm instead of l2

        phase = arg(statecoefs(i));
        
        arma::mat cells = exciton.generateCombinations(3, exciton.ndim, true);
        for(unsigned int n = 0; n < cells.n_rows; n++){
            arma::rowvec cell = cells.row(n)(0)*exciton.reciprocalLattice.row(0) + 
                                cells.row(n)(1)*exciton.reciprocalLattice.row(1);
            arma::rowvec displaced_k = exciton.kpoints.row(i) + cell;
            if(abs(displaced_k(0)) < boxLimit && abs(displaced_k(1)) < boxLimit){
                fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    displaced_k(0), displaced_k(1), displaced_k(2), module, phase);
            }
        }
        
    };
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
void Result::writeExtendedPhase(int stateindex, FILE* textfile){
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
void Result::writeRealspaceAmplitude(const arma::cx_vec& statecoefs, int holeIndex,
                                     const arma::rowvec& holeCell, FILE* textfile, int ncells){

    arma::rowvec holePosition = exciton.motif.row(holeIndex).subvec(0, 2) + holeCell;
    fprintf(textfile, "%11.8lf\t%11.8lf\t%14.11lf\n", holePosition(0), holePosition(1), 0.0);

    double radius = arma::norm(exciton.bravaisLattice.row(0)) * ncells;
    arma::mat cellCombinations = exciton.truncateSupercell(exciton.ncell, radius);
    arma::vec coefs = arma::zeros(cellCombinations.n_rows*exciton.motif.n_rows);
    int it = 0;

    // Compute probabilities
    for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.motif.n_rows; atomIndex++){
            //coefs(it) = atomCoefficientSquared(atomIndex, cell, holeCell, RScoefs);
            coefs(it) = realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
    }

    // Write probabilities to file
    it = 0;
    for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for(unsigned int atomIndex = 0; atomIndex < exciton.motif.n_rows; atomIndex++){
            arma::rowvec position = exciton.motif.row(atomIndex).subvec(0, 2) + cell;
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
void Result::writeRealspaceAmplitude(int stateindex, int holeIndex, 
                                     const arma::rowvec& holeCell, FILE* textfile, int ncells){

    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeRealspaceAmplitude(statecoefs, holeIndex, holeCell, textfile, ncells);
}

/** 
 * Method to write the eigenvalues in ascending order into a file. 
 * @param textfile Pointer to file
 * @param n Number of eigenvalues to write. If not specified, all eigenvalues are written.
 * @return void 
 */
void Result::writeEigenvalues(FILE* textfile, int n){

    if(n > exciton.excitonbasisdim || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }

    fprintf(textfile, "%d\t", exciton.excitonbasisdim);
    int maxEigval = (n == 0) ? exciton.excitonbasisdim : n;  
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
void Result::writeStates(FILE* textfile, int n){
    if(n > exciton.excitonbasisdim || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }
    // First write basis
    fprintf(textfile, "%d\n", exciton.excitonbasisdim);
    for(unsigned int i = 0; i < exciton.excitonbasisdim; i++){
        arma::irowvec state = exciton.basisStates.row(i);
        arma::rowvec kpoint = exciton.kpoints.row(state(2));
        int v = state(0);
        int c = state(1);
        fprintf(textfile, "%11.7lf\t%11.7lf\t%11.7lf\t%d\t%d\n", 
                kpoint(0), kpoint(1), kpoint(2), v, c);
    }

    int nstates = (n == 0) ? exciton.excitonbasisdim : n;  
    for(unsigned int i = 0; i < nstates; i++){
        for(unsigned int j = 0; j < exciton.excitonbasisdim; j++){
            fprintf(textfile, "%11.7lf\t%11.7lf\t", 
                    real(eigvec.col(i)(j)), imag(eigvec.col(i)(j)));
        }
        fprintf(textfile, "\n");
    }
}

/**
 * Method to compute and write the absorption spectra to a file.
 * @details This method computes both the single particle absorption, and the absorption
 * from the exciton spectrum. All the required parameters must be specified in a separate text file
 * named kubo_w.in
 * @return void 
 */
void Result::writeAbsorptionSpectrum(){

    int nR = exciton.unitCellList.n_rows;
    int norb = exciton.basisdim;
    int norb_ex = exciton.excitonbasisdim;
    int filling = exciton.filling;
    int nv = exciton.valenceBands.n_elem;
    int nc = exciton.conductionBands.n_elem;

    arma::mat Rvec = exciton.unitCellList;
    // Extend bravais lattice to 3x3 matrix
    arma::mat R = arma::zeros(3, 3);
    for (int i = 0; i < exciton.bravaisLattice.n_rows; i++){
        R.row(i) = exciton.bravaisLattice.row(i);
    }

    // arma::mat B = exciton.motif.cols(0, 2);
    arma::mat extendedMotif = arma::zeros(exciton.basisdim, 3);
    int it = 0;
    for(int i = 0; i < exciton.natoms; i++){
        arma::rowvec atom = exciton.motif.row(i).subvec(0, 2);
        int species = exciton.motif.row(i)(3);
        for(int j = 0; j < exciton.orbitals(species); j++){
            extendedMotif.row(it) = atom; 
            it++;
        }
    }
    arma::cx_cube hhop = exciton.hamiltonianMatrices;
    arma::cube shop(arma::size(hhop));
    if (exciton.overlapMatrices.empty()){
        for (int i = 0; i < hhop.n_slices; i++){
            shop.slice(i) = arma::eye(size(hhop.slice(i)));
        }
    }
    else{
        shop = arma::real(exciton.overlapMatrices);
    }
    int nk = exciton.nk;
    arma::vec rkx = exciton.kpoints.col(0);
    arma::vec rky = exciton.kpoints.col(1);
    arma::vec rkz = exciton.kpoints.col(2);

    arma::mat eigval_sp = exciton.eigvalKStack;
    arma::cx_cube eigvec_sp = exciton.eigvecKStack;

    skubo_w_(&nR, &norb, &norb_ex, &nv, &nc, &filling, 
             Rvec.memptr(), R.memptr(), extendedMotif.memptr(), hhop.memptr(), shop.memptr(), &nk, rkx.memptr(),
             rky.memptr(), rkz.memptr(), m_eigvec.memptr(), m_eigval.memptr(), eigval_sp.memptr(), eigvec_sp.memptr());
}

/**
 * Writes the total, electron and hole spin of the first n excitons.
 * @param n Number of excitons to compute and write spin.
 * @param textfile Textfile where the spins are written.
 */
void Result::writeSpin(int n, FILE* textfile){

    if(n > exciton.excitonbasisdim || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }

    int maxState = (n == 0) ? exciton.excitonbasisdim : n;  
    fprintf(textfile, "n\tSt\tSe\tSh\n");
    for(unsigned int i = 0; i < maxState; i++){
        auto spin = spinX(i);
        fprintf(textfile, "%d\t%11.7lf\t%11.7lf\t%11.7lf\n", i, real(spin(0)), real(spin(1)), real(spin(2)));
    }
}

/**
 * Method to compute the Fourier transform of the exciton envelope function A(k).
 * @param stateindex Index of exciton eigenstate.
 * @param electron_position Position of the electron where we evaluate the Fourier transform.
 * @param hole_position Position of the hole.
 * @return Fourier transform of the amplitude, evaluated at Re - Rh. 
 */
double Result::fourierTransformExciton(int stateindex, const arma::rowvec& electron_position, 
                                       const arma::rowvec& hole_position){

    arma::cx_vec coefs = eigvec.col(stateindex);

    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = exciton.basisdim;
    int nk = exciton.nk;
    arma::cx_double ft;
    std::complex<double> imag(0,1);
    // Matrix dimension is: #orbitals*#motif*#orbitals

    for (int i = 0; i < nk; i++){

        kBlockEnd += (exciton.excitonbasisdim/nk);
        
        arma::cx_vec A = coefs(arma::span(kBlock, kBlockEnd - 1));
        //arma::cx_vec v = exciton.eigvecKStack.slice(i).col();
        //arma::cx_vec c = exciton.eigvecKStack.slice(i).col();
        double summed_coefs = pow(arma::norm(A), 2);

        ft += summed_coefs * std::exp(-imag*arma::dot(exciton.kpoints.row(i), electron_position - hole_position));

        kBlock = kBlockEnd;
    };
    double result = std::real(ft);

    return result;
}

/**
 * Routine to determine the size of a box that contains the reciprocal unit cell as given
 * by the BZ mesh. Intended to use with full BZ meshes. 
 * @return Half the side of the box.
 */ 
double Result::boundingBoxBZ(){
    double max_x = arma::max(exciton.kpoints.col(0));
    double max_y = arma::max(exciton.kpoints.col(1));

    double value = (max_x > max_y) ? max_x : max_y;
    return value;
}

/**
 * Method to compute the real-space amplitude of an exciton state (not necessarily an eigenstate).
 * @details Used by writeRealSpaceAmplitude to write the probability density over several unit cells.
 * @param BSEcoefs State whose real-space amplitude we want to obtain.
 * @param electronIndex Index of the atom where we put the electron.
 * @param holeIndex Index of atom where we put the hole.
 * @param eCell Unit cell of the electron.
 * @param hCell Unit cell of the hole.
 * @return Real-space amplitude evaluated at those electron and hole positions.
 */
double Result::realSpaceWavefunction(const arma::cx_vec& BSEcoefs, int electronIndex, int holeIndex,
                             const arma::rowvec& eCell, const arma::rowvec& hCell){

    std::complex<double> imag(0, 1);
    double totalAmplitude = 0;
    arma::cx_vec eigvec = arma::cx_vec(BSEcoefs);
    int eOrbitals = exciton.orbitals(exciton.motif.row(electronIndex)(3));
    int hOrbitals = exciton.orbitals(exciton.motif.row(holeIndex)(3));

    // Compute index corresponding to electron and hole
    int eIndex = 0;
    int hIndex = 0;
    for(unsigned int i = 0; i < electronIndex; i++){
        eIndex += exciton.orbitals(exciton.motif.row(i)(3));
    }
    for(unsigned int i = 0; i < holeIndex; i++){
        hIndex += exciton.orbitals(exciton.motif.row(i)(3));
    }
    eigvec = addExponential(eigvec, eCell - hCell);

    for(int alpha = 0; alpha < eOrbitals; alpha++){
        for(int beta = 0; beta < hOrbitals; beta++){

        arma::cx_cube c = exciton.eigvecKQStack.tube(eIndex + alpha, exciton.valenceBands.n_elem, 
                                eIndex + alpha, exciton.valenceBands.n_elem + exciton.conductionBands.n_elem - 1);
        arma::cx_rowvec cFlat = arma::reshape(c, 1, c.n_elem, 1);
        arma::cx_rowvec cExtended = arma::kron(cFlat, arma::ones<arma::cx_rowvec>(exciton.valenceBands.n_elem));

        arma::cx_cube v = exciton.eigvecKStack.tube(hIndex + beta, 0, 
                                hIndex + beta, exciton.valenceBands.n_elem - 1);

        arma::cx_rowvec vFlat = arma::reshape(v, 1, v.n_elem, 1);
        arma::cx_rowvec vExtended = arma::zeros<arma::cx_rowvec>(vFlat.n_elem*exciton.conductionBands.n_elem);
        int blockSize = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
        for(unsigned int i = 0; i < exciton.nk; i++){
            vExtended.subvec(i*blockSize, (i + 1)*blockSize - 1) = arma::kron(arma::ones<arma::cx_rowvec>(exciton.conductionBands.n_elem), 
                                                        vFlat.subvec(i*exciton.valenceBands.n_elem, (i + 1)*exciton.valenceBands.n_elem - 1));
        }
        
        arma::cx_rowvec coefs = cExtended % arma::conj(vExtended);

        totalAmplitude += std::norm(arma::dot(coefs, eigvec));    
        }
    }

    return totalAmplitude;
};

/**
 * Method to add exponentials to some vector of coefficients.
 * @details Used in realSpaceWavefunction to compute the real-space exciton amplitudes.
 * Basically multiplies each coefficient by an exponential with phase ikR.
 * @param coefs Vector of electron-hole pair coefficients.
 * @param cell Unit cell used in the exponential.
 * @return Coefficients with the added exponential.
 */
arma::cx_vec Result::addExponential(arma::cx_vec& coefs, const arma::rowvec& cell){

    arma::vec product = exciton.kpoints * cell.t();
    std::complex<double> imag(0, 1);
    arma::cx_vec exponentials = arma::exp(imag*product);
    int nBandCombinations = exciton.valenceBands.n_elem*exciton.conductionBands.n_elem;
    exponentials = arma::kron(exponentials, arma::ones<arma::cx_vec>(nBandCombinations));

    coefs = coefs % exponentials;

    return coefs;
}


/** 
 * Routine to compute the density matrix for the excitons.
 * @details This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
 * Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object. Note: deprecated, not used.
 * @param exciton Exciton object.
 * @param BSEcoefs Exciton state.
 * @param eIndex Index of atom of electron.
 * @param hIndex Index of atom of hole.
 * @return Density matrix matrix element.
 */
std::complex<double> Result::densityMatrix(Exciton& exciton, const arma::cx_vec& BSEcoefs, 
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

/** 
 * Density matrix with k (?)
 * @details This routine requires having all eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
 * Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object. Note: Deprecated, not used.
 * @param kIndex Index of kpoint used.
 * @param exciton Exciton object.
 * @param BSEcoefs State to be used.
 * @param eIndex Index of atom of the electron.
 * @param hIndex Index of hole.
 * @return Matrix elements of density matrix evaluated at k.
 */
std::complex<double> Result::densityMatrixK(int kIndex, Exciton& exciton, const arma::cx_vec& BSEcoefs, 
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

}