#include "Result.hpp"
#include <complex>

double Result::kineticEnergy(int stateindex){
    arma::cx_vec coefs = eigvec.col(stateindex);
    std::complex<double> energy = arma::cdot(coefs, exciton.HK*coefs);
    return energy.real();
}

double Result::potentialEnergy(int stateindex){
    arma::cx_vec coefs = eigvec.col(stateindex);
    arma::cx_mat HV = exciton.HBS - exciton.HK;
    std::complex<double> energy = arma::cdot(coefs, HV*coefs);
    return energy.real();
}

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

// Routine to compute the gap from the bands based on the position of the centre
// of the exciton and the bands used. Beware: The gap is computed using only the 
// bands that are used in the exciton formation.
double Result::determineGap(){
    int stateindex = 0; // Ground state
    int kIndex = findExcitonPeak(stateindex);
    int valence = exciton.valenceBands.max();
    int conduction = exciton.conductionBands.max();
    double gap = exciton.eigvalKStack.col(kIndex)(conduction) - 
                 exciton.eigvalKStack.col(kIndex)(valence);
    return gap;
}

// Routine to find the k index where the exciton has the maximum amplitude, which
// usually corresponds with the band gap location
int Result::findExcitonPeak(int stateindex){
    int index = eigvec.col(stateindex).index_max();
    return index;
}

/* Routine to compute the expected Sz spin value of the electron
and hole that form a given exciton. */
arma::cx_vec Result::spinX(int stateindex){

    arma::cx_vec coefs = eigvec.col(stateindex);
    
    // Initialize Sz for both electron and hole to zero
    arma::cx_double electronSpin = 0;
    arma::cx_double holeSpin = 0;
    double totalSpin = 0;
    int dimX = exciton.basisStates.n_rows;

	arma::cx_vec spinEigvalues = {1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2, -1./2};
	arma::cx_vec spinVector = arma::kron(arma::ones(exciton.basisdim/exciton.norbitals, 1), spinEigvalues);
	arma::cx_vec eigvec, spinEigvec;

    // Initialize hole spin and electron spin operators
    int nbands = exciton.bandList.n_elem;
    int nbandsSq = nbands*nbands;

    arma::cx_mat spinHole = arma::zeros<arma::cx_mat>(dimX, dimX);
    arma::cx_mat spinElectron = arma::zeros<arma::cx_mat>(dimX, dimX);

    arma::cx_mat spinHoleReduced = arma::zeros<arma::cx_mat>(nbands, nbands);
    arma::cx_mat spinElectronReduced = arma::zeros<arma::cx_mat>(nbands, nbands);

    arma::cx_mat vMatrix = arma::eye<arma::cx_mat>(nbands, nbands);
    arma::cx_mat cMatrix = arma::eye<arma::cx_mat>(nbands, nbands);

    for(unsigned int k = 0; k < exciton.kpoints.n_rows; k++){

        for(int i = 0; i < nbands; i++){
            int vIndex = exciton.bandToIndex[exciton.valenceBands(i)];
            int cIndex = exciton.bandToIndex[exciton.conductionBands(i)];
            for(int j = 0; j < nbands; j++){
                int vIndex2 = exciton.bandToIndex[exciton.valenceBands(j)];
                int cIndex2 = exciton.bandToIndex[exciton.conductionBands(j)];
                eigvec = exciton.eigvecKStack.slice(k).col(vIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = exciton.eigvecKStack.slice(k).col(vIndex2);
                spinHoleReduced(i,j) = arma::cdot(eigvec, spinEigvec);

                eigvec = exciton.eigvecKQStack.slice(k).col(cIndex);
                spinEigvec = eigvec % spinVector;
                eigvec = exciton.eigvecKQStack.slice(k).col(cIndex2);
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
    
    arma::cx_vec results = {holeSpin, electronSpin, totalSpin};
    return results;
}

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

void Result::writeReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeReciprocalAmplitude(statecoefs, textfile);
};

/* Method to write the phase and module of each exciton coefficient. Note that this
routine is only intended to be used with excitons formed only with one electron-hole combination
for each k. */
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

void Result::writePhase(int stateindex, FILE* textfile){
    arma::cx_vec coefs = eigvec.col(stateindex);
    writePhase(coefs, textfile);
}

void Result::writeExtendedReciprocalAmplitude(const arma::cx_vec& statecoefs, FILE* textfile){
    fprintf(textfile, "kx\tky\tkz\tProb.\n");
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
            arma::rowvec cell = cells.row(n)(0)*exciton.reciprocalLattice.row(0) + 
                                cells.row(n)(1)*exciton.reciprocalLattice.row(1);
            arma::rowvec displaced_k = exciton.kpoints.row(i) + cell;
            if(abs(displaced_k(0)) < boxLimit && abs(displaced_k(1)) < boxLimit){
                fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    displaced_k(0), displaced_k(1), displaced_k(2), coef);
            }
        }
        
    };
    fprintf(textfile, "#\n");
}

void Result::writeExtendedReciprocalAmplitude(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeExtendedReciprocalAmplitude(statecoefs, textfile);
}

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

void Result::writeExtendedPhase(int stateindex, FILE* textfile){
    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeExtendedPhase(statecoefs, textfile);
}

void Result::writeRealspaceAmplitude(const arma::cx_vec& statecoefs, int holeIndex,
                                     const arma::rowvec& holeCell, FILE* textfile){

    std::cout << __LINE__ << std::endl;
    arma::rowvec holePosition = exciton.motif.row(holeIndex) + holeCell;
    std::cout << __LINE__ << std::endl;
    arma::cx_mat RScoefs = RScoefficientCalc(statecoefs, holeIndex);
    std::cout << __LINE__ << std::endl;
    fprintf(textfile, "%11.8lf\t%11.8lf\t%14.11lf\n", holePosition(0), holePosition(1), 0.0);

    double radius = arma::norm(exciton.bravaisLattice.row(0)) * exciton.cutoff;
    std::cout << __LINE__ << std::endl;
    arma::mat cellCombinations = exciton.truncateSupercell(exciton.ncell, radius);
    std::cout << __LINE__ << std::endl;
    arma::vec coefs = arma::zeros(cellCombinations.n_rows*exciton.motif.n_rows);
    int it = 0;
    double coefSum = 0;

    std::cout << __LINE__ << std::endl;
    // Compute probabilities
    for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.motif.n_rows; atomIndex++){
            coefs(it) = atomCoefficientSquared(atomIndex, cell, holeCell, RScoefs);
            coefSum += coefs(it);
            it++;
        }
    }
    std::cout << __LINE__ << std::endl;

    // Write probabilities to file
    it = 0;
    for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for(unsigned int atomIndex = 0; atomIndex < exciton.motif.n_rows; atomIndex++){
            arma::rowvec position = exciton.motif.row(atomIndex) + cell;
            fprintf(textfile, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            position(0), position(1), coefs(it)/coefSum);
        }
    }                                    

}

void Result::writeRealspaceAmplitude(int stateindex, int holeIndex, 
                                     const arma::rowvec& holeCell, FILE* textfile){

    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeRealspaceAmplitude(statecoefs, holeIndex, holeCell, textfile);

}

/* Method to write the eigenvalues into a file. Requires pointer to file, and has optional
argument n to specify the number of eigenvalues to write (ascending order). */
void Result::writeEigenvalues(FILE* textfile, int n){

    if(n > exciton.ncell || n < 0){
        throw std::invalid_argument("Optional argument n must be a positive integer equal or below basisdim");
    }

    fprintf(textfile, "%d\t", exciton.excitonbasisdim);
    int maxEigval = (n == 0) ? exciton.excitonbasisdim : n;  
    for(unsigned int i = 0; i < maxEigval; i++){
        fprintf(textfile, "%11.7lf\t", eigval(i));
    }
    fprintf(textfile, "\n");
}

/* Method to write eigenstates into a file. Has optional argument n to specify the number of states
to write (in ascending energy order) */
void Result::writeStates(FILE* textfile, int n){
    if(n > exciton.ncell || n < 0){
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

/* Routine to calculate the coefficients of the exciton real-space 
wavefunction, each one associated to each atom of the lattice.
Here we use the vectorized implementation.
NB: this function requires having the basis predefined, meaning
that we need valence and conduction arrays.
Input: cx_mat BSE coefficiens, vec kpoints,
int holePos (1D representation of positions, i.e. index),
int Ncell, int N, int nEdgeStates.
Output: cx_vec real-space coefficients; already summed over k */
arma::cx_mat Result::RScoefficientCalc(const arma::cx_vec& statecoefs, int holeIndex){

    arma::vec eigValTB;
    arma::cx_mat eigVecTB;
    int kBlock = 0;
    int kBlockEnd = 0;

    int dimTB = exciton.basisdim;
    int nk = exciton.nk;
    // Matrix dimension is: #orbitals*#motif*#orbitals
    arma::cx_mat RScoefs = arma::zeros<arma::cx_mat>(dimTB*exciton.norbitals, nk);

    for (int i = 0; i < nk; i++){

        kBlockEnd += (exciton.excitonbasisdim/nk);
        
        eigVecTB = exciton.eigvecKStack.slice(i);
        arma::cx_mat c = exciton.eigvecKQStack.slice(i);

        /* Precaution when picking the coeffs. associated to holePos, since
        now our basis splits H in two blocks for spin*/
        arma::cx_mat v = eigVecTB.rows(exciton.norbitals*holeIndex, 
                                        exciton.norbitals*(holeIndex + 1) - 1);
        /*arma::cx_mat v2 = eigVecTB.rows(dimTB/2 + exciton.norbitals/2*holePos, 
                                        dimTB/2 + exciton.norbitals/2*(holePos+1) - 1);
        arma::cx_mat v = arma::join_cols(v1, v2);*/

        v = v.cols(0, exciton.valenceBands.n_elem - 1);
        c = c.cols(exciton.valenceBands.n_elem, c.n_cols - 1);

        //arma::cx_mat M = arma::kron(c, arma::conj(v));
        arma::cx_mat M = arma::kron(c, v);
        arma::cx_vec A = statecoefs(arma::span(kBlock, kBlockEnd - 1));
        //arma::cout << "M: " << M << arma::endl;
        //arma::cout << "A: " << A << arma::endl;

        arma::cx_vec coefsK = M*A;
        RScoefs.col(i) = coefsK;

        kBlock = kBlockEnd;
    };

    return RScoefs;
};



/* Instead of using the s hydrogenic wave functions, we are simply
going to calculate the coefficient corresponding to each atom, since
initially we assumed orthonormality (i.e. they are punctual).
It is already modulus squared. 
Input: int n (1D atom position representation). 
int cell (corresponding unit cell) 
cx_mat RScoefs (matrix of real space coefficients as calculated previously)
vec kpoints.
Output: coefficient probability on atom (n, cell)*/
double Result::atomCoefficientSquared(int n, const arma::rowvec& cell, const arma::rowvec& hCell, 
                              const arma::cx_mat& RScoefs){

    std::complex<double> coef = 0;
    double tCoef = 0;
    std::complex<double> imag(0,1);
    int orbitalblock = exciton.norbitals*exciton.norbitals;

    // We need to separate in spin blocks
    for(int orb = 0; orb < orbitalblock; orb++){        

        coef = 0;
        // Spin up
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += RScoefs.col(i)(orbitalblock*n + orb)*exp(imag * arma::dot(exciton.kpoints.row(i), cell -  
            hCell));
        };
        tCoef += abs(coef)*abs(coef);
         
        // Spin down
        /*coef = 0;
        for(int i = 0; i < (int)RScoefs.n_cols; i++){
            coef += (RScoefs.col(i)(RScoefs.n_rows/2 + orbitalblock*n + orb)*exp(imag * arma::dot(exciton.kpoints.row(i), cell - hCell)));
        };
        tCoef += abs(coef)*abs(coef); */
    };
    
    //tCoef = real(coef);
    // arma::cout << real(coef) << "---" << tCoef << arma::endl;

    return tCoef;
};

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

// Routine to determine the size of a box that contains the reciprocal unit cell as given
// by the BZ mesh. Intended to use with full BZ meshes. Returns half the side of the box.
double Result::boundingBoxBZ(){
    double max_x = arma::max(exciton.kpoints.col(0));
    double max_y = arma::max(exciton.kpoints.col(1));

    double value = (max_x > max_y) ? max_x : max_y;
    return value;
}

/* ------------------------------ UNUSED ROUTINES (TO BE REMOVED) ------------------------------*/

double Result::realSpaceWavefunction(GExciton& exciton, const arma::cx_vec& BSEcoefs,
                             int electronIndex, int holeIndex,
                             const arma::rowvec& eCell, const arma::rowvec& hCell){

    std::complex<double> imag(0, 1);
    std::complex<double> amplitude = 0;
    std::complex<double> totalAmplitude = 0;
    arma::cx_vec BSEsquared = arma::conj(BSEcoefs) % BSEcoefs;
    // std::cout << arma::max(BSEsquared) << std::endl;
    double maxval = std::real(arma::max(BSEsquared));

    /* for (auto i : BSEsquared){
        if (std::real(i) > maxval/10){
            //std::cout << i << std::endl;
        }
    }*/

    for (int i = 0; i < exciton.excitonbasisdim; i++){
        for (int j = 0; j < exciton.excitonbasisdim; j++){


            amplitude = 0;
            arma::irowvec basisState = exciton.basisStates.row(i);
            arma::irowvec basisState2 = exciton.basisStates.row(j);

            int v1 = exciton.bandToIndex[basisState(0)];
            int v2 = exciton.bandToIndex[basisState2(0)];
            int c1 = exciton.bandToIndex[basisState(1)];
            int c2 = exciton.bandToIndex[basisState2(1)];
            int kIndex = basisState(2);
            int k2Index = basisState2(2);
            arma::rowvec kpoint = exciton.kpoints.row(kIndex);
            arma::rowvec kpoint2 = exciton.kpoints.row(k2Index);

            amplitude += std::exp(imag*arma::dot(kpoint - kpoint2, eCell - hCell))*
                         std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(electronIndex))*
                         exciton.eigvecKStack.slice(kIndex).col(v1)(holeIndex)*
                         std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(holeIndex))*
                         exciton.eigvecKQStack.slice(k2Index).col(c2)(electronIndex);

            if (false){
                amplitude += std::exp(imag*arma::dot(exciton.Q, eCell - hCell))*
                         std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(electronIndex))*
                         exciton.eigvecKQStack.slice(k2Index).col(c2)(holeIndex)*
                         std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(holeIndex))*
                         exciton.eigvecKStack.slice(kIndex).col(v1)(electronIndex);
            
                amplitude += std::exp(-imag*arma::dot(exciton.Q, eCell - hCell))*
                            std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(electronIndex))*
                            exciton.eigvecKStack.slice(kIndex).col(v1)(holeIndex)*
                            std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(holeIndex))*
                            exciton.eigvecKQStack.slice(k2Index).col(c2)(electronIndex);
                
                amplitude += std::exp(-imag*arma::dot(kpoint - kpoint2, eCell - hCell))*
                            std::conj(exciton.eigvecKStack.slice(k2Index).col(v2)(electronIndex))*
                            exciton.eigvecKQStack.slice(k2Index).col(c2)(holeIndex)*
                            std::conj(exciton.eigvecKQStack.slice(kIndex).col(c1)(holeIndex))*
                            exciton.eigvecKStack.slice(kIndex).col(v1)(electronIndex);
            }

            amplitude *= std::conj(BSEcoefs(i))*BSEcoefs(j);
            totalAmplitude += amplitude;     
        }
    }

    //arma::cout << totalAmplitude << "---" << std::abs(totalAmplitude) << arma::endl;
    if (std::imag(totalAmplitude) > 1E-5){
        arma::cout << "Warning: Imaginary part of amplitude is not" << arma::endl;
    }
    return std::abs(totalAmplitude);
};


// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
std::complex<double> Result::densityMatrix(GExciton& exciton, const arma::cx_vec& BSEcoefs, 
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

// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
std::complex<double> Result::densityMatrixK(int kIndex, GExciton& exciton, const arma::cx_vec& BSEcoefs, 
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