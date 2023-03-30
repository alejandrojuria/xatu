#include "xatu/Result.hpp"
#include <complex>

namespace xatu {


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
    int valence = exciton.bandToIndex[exciton.valenceBands.max()];
    int conduction = exciton.bandToIndex[exciton.conductionBands.min()];

    double gap = exciton.eigvalKStack.col(kIndex)(conduction) - 
                 exciton.eigvalKStack.col(kIndex)(valence);
    return gap;
}

// Routine to find the k index where the exciton has the maximum amplitude, which
// usually corresponds with the band gap location
int Result::findExcitonPeak(int stateindex){
    int index = eigvec.col(stateindex).index_max();
    int bandCombinations = exciton.valenceBands.n_elem*exciton.conductionBands.n_elem;
    index = (int)index/bandCombinations;
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

void Result::writeRealspaceAmplitude(int stateindex, int holeIndex, 
                                     const arma::rowvec& holeCell, FILE* textfile, int ncells){

    arma::cx_vec statecoefs = eigvec.col(stateindex);
    writeRealspaceAmplitude(statecoefs, holeIndex, holeCell, textfile, ncells);
}

/* Method to write the eigenvalues into a file. Requires pointer to file, and has optional
argument n to specify the number of eigenvalues to write (ascending order). */
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

/* Method to write eigenstates into a file. Has optional argument n to specify the number of states
to write (in ascending energy order) */
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

    skubo_w_(&nR, &norb, &norb_ex, &nv, &nc, &filling, 
             Rvec.memptr(), R.memptr(), extendedMotif.memptr(), hhop.memptr(), shop.memptr(), &nk, rkx.memptr(),
             rky.memptr(), rkz.memptr(), m_eigvec.memptr(), m_eigval.memptr());
}

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

arma::cx_vec Result::addExponential(arma::cx_vec& coefs, const arma::rowvec& cell){

    arma::vec product = exciton.kpoints * cell.t();
    std::complex<double> imag(0, 1);
    arma::cx_vec exponentials = arma::exp(imag*product);
    int nBandCombinations = exciton.valenceBands.n_elem*exciton.conductionBands.n_elem;
    exponentials = arma::kron(exponentials, arma::ones<arma::cx_vec>(nBandCombinations));

    coefs = coefs % exponentials;

    return coefs;
}


// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
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

// This routine requires having ALL eigenvectors from the Bloch Hamiltonian, otherwise it is not possible to compute.
// Meaning bool storeAllVectors must be set to TRUE when initializing the exciton object.
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