#include "Result.hpp"
#include <complex>

Result::Result(GExciton& exciton, arma::vec& eigval, arma::cx_mat& eigvec) : 
    exciton(exciton), eigval(eigval), eigvec(eigvec){};

double Result::kineticEnergy(int stateindex){
    std::complex<double> energy = arma::cdot(eigvec, exciton.HK*eigvec);
    return energy.real();
}

double Result::potentialEnergy(int stateindex){
    arma::cx_mat HV = exciton.HBS - exciton.HK;
    std::complex<double> energy = arma::cdot(eigvec, HV*eigvec);
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


void Result::writeReciprocalAmplitude(int stateindex, FILE* textfile){
    fprintf(textfile, "kx\tky\tkz\tProb.\n");
    arma::cx_vec state = eigvec.col(stateindex);
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    for (int i = 0; i < exciton.kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(state(nbandsCombinations*i + nband))*abs(state(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.kpoints.row(1) - exciton.kpoints.row(0)); // L2 norm instead of l2
        fprintf(textfile, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                    exciton.kpoints.row(i)(0), exciton.kpoints.row(i)(1), exciton.kpoints.row(i)(2), coef);
    };
    fprintf(textfile, "#\n");
};

void Result::writeExtendedReciprocalAmplitude(int stateindex, FILE* textfile){

}

void Result::writeRealspaceAmplitude(int stateindex, FILE* textfile){

}