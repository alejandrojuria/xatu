#include <complex>
#include <math.h>
#include <stdlib.h>
#include "xatu/SystemTB.hpp"
#include "xatu/utils.hpp"

namespace xatu {

// -------------------- Constructors and destructor --------------------

/**
 * Copy constructor.
 * @details Initializes one default System objects, and afterwards copies
 * its attributes to be the same as those of the given SystemTB object to copy.
 * @param system SystemTB object to copy.
 */ 
SystemTB::SystemTB(const SystemTB& system) : System(system) {

	isTriangular_ = system.isTriangular;
	isAU_     	  = system.isAU;
}

/**
 * Configuration constructor.
 * @details Constructor which takes in a SystemConfiguration object, i.e.
 * to init a SystemTB from a configuration file.
 * @param configuration SystemConfiguration object obtained from config. file.
 */
SystemTB::SystemTB(const SystemConfiguration& configuration) : System(configuration){
	this->isTriangular_ = checkIfTriangular(this->hamiltonianMatrices.slice(0));
};

// --------------------------------- Methods ---------------------------------

/**
 * Method to set the triangularity of the Fock matrices.
 * @param isTriangular Whether the Fock matrices are triangular or not. 
 */
void SystemTB::setTriangular(bool isTriangular){
	isTriangular_ = isTriangular;
};

/**
 * Method to set the type of system (CRYSTAL or not).
 * @details The isCRYSTAL attribute ensures proper handling of the energies (to be given in eV).
 * @param isCRYSTAL Whether the system is CRYSTAL or not.
 */
void SystemTB::setAU(bool isAU){
	isAU_ = isAU;
};

/**
 * Bloch Hamiltonian.
 * @details The Bloch Hamiltonian is constructed at each k point using the Fock
 * matrices of the system, H(R). The Fock matrices can be either complete or triangular.
 * @param k kpoint where we want to evaluate H(k).
 * @param isTriangular Whether the Fock matrices are triangular or complete.
 * @returns Bloch Hamiltonian matrix at k.
 */
arma::cx_mat SystemTB::hamiltonian(arma::rowvec k) const{

	arma::cx_mat h = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		h += hamiltonianMatrices.slice(i) * std::exp(imag*arma::dot(k, cell));
	};

	if (isTriangular){
		h.diag() -= h.diag()/2;
		h += h.t();
	}

	return h;
};

/**
 * Overlap matrix in reciprocal space.
 * @details The reciprocal overlap matrix S(k) is required to solve the generalized
 * eigenvalue problem, which appears with non-orthonormal basis. It is built from the
 * overlap matrices in real space S(R).
 * @param k kpoint where we want to evaluate S(k).
 * @param isTriangular To specify whether the overlap matrices are triangular.
 * @returns Reciprocal overlap matrix S(k).
 */
arma::cx_mat SystemTB::overlap(arma::rowvec k) const{

	arma::cx_mat s = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		s += overlapMatrices.slice(i) * std::exp(imag*arma::dot(k, cell));
	};

	if (isTriangular){
		s.diag() -= s.diag()/2;
		s += s.t();
	}

	return s;
}

/**
 * Method to obtain the energy bands and eigenvectors at a given k.
 * @details Depending on whether the Fock and overlap matrices are tringular
 * or not (if they are triangular also it is assumed that the basis is non-orthonormal),
 * the returning eigenvectors are orthonormalized.
 * @param k kpoint where we want the eigenvalues and eigenvectors.
 * @param eigval Vector to store the energies of the system.
 * @param eigvec Complex matrix to store eigenvectors.
 * @param triangular To specifiy if the matrices are triangular.
*/
void SystemTB::solveBands(arma::rowvec& k, arma::vec& eigval, arma::cx_mat& eigvec) const {
	arma::cx_mat h = hamiltonian(k);
	double auToEV = 27.2;

	if (!overlapMatrices.empty()){
		if (isAU){
			h *= auToEV;
		}
		orthogonalize_hamiltonian(k, h);	
	}

	arma::eig_sym(eigval, eigvec, h);
}

/**
 * Method to orthogonalize the basis.
 * @details This method acts directly over the eigenstates of the system, giving us
 * their coefficients if their were written in an orthonormal basis. The technique used
 * is Lowdin orthogonalization.
 * @param k kpoint of eigenvectors to orthonormalize.
 * @param hamiltonian Matrix of eigenstates.
 * @param Triangular To specify if the Hamiltonian is triangular.
*/
void SystemTB::orthogonalize_hamiltonian(const arma::rowvec& k, arma::cx_mat& hamiltonian) const {
	// First compute X
	arma::cx_mat s = overlap(k);
	arma::vec eigval;
	arma::cx_mat eigvec;
	arma::eig_sym(eigval, eigvec, s);

	try{
		eigval = 1./arma::sqrt(eigval);
	}
	catch(const std::exception& e){
		throw std::invalid_argument("Zero or negative overlap eigenvalues found, exiting...");
	}
	
	arma::cx_mat sRoot = arma::zeros<arma::cx_mat>(eigval.n_elem, eigval.n_elem);
	sRoot.diag() = arma::conv_to<arma::cx_vec>::from(eigval);
	sRoot = eigvec*sRoot*eigvec.t();

	// states = arma::inv_sympd(eigvec) * states;
	hamiltonian = sRoot * hamiltonian * sRoot;
}

// /*------------------ Utilities/Additional observables ------------------*/

/** 
 * Routine to calculate the coefficients corresponding to wavefunctions in the atomic gauge.
 * @param coefs Lattice gauge state coefficients on which we perform the gauge transformation.
 * @param k kpoint required to perform the transformation.
 * @return Atomic gauge state coefficients.
 */
arma::cx_vec SystemTB::latticeToAtomicGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

    arma::cx_vec phases(norbitals);
    std::complex<double> imag(0, 1);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.row(atomIndex)(3);
        for(int orbIndex = 0; orbIndex < orbitalsPerSpecies(species); orbIndex++){
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
arma::cx_vec SystemTB::atomicToLatticeGauge(const arma::cx_vec& coefs, const arma::rowvec& k){

    arma::cx_vec phases(norbitals);
    std::complex<double> imag(0, 1);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.row(atomIndex)(3);
        for(int orbIndex = 0; orbIndex < orbitalsPerSpecies(species); orbIndex++){
            arma::rowvec atomPosition = motif.row(atomIndex).subvec(0, 2);
            phases(it) = exp(-imag*arma::dot(k, atomPosition));
            it++;
        }
    }

    arma::cx_vec atomicCoefs = coefs % phases;
    return atomicCoefs;
}

/** 
 * Routine to calculate the expected value of the spin projection Sz of an eigenstate.
 * @details Note that this assumes a certain basis ordering of the system to produce
 * a correct spin value. The ordering assumed is {|1,up>, |1,down>, |2,up>, ...}, i.e.
 * alternating spin consecutively.
 * @param eigvec State.
 * @returns Expectation value of Sz.
 * */
double SystemTB::expectedSpinZValue(const arma::cx_vec& eigvec){

	arma::cx_vec spinEigvalues = {1./2, -1./2};
    arma::cx_vec spinVector = arma::kron(arma::ones(norbitals/2), spinEigvalues);
	arma::cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));
};

/**
 * Routine to calculate the expected value of the spin component Sy.
 * TODO: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sy.
 */
double SystemTB::expectedSpinYValue(const arma::cx_vec& eigvec){

	std::complex<double> i(0,1);
	arma::cx_mat operatorSy = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	operatorSy(0,1) = -i;
	operatorSy(1,0) = i;
	operatorSy.submat(2,5, 4,7) = -i*arma::eye<arma::cx_mat>(3,3);
	operatorSy.submat(5,2, 7,4) = i*arma::eye<arma::cx_mat>(3,3);
	
};

/**
 * Routine to calculate the expected value of the spin component Sx.
 * TODO: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sx.
 */
double SystemTB::expectedSpinXValue(const arma::cx_vec& eigvec){

	arma::cx_mat operatorSx = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	operatorSx(0,1) = 1;
	operatorSx(1,0) = 1;
	operatorSx.submat(2,5, 4,7) = arma::eye<arma::cx_mat>(3,3);
	operatorSx.submat(5,2, 7,4) = arma::eye<arma::cx_mat>(3,3);

};

/**
 * Routine to calculate the velocity matrix element between two bands.
 * @details This routine is intended to work with tight-binding systems only, and will produce
 * incorrect results otherwise.
 * TODO: Duplicity with velocitySingleParticle in the Result class, merge both routines.
 * @param k kpoint where we want to evaluate the velocity.
 * @param fBand Index of the first (row) band.
 * @param sBand Index of the second (column) band.
*/
arma::cx_vec SystemTB::velocity(const arma::rowvec k, int fBand, int sBand) const {
	arma::cx_mat h = hamiltonian(k);
	arma::vec eigval;
	arma::cx_mat eigvec;
	arma::eig_sym(eigval, eigvec, h);

	arma::cx_vec fBandEigvec = eigvec.col(fBand);
	arma::cx_vec sBandEigvec = eigvec.col(sBand);

    arma::cx_cube hkDerivative = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);
    arma::cx_cube iHt = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);

    std::complex<double> imag(0, 1);

    // First compute Hk derivative
    for (int j = 0; j < 3; j++){
        for (int i = 0; i < ncells; i++){
            arma::rowvec cell = unitCellList.row(i);
            hkDerivative.slice(j) += hamiltonianMatrices.slice(i) * 
                                     std::exp(imag*arma::dot(k, cell)) * cell(j) * imag;
	    };
    }

    // Next compute iH(t-t') matrix
    arma::cx_cube motifDifference = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);
    arma::cx_mat extendedMotif = arma::zeros<arma::cx_mat>(norbitals, 3);
    int currentIndex = 0;
    for (int i = 0; i < natoms; i++){
        int norb = orbitalsPerSpecies(motif.row(i)(3));
        extendedMotif.rows(currentIndex, currentIndex + norb - 1) = arma::kron(motif.row(i).subvec(0, 2),
                                                                         arma::ones<arma::cx_vec>(norb));
        currentIndex += norb;
    }

    arma::cx_mat blochHamiltonian = hamiltonian(k);
    for (int j = 0; j < 3; j++){
        motifDifference.slice(j) = arma::kron(extendedMotif.col(j), arma::ones<arma::cx_rowvec>(norbitals)) -
                                   arma::kron(extendedMotif.col(j).t(), arma::ones<arma::cx_vec>(norbitals));
        iHt.slice(j) = imag * blochHamiltonian % motifDifference.slice(j).t();
    }

    // Finally compute velocity matrix elements
    arma::cx_vec velocityMatrixElement = arma::zeros<arma::cx_vec>(3);
	for (int j = 0; j < 3; j++){
		velocityMatrixElement(j) = arma::cdot(fBandEigvec, (hkDerivative.slice(j) + iHt.slice(j)) * sBandEigvec);
	}

    return velocityMatrixElement;
};

}