#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "xatu/System.hpp"

namespace xatu {

// -------------------- Constructors and destructor --------------------

/**
 * Default constructor.
 */
System::System() : Crystal(){};

/**
 * Copy constructor.
 * @details Initializes one default System objects, and afterwards copies
 * its attributes to be the same as those of the given System object to copy.
 * @param system System object to copy.
 */ 
System::System(const System& system) : Crystal(system){

	systemName			 = system.systemName;
	orbitals_            = system.orbitals;
	hamiltonianMatrices_ = system.hamiltonianMatrices;
	overlapMatrices_     = system.overlapMatrices;
	filling_			 = system.filling;
	fermiLevel_			 = filling_ - 1;
	basisdim_ 			 = system.basisdim;
}

/**
 * Configuration constructor.
 * @details Constructor which takes in a SystemConfiguration object, i.e.
 * to init a System from a configuration file.
 * @param configuration SystemConfiguration object obtained from config. file.
 */
System::System(const SystemConfiguration& configuration) : Crystal(){

	initializeCrystalAttributes(configuration);
	initializeSystemAttributes(configuration);
};

// --------------------------------- Methods ---------------------------------
/**
 * Routine to extract the information contained in the SystemConfiguration object from
 * parsing the input text file.
 * @details To be called from the configuration constructor.
 * @param configuration SystemConfiguration object.
 */
void System::initializeSystemAttributes(const SystemConfiguration& configuration){
	
	orbitals_            = configuration.systemInfo.norbitals;
	hamiltonianMatrices_ = configuration.systemInfo.hamiltonian;
	overlapMatrices_     = configuration.systemInfo.overlap;
	filling_			 = configuration.systemInfo.filling;
	fermiLevel_			 = filling_ - 1;

    int basisdim = 0;
    for(int i = 0; i < natoms; i++){
        int species = this->motif.row(i)(3);
        basisdim += orbitals(species);
    }
	basisdim_   = basisdim;
}

/**
 * Bloch Hamiltonian.
 * @details The Bloch Hamiltonian is constructed at each k point using the Fock
 * matrices of the system, H(R). The Fock matrices can be either complete or triangular.
 * @param k kpoint where we want to evaluate H(k).
 * @param isTriangular Whether the Fock matrices are triangular or complete.
 * @returns Bloch Hamiltonian matrix at k.
 */
arma::cx_mat System::hamiltonian(arma::rowvec k, bool isTriangular) const{

	arma::cx_mat h = arma::zeros<arma::cx_mat>(basisdim, basisdim);
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
arma::cx_mat System::overlap(arma::rowvec k, bool isTriangular) const{

	arma::cx_mat s = arma::zeros<arma::cx_mat>(basisdim, basisdim);
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
 * Filling setter.
 * @details Sets both the filling and the Fermi level, which is defined as filling - 1.
 * Filling must be a positive integer.
 * @param filling Number of electrons of the system.
*/
void System::setFilling(int filling){
	if (filling > 0){
		filling_ = filling;
		fermiLevel_ = filling_ - 1;
	}
	else{
		std::cout << "Filling must be a positive integer" << std::endl;
	}
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
void System::solveBands(arma::rowvec& k, arma::vec& eigval, arma::cx_mat& eigvec, bool triangular) const {
	arma::cx_mat h = hamiltonian(k, triangular);
	double auToEV = 27.2;

	if (!overlapMatrices.empty()){
		h *= auToEV;
		orthogonalize_hamiltonian(k, h, triangular);	
	}

	arma::eig_sym(eigval, eigvec, h);
}

/**
 * Method to write to a file the energy bands evaluated on a set of kpoints specified on a file.
 * @details It creates a file with the name "[systemName].bands" where the bands are stores.
 * @param kpointsfile File with the kpoints where we want to obtain the bands.
 * @param triangular To specify if the Hamiltonian is triangular.
*/
void System::solveBands(std::string kpointsfile, bool triangular) const {
	std::ifstream inputfile;
	std::string line;
	double kx, ky, kz;
	arma::vec eigval;
	arma::cx_mat eigvec;
	std::string outputfilename = kpointsfile + ".bands";
	FILE* bandfile = fopen(outputfilename.c_str(), "w");
	try{
		inputfile.open(kpointsfile.c_str());
		while(std::getline(inputfile, line)){
			std::istringstream iss(line);
			iss >> kx >> ky >> kz;
			arma::rowvec kpoint{kx, ky, kz};
			solveBands(kpoint, eigval, eigvec, triangular);
			for (int i = 0; i < eigval.n_elem; i++){
				fprintf(bandfile, "%12.6f\t", eigval(i));
			}
			fprintf(bandfile, "\n");
		}
	}
	catch(const std::exception& e){
		std::cerr << e.what() << std::endl;
	}
	fclose(bandfile);
	arma::cout << "Done" << arma::endl;
}

/**
 * Method to orthogonalize the basis.
 * @details This method acts directly over the eigenstates of the system, giving us
 * their coefficients if their were written in an orthonormal basis. The technique used
 * is Lowdin orthogonalization.
 * @param k kpoint of eigenvectors to orthonormalize.
 * @param states Matrix of eigenstates.
 * @param Triangular To specify if the Hamiltonian is triangular.
*/
void System::orthogonalize(const arma::rowvec& k, arma::cx_mat& states, bool triangular) const {
	// First compute X
	arma::cx_mat s = overlap(k, triangular);
	arma::vec eigval;
	arma::cx_mat eigvec;
	arma::eig_sym(eigval, eigvec, s);

	eigval = arma::sqrt(eigval);
	arma::cx_mat sRoot = arma::zeros<arma::cx_mat>(eigval.n_elem, eigval.n_elem);
	sRoot.diag() = arma::conv_to<arma::cx_vec>::from(eigval);
	sRoot = eigvec*sRoot*eigvec.t();

	// // states = arma::inv_sympd(eigvec) * states;
	states = sRoot * states;
	for(int i = 0; i < states.n_cols; i++){
		states.col(i) = arma::normalise(states.col(i));
	}
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
void System::orthogonalize_hamiltonian(const arma::rowvec& k, arma::cx_mat& hamiltonian, bool triangular) const {
	// First compute X
	arma::cx_mat s = overlap(k, triangular);
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

	// // states = arma::inv_sympd(eigvec) * states;
	hamiltonian = sRoot * hamiltonian * sRoot;
}

// /*------------------ Utilities/Additional observables ------------------*/

/** 
 * Routine to calculate the expected value of the spin projection Sz of an eigenstate.
 * @details Note that this assumes a certain basis ordering of the system to produce
 * a correct spin value. The ordering assumed is {|1,up>, |1,down>, |2,up>, ...}, i.e.
 * alternating spin consecutively.
 * @param eigvec State.
 * @returns Expectation value of Sz.
 * */
double System::expectedSpinZValue(const arma::cx_vec& eigvec){

	arma::cx_vec spinEigvalues = {1./2, -1./2};
    arma::cx_vec spinVector = arma::kron(arma::ones(basisdim/2), spinEigvalues);
	arma::cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));
};

/**
 * Routine to calculate the expected value of the spin component Sy.
 * NB: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sy.
 */
double System::expectedSpinYValue(const arma::cx_vec& eigvec){

	std::complex<double> i(0,1);
	arma::cx_mat operatorSy = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	operatorSy(0,1) = -i;
	operatorSy(1,0) = i;
	operatorSy.submat(2,5, 4,7) = -i*arma::eye<arma::cx_mat>(3,3);
	operatorSy.submat(5,2, 7,4) = i*arma::eye<arma::cx_mat>(3,3);
	
};

/**
 * Routine to calculate the expected value of the spin component Sx.
 * NB: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sx.
 */
double System::expectedSpinXValue(const arma::cx_vec& eigvec){

	arma::cx_mat operatorSx = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	operatorSx(0,1) = 1;
	operatorSx(1,0) = 1;
	operatorSx.submat(2,5, 4,7) = arma::eye<arma::cx_mat>(3,3);
	operatorSx.submat(5,2, 7,4) = arma::eye<arma::cx_mat>(3,3);

};

}