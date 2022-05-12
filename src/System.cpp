#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "System.hpp"


// -------------------- Constructor and destructor --------------------

System::System(std::string filename) : Crystal(){

	SystemConfiguration configuration = SystemConfiguration(filename);
	initializeCrystalAttributes(configuration);
	initializeSystemAttributes(configuration);
	
	std::cout << "Correctly initiallized System object" << std::endl;
};

System::~System(){
	std::cout << "Destroying System object..." << std::endl;
};

// ----------------------------- Methods -----------------------------
/* Routine to extract the information contained in the SystemConfiguration object from
parsing the input text file */
void System::initializeSystemAttributes(const SystemConfiguration& configuration){
	
	orbitals            = configuration.systemInfo.norbitals;
	hamiltonianMatrices = configuration.systemInfo.hamiltonian;
	overlapMatrices     = configuration.systemInfo.overlap;
	filling_			= configuration.systemInfo.filling;

	norbitals_  = arma::sum(orbitals);
	basisdim_   = norbitals*natoms; // Wrong; valid only for single species solid
	fermiLevel_ = (int)(filling * basisdim) - 1;
}

/* Bloch hamiltonian for posterior diagonalization. Input: arma::vec k (wave number) */
arma::cx_mat System::hamiltonian(arma::rowvec k, bool isTriangular){

	arma::cx_mat h = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		h += hamiltonianMatrices.slice(i) * std::exp(imag*arma::dot(k, cell));
	};

	h += h.t();

	return h;
};

/* Overlap matrix in reciprocal space to solve generalized eigenvalue problem */
arma::cx_mat System::overlap(arma::rowvec k, bool isTriangular){

	arma::cx_mat s = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		s += overlapMatrices.slice(i) * std::exp(imag*arma::dot(k, cell));
	};

	s += s.t();
	return s;
}

void System::setFilling(double filling){
	if (filling > 0 && filling < 1){
		filling_ = filling;
		fermiLevel_ = (int)(filling * norbitals);
	}
	else{
		std::cout << "Filling must be between 0 and 1" << std::endl;
	}
}

// /*------------------ Utilities/Additional observables ------------------*/

// BEWARE: Basis order for spin value is NOT correct
/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinZValue(const arma::cx_vec& eigvec){

	arma::cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
	arma::cx_vec spinVector = arma::kron(arma::ones(basisdim/norbitals, 1), spinEigvalues);

	arma::cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinYValue(const arma::cx_vec& eigvec){

	std::complex<double> i(0,1);
	arma::cx_mat operatorSy = arma::zeros<arma::cx_mat>(8, 8);
	operatorSy(0,1) = -i;
	operatorSy(1,0) = i;
	operatorSy.submat(2,5, 4,7) = -i*arma::eye<arma::cx_mat>(3,3);
	operatorSy.submat(5,2, 7,4) = i*arma::eye<arma::cx_mat>(3,3);

	arma::cx_mat syMatrix = 0.5*arma::kron(arma::eye<arma::cx_mat>(basisdim/norbitals, 
                                           norbitals/norbitals), operatorSy);

	arma::cx_vec spinEigvec = syMatrix*eigvec;
	
	return abs(arma::cdot(eigvec, spinEigvec));
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinXValue(const arma::cx_vec& eigvec){

	arma::cx_mat operatorSx = arma::zeros<arma::cx_mat>(8, 8);
	operatorSx(0,1) = 1;
	operatorSx(1,0) = 1;
	operatorSx.submat(2,5, 4,7) = arma::eye<arma::cx_mat>(3,3);
	operatorSx.submat(5,2, 7,4) = arma::eye<arma::cx_mat>(3,3);

	arma::cx_mat sxMatrix = 0.5*arma::kron(arma::eye<arma::cx_mat>(basisdim/norbitals, basisdim/norbitals), 
                                                             operatorSx);

	arma::cx_vec spinEigvec = sxMatrix*eigvec;
	
	return abs(arma::cdot(eigvec, spinEigvec));
};