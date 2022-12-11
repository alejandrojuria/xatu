#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "xatu/System.hpp"

namespace xatu {

// -------------------- Constructors and destructor --------------------

// Default constructor
System::System() : Crystal(){};

// Copy constructor
System::System(const System& system) : Crystal(system){

	systemName			 = system.systemName;
	orbitals_            = system.orbitals;
	hamiltonianMatrices  = system.hamiltonianMatrices;
	overlapMatrices      = system.overlapMatrices;
	filling_			 = system.filling;
	fermiLevel_			 = filling_ - 1;
	basisdim_ 			 = system.basisdim;

	std::cout << "Correctly initiallized System object" << std::endl;
}

// Configuration constructor
System::System(const SystemConfiguration& configuration) : Crystal(){

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
	
	systemName			 = configuration.systemInfo.name;
	orbitals_            = configuration.systemInfo.norbitals;
	hamiltonianMatrices  = configuration.systemInfo.hamiltonian;
	overlapMatrices      = configuration.systemInfo.overlap;
	filling_			 = configuration.systemInfo.filling;
	fermiLevel_			 = filling_ - 1;

    int basisdim = 0;
    for(int i = 0; i < natoms; i++){
        int species = this->motif.row(i)(3);
        basisdim += orbitals(species);
    }
	basisdim_   = basisdim;
}

/* Bloch hamiltonian for posterior diagonalization. Input: arma::vec k (wave number) */
arma::cx_mat System::hamiltonian(arma::rowvec k, bool isTriangular){

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

/* Overlap matrix in reciprocal space to solve generalized eigenvalue problem */
arma::cx_mat System::overlap(arma::rowvec k, bool isTriangular){

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

void System::setFilling(int filling){
	if (filling > 0){
		filling_ = filling;
		fermiLevel_ = filling_;
	}
	else{
		std::cout << "Filling must be a positive integer" << std::endl;
	}
}

void System::solveBands(arma::rowvec& k, arma::vec& eigval, arma::cx_mat& eigvec, bool triangular){
	arma::cx_mat h = hamiltonian(k, triangular);
	double auToEV = 27.2;
	if (!overlapMatrices.empty()){
		arma::cx_vec auxEigval;
		arma::cx_mat s = overlap(k, triangular);
		arma::eig_pair(auxEigval, eigvec, h, s);
		arma::uvec sorted_indices = arma::sort_index(auxEigval);
		eigval = arma::sort(arma::real(auxEigval)*auToEV);
		eigvec = eigvec.cols(sorted_indices);
		for (int i = 0; i < eigvec.n_cols; i++){
			arma::cx_vec state = eigvec.col(i);
			double norm = real(arma::cdot(state, s*state));
			eigvec.col(i) /= std::sqrt(norm);
		}
		orthogonalize(k, eigvec, triangular);		
	}
	else{
		arma::eig_sym(eigval, eigvec, h);
	}
	
}

void System::solveBands(std::string kpointsfile, bool triangular){
	std::ifstream inputfile;
	std::string line;
	double kx, ky, kz;
	arma::vec eigval;
	arma::cx_mat eigvec;
	std::string outputfilename = systemName + ".bands";
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
}

void System::orthogonalize(const arma::rowvec& k, arma::cx_mat& states, bool triangular){
	// First compute X
	arma::cx_mat s = overlap(k, triangular);
	arma::vec eigval;
	arma::cx_mat eigvec;
	arma::eig_sym(eigval, eigvec, s);

	eigval = 1./arma::sqrt(eigval);
	arma::cx_mat sRoot = arma::zeros<arma::cx_mat>(eigval.n_elem, eigval.n_elem);
	sRoot.diag() = arma::conv_to<arma::cx_vec>::from(eigval);
	sRoot = eigvec*sRoot*eigvec.t();

	states = arma::inv_sympd(sRoot) * states;	
}

// /*------------------ Utilities/Additional observables ------------------*/

// BEWARE: Basis order for spin value is NOT correct
/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinZValue(const arma::cx_vec& eigvec){

	arma::cx_vec spinEigvalues = {1./2, -1./2};
    arma::cx_vec spinVector = arma::kron(arma::ones(basisdim/2), spinEigvalues);
	arma::cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));
};

// TODO: ADD CALCULATION OF X, Y SPIN COMPONENTS
/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinYValue(const arma::cx_vec& eigvec){

	std::complex<double> i(0,1);
	arma::cx_mat operatorSy = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	operatorSy(0,1) = -i;
	operatorSy(1,0) = i;
	operatorSy.submat(2,5, 4,7) = -i*arma::eye<arma::cx_mat>(3,3);
	operatorSy.submat(5,2, 7,4) = i*arma::eye<arma::cx_mat>(3,3);
	
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double System::expectedSpinXValue(const arma::cx_vec& eigvec){

	arma::cx_mat operatorSx = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	operatorSx(0,1) = 1;
	operatorSx(1,0) = 1;
	operatorSx.submat(2,5, 4,7) = arma::eye<arma::cx_mat>(3,3);
	operatorSx.submat(5,2, 7,4) = arma::eye<arma::cx_mat>(3,3);

};

}