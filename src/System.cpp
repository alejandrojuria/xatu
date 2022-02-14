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

	norbitals_ = arma::sum(orbitals);
	basisdim_  = norbitals*natoms;
}

/* Bloch hamiltonian for posterior diagonalization. Input: arma::vec k (wave number) */
arma::cx_mat System::hamiltonian(arma::rowvec k, bool isTriangular){

	arma::cx_mat h = arma::zeros<arma::cx_mat>(basisdim, basisdim);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		h += hamiltonianMatrices.slice(i) * std::exp(-imag*arma::dot(k, cell));
	};

	if (isTriangular){
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
		s += overlapMatrices.slice(i) * std::exp(-imag*arma::dot(k, cell));
	};

	s += s.t();
	return s;
}

// /*------------------ Utilities/Additional observables ------------------*/

// /* Routine to apply the inversion operator P over a eigenstate */
// cx_mat System::inversionOperator(const cx_vec& eigenvector){

// 	int dimTB = 2*(N+1)*8;
// 	cx_mat P = arma::zeros<cx_mat>(dimTB/8, dimTB/8);
// 	for(int i = 0; i < dimTB/8; i++){
// 		P(dimTB/8 - 1 - i, i) = 1.;
// 	};
// 	cx_mat spinOperator = arma::zeros<cx_mat>(8,8);
// 	spinOperator(0, 1) = -1;
// 	spinOperator(1, 0) = -1;
// 	spinOperator.submat(2,5, 4,7) = -arma::eye<cx_mat>(3,3);
// 	spinOperator.submat(5,2, 7,4) = -arma::eye<cx_mat>(3,3);
// 	P = arma::kron(P, spinOperator);

// 	return P*eigenvector;
// };

// /* Routine to calculate the expected value of the spin projection
// Sz of a TB eigenstate. */
// double System::expectedSpinZValue(const arma::cx_vec& eigvec){

// 	int dimTB = 2*(N+1)*8;
// 	cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
// 	cx_vec spinVector = arma::kron(arma::ones(dimTB/8, 1), spinEigvalues);

// 	cx_vec spinEigvec = eigvec % spinVector;

// 	return real(arma::cdot(eigvec, spinEigvec));
// };

// /* Routine to calculate the expected value of the spin projection
// Sz of a TB eigenstate. */
// double System::expectedSpinYValue(const arma::cx_vec& eigvec){

// 	int dimTB = 2*(N+1)*8;
// 	std::complex<double> i(0,1);
// 	cx_mat operatorSy = arma::zeros<cx_mat>(8, 8);
// 	operatorSy(0,1) = -i;
// 	operatorSy(1,0) = i;
// 	operatorSy.submat(2,5, 4,7) = -i*arma::eye<cx_mat>(3,3);
// 	operatorSy.submat(5,2, 7,4) = i*arma::eye<cx_mat>(3,3);

// 	cx_mat syMatrix = 0.5*arma::kron(arma::eye<cx_mat>(dimTB/8, dimTB/8), operatorSy);

// 	cx_vec spinEigvec = syMatrix*eigvec;
	
// 	return abs(arma::cdot(eigvec, spinEigvec));
// };

// /* Routine to calculate the expected value of the spin projection
// Sz of a TB eigenstate. */
// double System::expectedSpinXValue(const arma::cx_vec& eigvec){

// 	int dimTB = 2*(N+1)*8;
// 	cx_mat operatorSx = arma::zeros<cx_mat>(8, 8);
// 	operatorSx(0,1) = 1;
// 	operatorSx(1,0) = 1;
// 	operatorSx.submat(2,5, 4,7) = arma::eye<cx_mat>(3,3);
// 	operatorSx.submat(5,2, 7,4) = arma::eye<cx_mat>(3,3);

// 	cx_mat sxMatrix = 0.5*arma::kron(arma::eye<cx_mat>(dimTB/8, dimTB/8), operatorSx);

// 	cx_vec spinEigvec = sxMatrix*eigvec;
	
// 	return abs(arma::cdot(eigvec, spinEigvec));
// };


/* Definition of non-interacting retarded Green function */
std::complex<double> System::rGreenF(double energy, double delta, double eigEn){

	std::complex<double> i(0,1);
	return 1./((energy + i*delta) - eigEn);
};

// /* Routine to calcule the density of states at a given energy,
// associated to a given set of eigenvalues (e.g. bulk or edge).
// NB: It is NOT normalized. */
// double System::densityOfStates(double energy, double delta, const mat& energies){

//  	double dos = 0;
//  	for(int i = 0; i < (int)energies.n_rows; i++){
//  		for(int j = 0; j < (int)energies.n_cols; j++){
//  			double eigEn = energies(i,j);
//  			dos += -PI*imag(rGreenF(energy, delta, eigEn));
//  		};
//  	};
//  	dos /= energies.n_cols*a; // Divide by number of k's and length a

//  	return dos;
// }

// /* Routine to calculate and write the density of states associated 
// to a given set of eigenenergies. Density of states is normalized to 1
// (integral of DOS over energies equals 1)
// Input: mat energies (eigenvalues), double delta (convergence parameter),
// FILE* dosfile output file
// Output: void (write results to output file) */
// void System::writeDensityOfStates(const mat& energies, double delta, FILE* dosfile){

// 	double minE = energies.min();
// 	double maxE = energies.max();
// 	int nE = 2000; // Number of points in energy mesh

// 	vec energyMesh = arma::linspace(minE - 0.5, maxE + 0.5, nE);

// 	// Main loop
// 	double totalDos = 0;
// 	double dos = 0;
// 	// First loop over energies to normalize
// 	for(int n = 0; n < (int)energyMesh.n_elem; n++){
// 		double energy = energyMesh(n);
// 		double deltaE = energyMesh(1) - energyMesh(0);

// 		dos = densityOfStates(energy, delta, energies);
// 		totalDos += dos*deltaE;
// 	};
// 	// Main loop to write normalized DOS
// 	for(int n = 0; n < (int)energyMesh.n_elem; n++){
// 		double dos = 0;
// 		double energy = energyMesh(n);

// 		dos = densityOfStates(energy, delta, energies);
// 		dos = dos/totalDos; // Normallise
// 		fprintf(dosfile, "%lf\t%lf\n", energy, dos);
// 	};
// 	return;
// };
