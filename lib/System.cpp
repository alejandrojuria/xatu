#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "System.hpp"

using namespace arma;
using namespace std::chrono;

System::System(std::string filename){

	readConfigurationFile(filename);
	cout << "Correctly initiallized Zigzag object" << endl;

};

System::~System(){
	cout << "Destroying Zigzag object..." << endl;
};

/* Routine to extract information of the system in which we want to calculate
   excitons. Config. file format is tight-binder output style */
void System::readConfigurationFile(std::string filename){
    std::ifstream configfile (filename.c_str());
    std::string line;
    if (configfile.is_open()){

		getline(configfile, line);
		std::istringstream iss(line);
		iss >> ndim >> nmotif >> norbitals >> ncells;

		bravais_lattice(ndim, 3);
		motif(nmotif, 3);
		unitCellList(ncells, 3);
		hamiltonianMatrices(norbitals, norbitals, ncells);
		
		// Read bravais lattice vectors
        for(int i = 0; i < ndim; i++){
			double value;
			getline(configfile, line);
			std::istringstream iss(line);
			std::vector<double> array;
			
			while (iss >> value){
				array.push_back(value);
			}
			bravais_lattice.row(i) = arma::rowvec(array);
		};

		// Read motif
		for(int i = 0; i < nmotif; i++){
			double value;
			getline(configfile, line);
			std::istringstream iss(line);
			std::vector<double> array;
			
			while (iss >> value){
				array.push_back(value);
			}
			motif.row(i) = arma::rowvec(array);
		};

		double x, y, z;
		double re, im;
		char placeholder;
		std::string strValue;
		for(int i = 0; i < ncells; i++){
			getline(configfile, line);
			std::istringstream iss(line);
			iss >> x >> y >> z;
			unitCellList.row(i) = arma::rowvec{x, y, z};

			for(int j = 0; j < norbitals; j++){
				arma::cx_rowvec array(norbitals);
				getline(configfile, line);
				std::istringstream iss(line);
				
				int it = 0;
				while (iss >> strValue){
					std::istringstream iss2(strValue);
					// Read '(a + bi)' -> re=a, im=b
					iss2 >> placeholder >> re >> im >> placeholder;
					std::complex<double> value(re, im);
					array(it) = value;
					it++;
				};
				hamiltonianMatrices.slice(i).row(j) = array;
			}
			getline(configfile, line); // Skip # lines
		}

    }
	else cout << "Unable to open config. file";
};

int main(){
	System("model_bi");
}

/* Initialize Bloch hamiltonian for posterior diagonalization. Input: arma::vec k (wave number) */
cx_mat System::hamiltonian(arma::vec k) {

	cx_mat h = arma::zeros<cx_mat>(norbitals, norbitals);
	std::complex<double> i(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec a = unitCellList.row(i);
		h += hamiltonianMatrices.slice(i) * std::exp(-i*arma::dot(k, a));
	};

	return h;
};

/* Routine to generate the reciprocal lattice basis vectors from the bravais lattice basis. 
	The algorithm is based on the fact that
    a_i\dot b_j=2PI\delta_ij, which can be written as a linear system of
    equations to solve for b_j. Resulting vectors have 3 components independently
    of the dimension of the vector space they span.*/
mat System::calculate_reciprocal_lattice(){
	reciprocal_lattice = arma::zeros(ndim, 3);
	arma::mat coefficient_matrix = bravais_lattice;

	if (ndim == 1){
		reciprocal_lattice = 2*PI*coefficient_matrix / pow(arma::norm(coefficient_matrix), 2);
	}
	else{
		coefficient_matrix = coefficient_matrix.cols(0, ndim - 1);
		
		for (int i = 0; i < ndim; i++){
			arma::vec coefficient_vector = arma::zeros(ndim);
			coefficient_vector(i) = 2*PI;
			arma::rowvec reciprocal_vector;
			try{
				reciprocal_vector = arma::solve(coefficient_matrix, coefficient_vector);
			}
			catch (std::runtime_error e) {
				cout << "Failed to obtain reciprocal lattice vectors" << endl;
			};

			reciprocal_lattice.row(i) = reciprocal_vector;
		};
	};
};

/* Routine to compute a mesh of the first Brillouin zone using the
        Monkhorst-Pack algorithm. Returns a list of k vectors. 
		Takes two parameters:
		int n: Number of points to generate along one direction.
		int ny = 0. Defaults to 0. If given, generates*/
arma::mat System::brillouin_zone_mesh(int n){

	int nk = pow(n, ndim);
	arma::mat kpoints(nk, 3);
	arma::mat combinations = generate_combinations(n, ndim);

	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint(3);
		for (int j = 0; j < ndim; j++){
			kpoint += (2.*combinations.row(i)(j) - n)/(2*n)*reciprocal_lattice.row(j);
		}
		kpoints.row(i) = kpoint;
	}
	return kpoints;
}

arma::mat generate_combinations(int n, int ndim){
	int ncombinations = pow(n, ndim);
	arma::mat combinations(ncombinations, ndim);
	int it = 0;
	for (int i = 0; i < n; i++){
		for (int j = 0; j < n; j++){
			for (int k = 0; k < n; k++){
				combinations.row(it) = arma::rowvec{(float)i, (float)j, (float)k};
			}
		}
	}
	return combinations;
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


// /* Definition of non-interacting retarded Green function */
// std::complex<double> System::rGreenF(double energy, double delta, double eigEn){

// 	std::complex<double> i(0,1);
// 	return 1./((energy + i*delta) - eigEn);
// };

// /* Routine to calcule the density of states at a given energy,
// associated to a given set of eigenvalues (e.g. bulk or edge).
// NB: It is NOT normalized. */
// double System::densityOfStates(double energy, double delta, const mat& energies){

// 		double dos = 0;
// 		for(int i = 0; i < (int)energies.n_rows; i++){
// 			for(int j = 0; j < (int)energies.n_cols; j++){
// 				double eigEn = energies(i,j);
// 				dos += -PI*imag(rGreenF(energy, delta, eigEn));
// 			};
// 		};
// 		dos /= energies.n_cols*a; // Divide by number of k's and length a

// 		return dos;
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
