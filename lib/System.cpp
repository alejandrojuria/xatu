#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include "System.hpp"

using namespace arma;
using namespace std::chrono;

System::System(std::string filename, std::string overlapFile, bool isReal){

	readConfigurationFile(filename, isReal);
	if (!overlapFile.empty()){
		readOverlapFile(overlapFile, isReal);
	}
	extractLatticeParameters();
	calculate_reciprocal_lattice();
	cout << "Correctly initiallized System object" << endl;

};

System::~System(){
	cout << "Destroying System object..." << endl;
};

/* Routine to extract information of the system in which we want to calculate
   excitons. Config. file format is tight-binder output style */
void System::readConfigurationFile(std::string filename, bool isReal){
    std::ifstream configfile (filename.c_str());
    std::string line;
    if (configfile.is_open()){
		getline(configfile, line);
		std::istringstream iss(line);
		try{
			iss >> ndim >> natoms >> norbitals >> ncells;
			if (ndim > 3){
				throw "Error: ndim must be lower than three";
			};
		}
		catch (std::string e){
			std::cerr << e;
		};
		// basisdim = natoms * norbitals;
		basisdim = norbitals; // Monkeypatch for WSe2
		bravais_lattice = arma::zeros(ndim, 3);
		motif = arma::zeros(natoms, 3);
		unitCellList = arma::zeros(ncells, 3);
		hamiltonianMatrices = arma::cx_cube(basisdim, basisdim, ncells);
		
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
		for(int i = 0; i < natoms; i++){
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

			for(int j = 0; j < basisdim; j++){
				arma::cx_rowvec array = arma::zeros<cx_rowvec>(basisdim);
				getline(configfile, line);
				std::istringstream iss(line);

				int it = 0;
				while (iss >> strValue){
					std::istringstream iss2(strValue);
					if (isReal){
						iss2 >> re;
						im = 0.0;
					}
					else{
						// Read '(a + bi)' -> re=a, im=b
						iss2 >> placeholder >> re >> im >> placeholder;
					}
					if (it == j){
						re /= 2.;
					}
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

/* In case that the basis used is not orthogonal, we need to the overlaps 
   as well to solve the generalized eigenvalue problem */
void System::readOverlapFile(std::string filename, bool isReal){
	overlapMatrices = arma::cx_cube(basisdim, basisdim, ncells);
	std::ifstream configfile (filename.c_str());
    std::string line;
    if (configfile.is_open()){
		double x, y, z;
		double re, im;
		char placeholder;
		std::string strValue;
		for(int i = 0; i < ncells; i++){
			getline(configfile, line);
			std::istringstream iss(line);
			iss >> x >> y >> z;
			unitCellList.row(i) = arma::rowvec{x, y, z};

			for(int j = 0; j < basisdim; j++){
				arma::cx_rowvec array = arma::zeros<cx_rowvec>(basisdim);
				getline(configfile, line);
				std::istringstream iss(line);

				int it = 0;
				while (iss >> strValue){
					std::istringstream iss2(strValue);
					if (isReal){
						iss2 >> re;
						im = 0.0;
					}
					else{
						// Read '(a + bi)' -> re=a, im=b
						iss2 >> placeholder >> re >> im >> placeholder;
					}
					if (it == j){
						re /= 2.;
					}
					std::complex<double> value(re, im);
						
					array(it) = value;
					it++;
				};
				overlapMatrices.slice(i).row(j) = array;
			}
			getline(configfile, line); // Skip # lines
		}
		
    }
	else cout << "Unable to open config. file";
}

/* Bloch hamiltonian for posterior diagonalization. Input: arma::vec k (wave number) */
arma::cx_mat System::hamiltonian(arma::rowvec k, bool isTriangular){

	cx_mat h = arma::zeros<cx_mat>(basisdim, basisdim);
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

	cx_mat s = arma::zeros<cx_mat>(basisdim, basisdim);
	std::complex<double> imag(0, 1);
	for (int i = 0; i < ncells; i++){
		arma::rowvec cell = unitCellList.row(i);
		s += overlapMatrices.slice(i) * std::exp(-imag*arma::dot(k, cell));
	};

	s += s.t();
	return s;
}

/* Routine to generate the reciprocal lattice basis vectors from the bravais lattice basis. 
	The algorithm is based on the fact that
    a_i\dot b_j=2PI\delta_ij, which can be written as a linear system of
    equations to solve for b_j. Resulting vectors have 3 components independently
    of the dimension of the vector space they span.*/
void System::calculate_reciprocal_lattice(){
	reciprocal_lattice = arma::zeros(ndim, 3);
	arma::mat coefficient_matrix = bravais_lattice;

	if (ndim == 1){
		reciprocal_lattice = 2.*PI*coefficient_matrix / pow(arma::norm(coefficient_matrix), 2);
	}
	else{
		coefficient_matrix = coefficient_matrix.cols(0, ndim - 1);
		
		for (int i = 0; i < ndim; i++){
			arma::vec coefficient_vector = arma::zeros(ndim);
			coefficient_vector(i) = 2*PI;
			arma::vec reciprocal_vector = arma::zeros(3, 1);
			try{
				reciprocal_vector.rows(0, ndim - 1) = arma::solve(coefficient_matrix, coefficient_vector);
			}
			catch (std::runtime_error e) {
				cout << "Failed to obtain reciprocal lattice vectors" << endl;
			};
			reciprocal_lattice.row(i) = reciprocal_vector.t();
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
	arma::mat kpoints(pow(n, ndim), 3);
	arma::mat combinations = generate_combinations(n, ndim);
	int it = 0;
	bool removeBoundary = false;

	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = arma::zeros<rowvec>(3);
		if(removeBoundary){
			if(!arma::all(combinations.row(i))){
				continue;
			};
		}
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n)*reciprocal_lattice.row(j);
		}
		kpoints.row(it) = kpoint;
		it++;
	}
	return kpoints;
}

/* Routine to generate a mesh for the BZ that preserves the C3 
symmetry of the hexagonal lattice */
arma::mat System::C3_BZ_Mesh(int n){

	int nk = pow(n, ndim);
	nk = nk - (2*n - 1);
	int it = 0;
	arma::mat kpoints_block(nk, 3);
	arma::mat kpoints(3*nk + 3*n - 2, 3);
	double norm = arma::norm(reciprocal_lattice.row(0));
    arma::rowvec K = norm/sqrt(3)*(reciprocal_lattice.row(0)/2. -
                                	reciprocal_lattice.row(1)/2.)/arma::norm(
                                    reciprocal_lattice.row(0)/2. -
                                    reciprocal_lattice.row(1)/2.)/2;
	arma::rowvec K_rotated = rotateC3(K);

	arma::mat combinations = generate_combinations(n, ndim);

	for (int i = 0; i < combinations.n_rows; i++){
		arma::rowvec kpoint = arma::zeros<arma::rowvec>(3);
		if(combinations.row(i)(0) == 0 || combinations.row(i)(1) == 0){
			continue;
		}
		kpoint =  combinations.row(i)(0)/(n-1)*K + 
				  combinations.row(i)(1)/(n-1)*K_rotated;
		kpoints_block.row(it) = kpoint;
		it++;
	}
	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = kpoints_block.row(i);
		arma::rowvec kpoint_rotated = rotateC3(kpoint);
		arma::rowvec kpoint_rotated_twice = rotateC3(kpoint_rotated);

		kpoints.row(i) = kpoint;
		kpoints.row(nk + i) = kpoint_rotated;
		kpoints.row(2*nk + i) = kpoint_rotated_twice;
	}
	for (int i = 1; i < n; i++){
		arma::rowvec kpoint = (double)i/(n-1)*K;
		arma::rowvec kpoint_rotated = rotateC3(kpoint);
		arma::rowvec kpoint_rotated_twice = rotateC3(kpoint_rotated);

		kpoints.row(3*nk + i - 1) = kpoint;
		kpoints.row(3*nk + n + i - 2) = kpoint_rotated;
		kpoints.row(3*nk + 2*n + i - 3) = kpoint_rotated_twice;
	}
	kpoints.row(kpoints.n_rows - 1) = arma::rowvec{0, 0, 0};

	return kpoints;
}

void System::extractLatticeParameters(){

	try{
		if (motif.is_empty() || bravais_lattice.is_empty()){
			throw "Error: Can not obtain lattice parameters (no Bravais lattice or motif)";
		}
	}
	catch (std::string e){
			std::cerr << e;
	}
	this->a = arma::norm(bravais_lattice.row(0));

	double reference_height = motif.row(0)(2);
	double c = 0;
	for (int i = 0; i < motif.n_rows; i++){
		double diff = abs(motif.row(i)(2) - reference_height);
		if (diff > c){
			c = diff;
		}
	}
	if (c == 0){
		c = 1;
	}
	this->c = c;
}

arma::mat System::wigner_seitz_supercell(int Ncell){

	// Generate combinations of [-1,0,1] to determine relevant lattice
	// vectors
	arma::mat lattice_combinations = generate_combinations_gamma(3, ndim);
	double norm = arma::norm(bravais_lattice.row(0));
	std::vector<arma::rowvec> lattice_vectors;

	for (int i = 0; i < lattice_combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += lattice_combinations.row(i)(j) * bravais_lattice.row(j);
		}
		if (abs(arma::norm(lattice_vector) - norm) < 0.001){
			lattice_vectors.push_back(lattice_vector);
		}
	}
	// Obtain midpoints of said supercell lattice vectors
	arma::mat midpoints(lattice_vectors.size(), 3);
	int it = 0;
	for (const auto& vector : lattice_vectors){
		arma::rowvec large_lattice_vector = vector * Ncell;
		midpoints.row(it) = large_lattice_vector/2;
		it++;
	}
	// Calculate angle of each point with respect to origin
	arma::rowvec angles(midpoints.n_rows);
	angles = arma::atan2(midpoints.col(1), midpoints.col(0)).t();

	// Determine perpendicular planes to each midpoint
	arma::mat planes(midpoints.n_rows, 3);
	for (int i = 0; i < midpoints.n_rows; i++){
		double A = midpoints.row(i)(0);
		double B = midpoints.row(i)(1);
		double d = -A*A - B*B;
		planes.row(i) = arma::rowvec{A, B, d};
	}

	cout << midpoints << endl;
	cout << planes << endl;
	// Generate standard supercell
	arma::mat standard_supercell_coefs = generate_combinations(Ncell, ndim);
	arma::mat cells = arma::zeros(standard_supercell_coefs.n_rows, 3);
	arma::mat combinations = generate_combinations_gamma(3, ndim);
	for (int i = 0; i < standard_supercell_coefs.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += standard_supercell_coefs.row(i)(j) * bravais_lattice.row(j);
		}
		// Check is lattice vector plus some lattice vector is within the WS cell
		for(int n = 0; n < combinations.n_rows; n++){
			arma::rowvec translation = arma::zeros<arma::rowvec>(3);
			for (int m = 0; m < ndim; m++){
				translation += combinations.row(n)(m) * bravais_lattice.row(m) * Ncell;
			}
			cout << lattice_vector << endl;
			cout << translation << endl;
			arma::rowvec translated_vector = lattice_vector + translation;
			cout << translated_vector << endl;
			if (is_inside_ws_cell(translated_vector, planes, angles)){
				cells.row(i) = translated_vector;
				break;
			}
		}
	}	

	return cells;
}

bool System::is_inside_ws_cell(const arma::rowvec& point, 
							   const arma::mat& planes, const arma::rowvec& angles){

	arma::rowvec checks(angles.n_elem);
	for (int i = 0; i < angles.n_elem; i++){
		double angle = angles(i);
		arma::rowvec plane = planes.row(i);
		double side = plane(0)*point(0) + plane(1)*point(1) + plane(2);
		cout << angle << endl;
		cout << side << endl;
		cout << "---------" << endl;
		if (side <= 0){
			checks(i) = 1;
		}
		else{
			checks(i) = 0;
		}
	}
	cout<< checks << endl;
	bool is_inside = false;
	if (arma::all(checks)){
		is_inside = true;
	}

	return is_inside;
};

arma::mat System::generate_combinations(int nvalues, int ndim){
	int ncombinations = pow(nvalues, ndim);
	arma::vec ones = arma::ones(nvalues);
	arma::mat combinations(ncombinations, ndim);
	arma::vec auxvector;
	arma::rowvec combination(ndim);
	for(int n = 0; n < ndim; n++){
		arma::vec values = arma::regspace(0, nvalues - 1);
		for(int i = 0; i < ndim - n - 1; i++){
			values = arma::kron(ones, values);
		}
		for(int j = 0; j < n; j++){
			values = arma::kron(values, ones);
		}
		combinations.col(n) = values;
	}

	return combinations;
}

arma::mat System::generate_combinations_gamma(int nvalues, int ndim){
	if (nvalues%2 == 0){
		nvalues++;
	}
	int ncombinations = pow(nvalues, ndim);
	arma::vec ones = arma::ones(nvalues);
	arma::mat combinations(ncombinations, ndim);
	arma::vec auxvector;
	arma::rowvec combination(ndim);
	for(int n = 0; n < ndim; n++){
		arma::vec values = arma::regspace(-(int)nvalues/2, (int)nvalues/2);
		for(int i = 0; i < ndim - n - 1; i++){
			values = arma::kron(ones, values);
		}
		for(int j = 0; j < n; j++){
			values = arma::kron(values, ones);
		}
		combinations.col(n) = values;
	}

	return combinations;
}

arma::mat System::truncate_supercell(int Ncell, double radius){

	arma::mat combinations = generate_combinations_gamma(Ncell, ndim);
	std::vector<arma::rowvec> cells_vector;
	for (int i = 0; i < combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * bravais_lattice.row(j);
		};
		double distance = arma::norm(lattice_vector);
		if (distance < radius + 1E-5){
			cells_vector.push_back(lattice_vector);
		}
	}
	int total_cells = cells_vector.size();
	arma::mat cells = arma::zeros(total_cells, 3);
	for (int i = 0; i < total_cells; i++){
		cells.row(i) = cells_vector[i];
	}

	return cells;
}

/* Routine to rotate a position by 2pi/3, either on real space
or on reciprocal space to enforce C3 symmetry */
arma::rowvec System::rotateC3(const arma::rowvec& position){
	double theta = 2*PI/3;
	arma::mat C3rotation = {{cos(theta), -sin(theta), 0},
							{sin(theta), cos(theta) , 0},
							{         0,		   0, 1}};
	arma::vec rotated_position = C3rotation*position.t();

	return rotated_position.t();
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
