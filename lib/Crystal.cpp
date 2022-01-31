#include "Crystal.hpp"


/* Initialize Crystal attributes from SystemConfiguration object */
void Crystal::initializeCrystalAttributes(const SystemConfiguration& configuration){
    ndim           = configuration.systemInfo.ndim;
    bravaisLattice = configuration.systemInfo.bravaisLattice;
    motif          = configuration.systemInfo.motif;
    atomToIndex    = configuration.systemInfo.atomToIndex;
    
    natoms = motif.n_rows;

    calculateReciprocalLattice();
    extractLatticeParameters();
}

/* Routine to generate the reciprocal lattice basis vectors from the bravais lattice basis. 
	The algorithm is based on the fact that
    a_i\dot b_j=2PI\delta_ij, which can be written as a linear system of
    equations to solve for b_j. Resulting vectors have 3 components independently
    of the dimension of the vector space they span.*/
void Crystal::calculateReciprocalLattice(){
	reciprocalLattice = arma::zeros(ndim, 3);
	arma::mat coefficient_matrix = bravaisLattice;

	if (ndim == 1){
		reciprocalLattice = 2.*PI*coefficient_matrix / pow(arma::norm(coefficient_matrix), 2);
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
				std::cout << "Failed to obtain reciprocal lattice vectors" << std::endl;
				throw;
			};
			reciprocalLattice.row(i) = reciprocal_vector.t();
		};
	};
};

/* Routine to compute a mesh of the first Brillouin zone using the
        Monkhorst-Pack algorithm. Returns a list of k vectors. 
		Takes two parameters:
		int n: Number of points to generate along one direction.
		int ny = 0. Defaults to 0. If given, generates*/
arma::mat Crystal::brillouinZoneMesh(int n){

	int nk = pow(n, ndim);
	arma::mat kpoints(pow(n, ndim), 3);
	arma::mat combinations = generateCombinations(n, ndim);
	int it = 0;
	bool removeBoundary = false;

	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = arma::zeros<arma::rowvec>(3);
		if(removeBoundary){
			if(!arma::all(combinations.row(i))){
				continue;
			};
		}
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n)*reciprocalLattice.row(j);
		}
		kpoints.row(it) = kpoint;
		it++;
	}
	return kpoints;
}

/* Routine to generate a mesh for the BZ that preserves the C3 
symmetry of the hexagonal lattice */
arma::mat Crystal::c3BzMesh(int n){

	int nk = pow(n, ndim);
	nk = nk - (2*n - 1);
	int it = 0;
	arma::mat kpoints_block(nk, 3);
	arma::mat kpoints(3*nk + 3*n - 2, 3);
	double norm = arma::norm(reciprocalLattice.row(0));
    arma::rowvec K = norm/sqrt(3)*(reciprocalLattice.row(0)/2. -
                                	reciprocalLattice.row(1)/2.)/arma::norm(
                                    reciprocalLattice.row(0)/2. -
                                    reciprocalLattice.row(1)/2.)/2;
	arma::rowvec K_rotated = rotateC3(K);

	arma::mat combinations = generateCombinations(n, ndim);

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

void Crystal::extractLatticeParameters(){

	try{
		if (motif.is_empty() || bravaisLattice.is_empty()){
			throw "Error: Can not obtain lattice parameters (no Bravais lattice or motif)";
		}
	}
	catch (std::string e){
			std::cerr << e;
	}
	this->a = arma::norm(bravaisLattice.row(0));

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

arma::mat Crystal::wignerSeitzSupercell(int Ncell){

	// Generate combinations of [-1,0,1] to determine relevant lattice
	// vectors
	arma::mat lattice_combinations = generateCombinationsGamma(3, ndim);
	double norm = arma::norm(bravaisLattice.row(0));
	std::vector<arma::rowvec> lattice_vectors;

	for (int i = 0; i < lattice_combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += lattice_combinations.row(i)(j) * bravaisLattice.row(j);
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

	std::cout << midpoints << std::endl;
	std::cout << planes << std::endl;
	// Generate standard supercell
	arma::mat standard_supercell_coefs = generateCombinations(Ncell, ndim);
	arma::mat cells = arma::zeros(standard_supercell_coefs.n_rows, 3);
	arma::mat combinations = generateCombinations(3, ndim);

    for (int i = 0; i < standard_supercell_coefs.n_rows; i++){
        arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
        for (int j = 0; j < ndim; j++){
            lattice_vector += standard_supercell_coefs.row(i)(j) * bravaisLattice.row(j);
        }
        // Check is lattice vector plus some lattice vector is within the WS cell
        for(int n = 0; n < combinations.n_rows; n++){
            arma::rowvec translation = arma::zeros<arma::rowvec>(3);
            for (int m = 0; m < ndim; m++){
                translation += combinations.row(n)(m) * bravaisLattice.row(m) * Ncell;
            }
            std::cout << lattice_vector << std::endl;
            std::cout << translation << std::endl;
            arma::rowvec translated_vector = lattice_vector + translation;
            std::cout << translated_vector << std::endl;
            if (isInsideWsCell(translated_vector, planes, angles)){
                cells.row(i) = translated_vector;
                break;
            }
        }
    }	

	return cells;
}

bool Crystal::isInsideWsCell(const arma::rowvec& point, 
							   const arma::mat& planes, const arma::rowvec& angles){

	arma::rowvec checks(angles.n_elem);
	for (int i = 0; i < angles.n_elem; i++){
		double angle = angles(i);
		arma::rowvec plane = planes.row(i);
		double side = plane(0)*point(0) + plane(1)*point(1) + plane(2);
		std::cout << angle << std::endl;
		std::cout << side << std::endl;
		std::cout << "---------" << std::endl;
		if (side <= 0){
			checks(i) = 1;
		}
		else{
			checks(i) = 0;
		}
	}
	std::cout<< checks << std::endl;
	bool is_inside = false;
	if (arma::all(checks)){
		is_inside = true;
	}

	return is_inside;
};

arma::mat Crystal::generateCombinations(int nvalues, int ndim){
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

arma::mat Crystal::generateCombinationsGamma(int nvalues, int ndim){
	if (nvalues%2 == 0){
		nvalues++;
	}
	int ncombinations = pow(nvalues, ndim);
	arma::vec ones = arma::ones(nvalues);
	arma::mat combinations(ncombinations, ndim);
	arma::vec auxvector;
	arma::rowvec combination(ndim);
	std::cout << __LINE__ << std::endl;
	for(int n = 0; n < ndim; n++){
		arma::vec values = arma::regspace(-(int)nvalues/2, (int)nvalues/2);
		for(int i = 0; i < ndim - n - 1; i++){
			values = arma::kron(ones, values);
		}
		std::cout << __LINE__ << std::endl;
		for(int j = 0; j < n; j++){
			values = arma::kron(values, ones);
		}
		std::cout << __LINE__ << std::endl;
		combinations.col(n) = values;
		std::cout << __LINE__ << std::endl;
	}

	std::cout << __LINE__ << std::endl;
	return combinations;
}

arma::mat Crystal::truncateSupercell(int Ncell, double radius){

	arma::mat combinations = generateCombinations(Ncell, ndim);
	std::vector<arma::rowvec> cells_vector;
	for (int i = 0; i < combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * bravaisLattice.row(j);
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
arma::rowvec Crystal::rotateC3(const arma::rowvec& position){
	double theta = 2*PI/3;
	arma::mat C3rotation = {{cos(theta), -sin(theta), 0},
							{sin(theta), cos(theta) , 0},
							{         0,		   0, 1}};
	arma::vec rotated_position = C3rotation*position.t();

	return rotated_position.t();
}