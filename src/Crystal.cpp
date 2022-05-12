#include "Crystal.hpp"
#include <numeric>


/* Initialize Crystal attributes from SystemConfiguration object */
void Crystal::initializeCrystalAttributes(const SystemConfiguration& configuration){
    ndim_           = configuration.systemInfo.ndim;
    bravaisLattice_ = configuration.systemInfo.bravaisLattice;
    motif_          = configuration.systemInfo.motif;
    atomToIndex     = configuration.systemInfo.atomToIndex;
	unitCellList_   = configuration.systemInfo.bravaisVectors;
    
    natoms_ = motif.n_rows;
	ncells_ = unitCellList.n_rows;
	nk_ = 0; // Mesh has to be explicitly initialized

    calculateReciprocalLattice();
    extractLatticeParameters();
}

/* Routine to generate the reciprocal lattice basis vectors from the bravais lattice basis. 
	The algorithm is based on the fact that
    a_i\dot b_j=2PI\delta_ij, which can be written as a linear system of
    equations to solve for b_j. Resulting vectors have 3 components independently
    of the dimension of the vector space they span.*/
void Crystal::calculateReciprocalLattice(){
	reciprocalLattice_ = arma::zeros(ndim, 3);
	arma::mat coefficient_matrix = bravaisLattice;

	if (ndim == 1){
		reciprocalLattice_ = 2.*PI*coefficient_matrix / pow(arma::norm(coefficient_matrix), 2);
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
			reciprocalLattice_.row(i) = reciprocal_vector.t();
		};
	};
};

/* Routine to compute a mesh of the first Brillouin zone using the
        Monkhorst-Pack algorithm. Returns a list of k vectors. 
		Always returns a mesh centered on the Gamma point, for both n even
		or odd. 
		int n: Number of points to generate along one direction. */
void Crystal::brillouinZoneMesh(int n){

	std::cout << "Creating BZ mesh... " << std::flush;

	int nk = pow(n, ndim);
	arma::mat kpoints(pow(n, ndim), 3);
	arma::mat combinations = generateCombinations(n, ndim);
	if (n % 2 == 1){
		combinations += 1./2;
	}
	
	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n)*reciprocalLattice_.row(j);
		}
		kpoints.row(i) = kpoint;
	}
	kpoints_ = kpoints;
	nk_ = kpoints.n_rows;
	std::cout << "Done. Number of k points: " << nk << std::endl;
}

arma::mat Crystal::brillouinZoneMeshCrystalCoordinates(int n){
	int nk = pow(n, ndim);
	arma::mat kpoints(pow(n, ndim), 3);
	arma::mat combinations = generateCombinations(n, ndim);
	if (n % 2 == 1){
		combinations += 1./2;
	}
	combinations = (2*combinations - n)/(2*n);

	return combinations;
}

/* Routine to obtain a kpoint mesh which is a subset of the full BZ mesh. 
If n is even, we substract one so that the mesh is symmetric under inversion */
void Crystal::reducedBrillouinZoneMesh(int n, int ncell){

	std::cout << "Creating BZ mesh... " << std::flush;

	n = (n % 2 == 0) ? n - 1 : n;
	int nk = pow(n, ndim);
	arma::mat kpoints(pow(n, ndim), 3);
	arma::mat combinations = generateCombinations(n, ndim);
	combinations += 1./2;
	
	int it = 0;
	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = arma::zeros<arma::rowvec>(3);
		// Remove "boundary" k so that mesh is symmery under k <-> -k
		if(!arma::all(combinations.row(i))){
			arma::cout << combinations.row(i) << arma::endl;
			continue;
		};
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*ncell)*reciprocalLattice_.row(j);
		}
		kpoints.row(it) = kpoint;
		it++;
	}
	kpoints_ = kpoints;
	nk_ = kpoints.n_rows;
	std::cout << "Done. Number of k points: " << nk << std::endl;
}


/* Routine to shift the center of the BZ mesh to a given point */
void Crystal::shiftBZ(const arma::rowvec& shift){
	arma::cout << "Shifting BZ mesh by vector: " << shift << arma::endl;
	if(shift.n_elem != 3){
		std::cout << "shift vector must be 3d" << std::endl;
		return; 
	}
	if(kpoints.empty()){
		std::cout << "To call this method kpoints must be initiallized first" << std::endl;
		return;
	}
	for(int i = 0; i < kpoints.n_rows; i++){
		kpoints_.row(i) += shift;
	}
}



/* Routine to restrict the kpoint mesh so that it has C3 symmetry. 
Note that this routine is intended to be used with systems with C3 symmetry */
void Crystal::preserveC3(){
	arma::vec norms = arma::diagvec((kpoints * kpoints.t()));
	std::vector<int> indices(norms.n_elem), coincidences, indicesToRemove;
	int coincidence;
	std::iota(indices.begin(), indices.end(), 1);
	
	for(const double norm : norms){
		coincidence = 0;
		if (norm == 0){ continue; } // Skip kpoint zero
		for(const int index : indices){
			if (norm == norms[index]){
				coincidence++;
				coincidences.push_back(index);
			}
		}
		for (const int index : coincidences){
			// Pop coincidence indices from indices
			// .erase() requires an iterator as parameter
			indices.erase(indices.begin() + index);
			if(coincidence < 3){
				indicesToRemove.push_back(index);		
			}
		}
	}
	// Generate kpoint matrix without the non-conserving C3 points
	std::vector<arma::uword> complementaryIndices;
	bool isRemoved;
	for(arma::uword i = 0; i < norms.n_elem; i++){
		isRemoved = false;
		for(const int index : indicesToRemove){
			isRemoved = (i == index) ? true : false;
		}
		if(!isRemoved){
			complementaryIndices.push_back(i);
		}
	}
	kpoints_ = kpoints_.rows(arma::uvec(complementaryIndices));
	nk_ = kpoints.n_rows;
}


/* Routine to generate a mesh for the BZ that preserves the C3 
symmetry of the hexagonal lattice */
void Crystal::brillouinZoneC3Mesh(int n){

	int nk = pow(n, ndim);
	nk = nk - (2*n - 1);
	int it = 0;
	arma::mat kpoints_block(nk, 3);
	arma::mat kpoints(3*nk + 3*n - 2, 3);
	double norm = arma::norm(reciprocalLattice_.row(0));
    arma::rowvec K = norm/sqrt(3)*(reciprocalLattice_.row(0)/2. -
                                	reciprocalLattice_.row(1)/2.)/arma::norm(
                                    reciprocalLattice_.row(0)/2. -
                                    reciprocalLattice_.row(1)/2.)/2;
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
    this->kpoints_ = kpoints;
    this->nk_ = kpoints.n_rows;
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
	this->a_ = arma::norm(bravaisLattice.row(0));

	double reference_height = motif.row(0)(2);
	double c = 0;
	for (arma::uword i = 0; i < motif.n_rows; i++){
		double diff = abs(motif.row(i)(2) - reference_height);
		if (diff > c){
			c = diff;
		}
	}
	if (c == 0){
		c = 1;
	}
	this->c_ = c;
}

arma::mat Crystal::wignerSeitzSupercell(int Ncell){

	// Generate combinations of [-1,0,1] to determine relevant lattice
	// vectors
	arma::mat lattice_combinations = generateCombinations(3, ndim, true);
	double norm = arma::norm(bravaisLattice.row(0));
	std::vector<arma::rowvec> lattice_vectors;

	for (arma::uword i = 0; i < lattice_combinations.n_rows; i++){
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
	for (arma::uword i = 0; i < midpoints.n_rows; i++){
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

arma::mat Crystal::generateCombinations(int nvalues, int ndim, bool centered){
	int ncombinations = pow(nvalues, ndim);
	arma::vec ones = arma::ones(nvalues);
	arma::mat combinations(ncombinations, ndim);
	arma::vec auxvector;
	arma::rowvec combination(ndim);
	int shift = centered ? (int)nvalues/2 : 0;
	for(int n = 0; n < ndim; n++){
		arma::vec values = arma::regspace(0, nvalues - 1) - shift;
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

arma::mat Crystal::supercellCutoff(int cutoff){
	arma::mat combinations = generateCombinations(cutoff, ndim, true);
	arma::mat cells = arma::zeros(combinations.n_rows, 3);

	for (int i = 0; i < combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * bravaisLattice.row(j);
		};
		cells.row(i) = lattice_vector;
	};

	return cells;
}

arma::mat Crystal::truncateSupercell(int ncell, double radius){

	arma::mat combinations = generateCombinations(ncell, ndim, true);
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

arma::mat Crystal::truncateReciprocalSupercell(int ncell, double radius){

	arma::mat combinations = generateCombinations(ncell, ndim, true);
	std::vector<arma::rowvec> cells_vector;
	for (int i = 0; i < combinations.n_rows; i++){
		arma::rowvec lattice_vector = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * reciprocalLattice.row(j);
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
							{sin(theta),  cos(theta), 0},
							{         0,		   0, 1}};
	
	arma::vec rotated_position = arma::inv(C3rotation)*(position.t());

	return rotated_position.t();
};

void Crystal::calculateInverseReciprocalMatrix(){
	arma::mat coefs = arma::zeros(ndim, ndim);
	coefs = reciprocalLattice * reciprocalLattice.t();
	arma::mat inverse;
	try{
		inverse = arma::inv(coefs);	
	}
	catch(std::runtime_error e){
		std::cout << "Unable to compute inverse reciprocal coefficients" << std::endl;
		throw;
	}
	this->inverseReciprocalMatrix = inverse;
}

// Given a kpoint, find its reciprocal coordinates to determine if it is outside
// of the BZ mesh. If it is, shift it by a reciprocal lattice vector so it falls
// on a point of the mesh. Intended to be used only with mesh of full BZ.
int Crystal::findEquivalentPointBZ(const arma::rowvec& kpoint, int ncell){
	if(inverseReciprocalMatrix.empty()){
		throw std::logic_error("Inverse reciprocal matrix has to be computed first.");
	}
	arma::vec independentTerm = reciprocalLattice * kpoint.t();
	arma::vec coefs = inverseReciprocalMatrix * independentTerm * 2*ncell;
	coefs = (ncell % 2 == 1) ? coefs - 1/2 : coefs; 
	coefs = arma::round(coefs);

	for(int i = 0; i < coefs.n_elem; i++){
		if (coefs(i) >= ncell){
			coefs(i) -= 2*ncell;
		}
		else if(coefs(i) < -ncell){
			coefs(i) += 2*ncell;
		}
		coefs(i) += ncell;
		coefs(i) /= 2;
	}
	int index = coefs(0) + coefs(1)*ncell;

	return index;
}