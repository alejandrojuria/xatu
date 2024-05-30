#include "xatu/Lattice.hpp"

namespace xatu {

/**
 * Copy constructor.
 * 
 * Creates a Lattice object using another Lattice object.
 * @param lattice Lattice to be copied. 
 */
Lattice::Lattice(const Lattice& lattice){

	ndim_ 			= lattice.ndim;
	bravaisLattice_ = lattice.bravaisLattice;
	motif_          = lattice.motif;
	unitCellList_   = lattice.unitCellList;
    
    natoms_ = motif.n_rows;
	ncells_ = unitCellList.n_rows;

	calculateReciprocalLattice();
	extractLatticeParameters();
}

/**
 * Method to initialize Lattice attributes from SystemConfiguration object.
 * 
 * @param configuration SystemConfiguration object. 
 */
void Lattice::initializeLatticeAttributes(const SystemConfiguration& configuration){

    ndim_           = configuration.systemInfo.ndim;	
    bravaisLattice_ = configuration.systemInfo.bravaisLattice;
    motif_          = configuration.systemInfo.motif;
	unitCellList_   = configuration.systemInfo.bravaisVectors;
    
    natoms_ = motif.n_rows;
	ncells_ = unitCellList.n_rows;

	nk_ = 0; // Mesh has to be explicitly initialized

    calculateReciprocalLattice();
    extractLatticeParameters();
	computeUnitCellVolume();
}

/**
 *  Routine to generate the reciprocal lattice basis vectors from the bravais lattice basis. 
 *	@details The algorithm is based on the fact that
 *  a_i\dot b_j=2PI\delta_ij, which can be written as a linear system of
 *  equations to solve for b_j. Resulting vectors have 3 components independently
 *  of the dimension of the vector space they span.
 */
void Lattice::calculateReciprocalLattice(){
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

/**
 * Method to obtain the lattice parameters of 2D lattices.
 * @details For simplicity, it takes a as the norm of the first Bravais vector
 * and c as the height of the 2D crystal, taking as reference the hexagonal lattice.
*/
void Lattice::extractLatticeParameters(){

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

/**
 * Method to compute the unit cell volume, area or length of the crystal.
 * @details Method adapted to compute the relevant quantity depending on the dimensionality
 * of the problem (1D = length, 2D = area, 3D = volume).
 * @return void
 */
void Lattice::computeUnitCellVolume(){
	if(ndim == 1){
		this->unitCellVolume_ = arma::norm(bravaisLattice.row(0));
	}
	else if(ndim == 2){
		arma::rowvec crossProduct = arma::cross(bravaisLattice.row(0), bravaisLattice.row(1));
		this->unitCellVolume_ = arma::norm(crossProduct);
	}
	else if(ndim == 3){
		this->unitCellVolume_ = std::abs(arma::det(bravaisLattice));
	}
}

/* --------------------------- Brillouin zone methods --------------------------- */

/**
 * Routine to compute a mesh of the first Brillouin zone using the Monkhorst-Pack algorithm.
 * @details Always returns a mesh centered at the Gamma point, for both n even or odd. 
 * @param n Number of points along one of the axis. 
 */
void Lattice::brillouinZoneMesh(int n){

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
	meshBZ_ = kpoints;
	nk_ = kpoints.n_rows;
	std::cout << "Done. Number of k points in BZ mesh: " << nk << std::endl;
}

/**
 * Method to generate a kpoint mesh which is a subset of the full BZ mesh. 
 * @details If n is even, we substract one so that the mesh is symmetric under inversion.
 * The reduction factor specifies is used to create a mesh of the full BZ of factor*n points,
 * from which we extract the corresponding mesh for n points centered at Gamma.
 * @param n Number of points along one axis.
 * @param factor Reduction factor of the mesh.
 * */
void Lattice::reducedBrillouinZoneMesh(int n, int factor){

	// First create mesh of whole BZ
	brillouinZoneMesh(n*factor);

	// Now create submesh
	int nk = pow(n, ndim);

	arma::mat kpoints(nk, 3);
	arma::mat combinations = generateCombinations(n, ndim);
	if (n % 2 == 1){
		combinations += 1./2;
	}
	
	for (int i = 0; i < nk; i++){
		arma::rowvec kpoint = arma::zeros<arma::rowvec>(3);
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n*factor)*reciprocalLattice_.row(j);
		}
		kpoints.row(i) = kpoint;
	}
	kpoints_ = kpoints;
	nk_ = nk;
	factor_ = factor;
	std::cout << "Number of k points in submesh: " << nk << std::endl;
}

/**
 * Method to shift the center of the BZ mesh to a given point.
 * @param shift Vector to shift center of BZ mesh. 
 */
void Lattice::shiftBZ(const arma::rowvec& shift){
	std::cout << std::left << std::setw(30) << "Shifting BZ mesh by vector: ";
	for (auto si : shift){
		std::cout << si << "  ";
	}
	std::cout << std::endl;
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

/**
 * Auxiliary method to compute the inverse reciprocal matrix.
 * @details The reciprocal matrix is defined as the R*R.t, where R is the matrix
 * containing the reciprocal lattice vectors as row. This inverted matrix is required
 * to map any kpoint back to the original Monkhorst-Pack mesh.
*/
void Lattice::calculateInverseReciprocalMatrix(){
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

/**
 * Method to determine a kpoint equivalent to another within the BZ mesh.
 * @details Given a kpoint, this function finds its reciprocal coordinates to determine if 
 * it is outside of the BZ mesh. If it is, shift it by a reciprocal lattice vector so it falls
 * on a point of the mesh. Intended to be used only with mesh of full BZ (i.e. Monkhorst-Pack).
 * @param kpoint kpoint to be mapped back.
 * @param ncell Number of points used in the original BZ mesh, usually equivalent to the number of cells.
 * @returns Index (row) of the equivalent kpoint from the BZ mesh matrix.
 */ 
int Lattice::findEquivalentPointBZ(const arma::rowvec& kpoint, int ncell){
	if(inverseReciprocalMatrix.empty()){
		calculateInverseReciprocalMatrix();
	}
	ncell = ncell * factor_;
	arma::vec independentTerm = reciprocalLattice * kpoint.t();
	arma::vec coefs = inverseReciprocalMatrix * independentTerm * 2*ncell;
	coefs = (ncell % 2 == 1) ? coefs - 1 : coefs; 
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
	int index = 0;
	
	std::vector<int> cells_array = {1, ncell, ncell*ncell}; // Auxiliar array to avoid using std::pow
	for(int i = 0; i < coefs.n_elem; i++){
		index += coefs(i)*cells_array[i];
	}

	return index;
};

/* --------------------------- Supercell methods --------------------------- */

/**
 * Method to generate a list with combinations of Bravais vectors.
 * @details Each row of the list stores the integers which correspond
 * to the coefficients of the linear combinations of Bravais basis vectors.
 * @param nvalues Number of cells in one direction.
 * @param ndim Dimension of the system.
 * @param centered If true, the combinations are centered at the origin cell.
 * 				   If false, all combinations have positive coefficients.
 * @returns List of cell combinations.
*/
arma::mat Lattice::generateCombinations(int nvalues, int ndim, bool centered){
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

/**
 * Method to generate a list of cells within a sphere of specified radius.
 * @param ncell List of cells along one axis.
 * @param radius Cutoff radius of the sphere.
 * @returns List of cells within the sphere.
*/
arma::mat Lattice::truncateSupercell(int ncell, double radius){

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

/**
 * Method to generate a list of reciprocal cells contained within a sphere of specified radius.
 * @param ncell List of cells in one direction.
 * @param radius Radius of the cutoff sphere.
 * @returns List of reciprocal cells in cartesian coordinates.
*/
arma::mat Lattice::truncateReciprocalSupercell(int ncell, double radius){

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

/* --------------------------- Other methods --------------------------- */

/**
 * Routine to rotate a position by 2pi/3, either on real space
 * or on reciprocal space to enforce C3 symmetry.
 * @param position Vector to rotate.
 * @returns Rotated vector.
 */
arma::rowvec Lattice::rotateC3(const arma::rowvec& position){
	double theta = 2*PI/3;
	arma::mat C3rotation = {{cos(theta), -sin(theta), 0},
							{sin(theta),  cos(theta), 0},
							{         0,		   0, 1}};
	
	arma::vec rotated_position = arma::inv(C3rotation)*(position.t());

	return rotated_position.t();
};


}

