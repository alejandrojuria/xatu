#include "xatu/SystemTB.hpp"

namespace xatu {

/**
 * Configuration constructor.
 * @details Constructor which takes in a ConfigurationSystem object, i.e.
 * to init a SystemTB from a configuration file.
 * @param SystemConfig ConfigurationSystem object obtained from config. file.
 */
SystemTB::SystemTB(const ConfigurationSystem& SystemConfig) : System(SystemConfig){

	this->isTriangular_ 	    = checkIfTriangular(SystemConfig.hamiltonianMatrices.slice(0));
	this->motif_                = SystemConfig.motif;
	this->natoms_               = motif.n_cols;
	this->orbitalsPerSpecies_   = SystemConfig.orbitalsPerSpecies;

	extractLatticeParameters();

}

/**
 * Method to set the type of system (CRYSTAL or not).
 * @details The isCRYSTAL attribute ensures proper handling of the energies (to be converted from atomic units to eV).
 * @param isCRYSTAL Whether the system is CRYSTAL or not.
 * @return void.
 */
void SystemTB::setCRYSTAL(const bool isCRYSTAL){

	this->isCRYSTAL_ = isCRYSTAL;

}

/**
 * Filling setter.
 * @param filling Number of filled bands.
 * @return void.
*/
void SystemTB::setFilling(const int filling){

	if (filling > 0){
		highestValenceBand_ = filling - 1;
	}
	else{
		std::cout << "Filling must be a positive integer" << std::endl;
	}

}

// /*------------------------------- Bloch Hamiltonian -------------------------------*/

/**
 * Bloch Hamiltonian.
 * @details The Bloch Hamiltonian is constructed at each k point using the Fock
 * matrices of the system, H(R). The Fock matrices can be either complete or triangular.
 * @param k k-point where H(k) is evaluated.
 * @return cx_mat Bloch Hamiltonian matrix at k.
 */
arma::cx_mat SystemTB::hamiltonian(const arma::colvec& k) const{

	arma::cx_mat h = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	for (int i = 0; i < ncells; i++){
		arma::colvec cell = Rlist.col(i);
		h += (*ptr_hamiltonianMatrices).slice(i) * std::exp(imag*arma::dot(k, cell));
	}

	if (isTriangular){
		h.diag() -= h.diag()/2;
		h += h.t();
	}

	return h;

}

/**
 * Overlap matrix in reciprocal space.
 * @details The reciprocal overlap matrix S(k) is required to solve the generalized
 * eigenvalue problem, which appears with non-orthonormal basis. It is built from the
 * overlap matrices in real space S(R).
 * @param k k-point where S(k) is evaluated.
 * @return cx_mat Reciprocal overlap matrix S(k).
 */
arma::cx_mat SystemTB::overlap(const arma::colvec& k) const{

	arma::cx_mat s = arma::zeros<arma::cx_mat>(norbitals, norbitals);
	for (int i = 0; i < ncells; i++){
		arma::colvec cell = Rlist.col(i);
		s += (*ptr_overlapMatrices).slice(i) * std::exp(imag*arma::dot(k, cell));
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
 * or not (if they are triangular it is also assumed that the basis is non-orthonormal),
 * the returning eigenvectors are orthonormalized.
 * @param k k-point where the eigenvalues and eigenvectors are computed.
 * @param eigval Vector to store the energies of the system.
 * @param eigvec Complex matrix to store eigenvectors.
 * @return void.
*/
void SystemTB::solveBands(const arma::colvec& k, arma::vec& eigval, arma::cx_mat& eigvec) const {

	arma::cx_mat h = hamiltonian(k);
	if (!(*ptr_overlapMatrices).empty()){ 
		if (isCRYSTAL){
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
 * @return void.
*/
void SystemTB::orthogonalize_hamiltonian(const arma::colvec& k, arma::cx_mat& hamiltonian) const {

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
 * @param k k-point required to perform the transformation.
 * @return arma::cx_vec Atomic gauge state coefficients.
 */
arma::cx_vec SystemTB::latticeToAtomicGauge(const arma::cx_vec& coefs, const arma::colvec& k){

    arma::cx_vec phases(norbitals);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.col(atomIndex)(3);
        for(uint orbIndex = 0; orbIndex < orbitalsPerSpecies(species); orbIndex++){
            arma::colvec atomPosition = motif.col(atomIndex).subvec(0, 2);
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
 * @param k k-point used in the transformation.
 * @return arma::cx_vec Lattice gauge coefficients.
 */
arma::cx_vec SystemTB::atomicToLatticeGauge(const arma::cx_vec& coefs, const arma::colvec& k){

    arma::cx_vec phases(norbitals);
    int it = 0;
    for(int atomIndex = 0; atomIndex < natoms; atomIndex++){
        int species = motif.col(atomIndex)(3);
        for(uint orbIndex = 0; orbIndex < orbitalsPerSpecies(species); orbIndex++){
            arma::colvec atomPosition = motif.col(atomIndex).subvec(0, 2);
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
 * alternating spin consecutively, or kron(orbital,spin).
 * @param eigvec State.
 * @returns double Expectation value of Sz.
 * */
double SystemTB::expectedSpinZValue(const arma::cx_vec& eigvec){

	arma::cx_vec spinEigvalues = {1./2, -1./2};
    arma::cx_vec spinVector = arma::kron(arma::ones(norbitals/2), spinEigvalues);
	arma::cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));

}

/**
 * Routine to calculate the expected value of the spin component Sy.
 * TODO: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sy.
 */
// double SystemTB::expectedSpinYValue(const arma::cx_vec& eigvec){

// 	std::complex<double> i(0,1);
// 	arma::cx_mat operatorSy = arma::zeros<arma::cx_mat>(norbitals, norbitals);
// 	operatorSy(0,1) = -i;
// 	operatorSy(1,0) = i;
// 	operatorSy.submat(2,5, 4,7) = -i*arma::eye<arma::cx_mat>(3,3);
// 	operatorSy.submat(5,2, 7,4) = i*arma::eye<arma::cx_mat>(3,3);
	
// }

/**
 * Routine to calculate the expected value of the spin component Sx.
 * TODO: Not correctly implemented yet (incorrect basis ordering).
 * @param eigvec State.
 * @returns Expectation value of Sx.
 */
// double SystemTB::expectedSpinXValue(const arma::cx_vec& eigvec){

// 	arma::cx_mat operatorSx = arma::zeros<arma::cx_mat>(norbitals, norbitals);
// 	operatorSx(0,1) = 1;
// 	operatorSx(1,0) = 1;
// 	operatorSx.submat(2,5, 4,7) = arma::eye<arma::cx_mat>(3,3);
// 	operatorSx.submat(5,2, 7,4) = arma::eye<arma::cx_mat>(3,3);

// }

/**
 * Routine to calculate the velocity matrix element between two bands.
 * @details This routine is intended to work with tight-binding systems only, and will produce
 * incorrect results otherwise.
 * TODO: Duplicity with velocitySingleParticle in the Result class, merge both routines.
 * @param k kpoint where we want to evaluate the velocity.
 * @param fBand Index of the first (row) band.
 * @param sBand Index of the second (column) band.
 * @return cx_vec Vector with the 3 spatial components of the velocity matrix elements between the selected bands.
*/
arma::cx_vec SystemTB::velocity(const arma::colvec& k, int fBand, int sBand) const {
	arma::cx_mat h = hamiltonian(k);
	arma::vec eigval;
	arma::cx_mat eigvec;
	arma::eig_sym(eigval, eigvec, h);

	arma::cx_vec fBandEigvec = eigvec.col(fBand);
	arma::cx_vec sBandEigvec = eigvec.col(sBand);

    arma::cx_cube hkDerivative = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);
    arma::cx_cube iHt = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);

    // First compute Hk derivative
    for (int j = 0; j < 3; j++){
        for (int i = 0; i < ncells; i++){
            arma::colvec cell = Rlist.col(i);
            hkDerivative.slice(j) += (*ptr_hamiltonianMatrices).slice(i) * 
                                     std::exp(imag*arma::dot(k, cell)) * cell(j) * imag;
	    }
    }

    // Next compute iH(t-t') matrix
    arma::cx_cube motifDifference = arma::zeros<arma::cx_cube>(norbitals, norbitals, 3);
    arma::cx_mat extendedMotif = arma::zeros<arma::cx_mat>(norbitals, 3);
    int currentIndex = 0;
    for (int i = 0; i < natoms; i++){
        int norb = orbitalsPerSpecies(motif.col(i)(3));
        extendedMotif.rows(currentIndex, currentIndex + norb - 1) = arma::kron((motif.col(i).subvec(0, 2)).t(),
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
}

// /*------------------------------------- Modifiers -------------------------------------*/

/**
 * Method to add a Zeeman (on-site) term to the Hamiltonian. To be used in the TB mode. 
 * @details This method assumes that the Hamiltonian incorporates spin in the following way: 
 * |i1,up>,|i1,down>,...,|in,up>,|in,down>, where i runs over orbitals. This ordering is exclusive to TB mode.
 * IMPORTANT: since no non-const copies of the Hamiltonian are now made in System, the attribute hamiltonianMatrices 
 * can no longer be modified to include the Zeeman term. Instead, a new attribute for the Zeeman-modified Hamiltonian 
 * has been created, and it is initialized only by this method. Care must be taken in using this attribute instead
 * @param amplitude Strength of the Zeeman term.
 * @return void.
*/
void SystemTB::addZeeman(const double amplitude){

	arma::cx_vec zeeman_values = {amplitude, -amplitude};
	arma::cx_mat zeeman_matrix = arma::diagmat<arma::cx_mat>(arma::kron(arma::ones<arma::cx_vec>(norbitals/2), zeeman_values));

	// Identify hamiltonian slice for R=0
	uint idx = 0;
	for (unsigned int i = 0; i < Rlist.n_cols; i++){
		if (arma::norm(Rlist.col(i)) < 1E-5){
			arma::cout << Rlist.col(i) << arma::endl;
			idx = i;
			break;
		}
	}

	arma::cx_cube hamiltonianMatrices_Zeeman = *ptr_hamiltonianMatrices;
	arma::uvec indices = {idx};
	hamiltonianMatrices_Zeeman.each_slice(indices) += zeeman_matrix;

	this->hamiltonianMatrices_Zeeman_ = hamiltonianMatrices_Zeeman;

}

// /*------------------ Brillouin zone meshing, supercells & utilities ------------------*/

/**
 * Routine to compute a mesh of the first Brillouin zone using the Monkhorst-Pack algorithm.
 * @details Always returns a mesh centered at the Gamma point, for both n even or odd. 
 * @param n Number of points along one of the axis. 
 * @return void.
 */
void SystemTB::brillouinZoneMesh(int n){

	std::cout << "Creating BZ mesh... " << std::flush;

	uint32_t nk = pow(n, ndim);
	arma::mat kpoints(3, nk);
	arma::mat combinations = generateCombinations(n);
	if (n % 2 == 1){
		combinations += 1./2;
	}
	
	for (uint32_t i = 0; i < nk; i++){
		arma::colvec kpoint = arma::zeros<arma::colvec>(3);
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n)*Gbasis_.col(j);
		}
		kpoints.col(i) = kpoint;
	}

	this->kpoints_ = kpoints;
	this->nk_      = kpoints.n_cols;
	std::cout << "Done. Number of k points in BZ mesh: " << nk << std::endl;
	this->meshBZ_  = kpoints;

}

/**
 * Method to generate a k-point mesh which is a subset of the full BZ mesh. 
 * @details If n is even, we substract one so that the mesh is symmetric under inversion.
 * The reduction factor specifies is used to create a mesh of the full BZ of factor*n points,
 * from which we extract the corresponding mesh for n points centered at Gamma.
 * @param n Number of points along one axis.
 * @param factor Reduction factor of the mesh.
 * @return void.
 * */
void SystemTB::reducedBrillouinZoneMesh(const int n, const int factor){

	// First create mesh of whole BZ
	brillouinZoneMesh(n*factor);

	// Now create submesh
	uint32_t nk = pow(n, ndim);

	arma::mat kpoints(3, nk);
	arma::mat combinations = generateCombinations(n);
	if (n % 2 == 1){
		combinations += 1./2;
	}
	
	for (uint32_t i = 0; i < nk; i++){
		arma::colvec kpoint = arma::zeros<arma::colvec>(3);
		for (int j = 0; j < ndim; j++){
			kpoint += (2*combinations.row(i)(j) - n)/(2*n*factor)*Gbasis_.col(j);
		}
		kpoints.col(i) = kpoint;
	}

	this->kpoints_ = kpoints;
	this->nk_      = nk;
	this->factor_  = factor;
	std::cout << "Number of k points in submesh: " << nk << std::endl;

}

/**
 * Method to shift the center of the BZ mesh to a given point.
 * @param shift Vector to shift center of BZ mesh. 
 * @return void.
 */
void SystemTB::shiftBZ(const arma::colvec& shift){

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
	kpoints_.each_col() += shift;

}

/**
 * Auxiliary method to compute the inverse reciprocal matrix.
 * @details The reciprocal matrix is defined as the R*R.t, where R is the matrix
 * containing the reciprocal lattice vectors as row. This inverted matrix is required
 * to map any kpoint back to the original Monkhorst-Pack mesh.
 * @return void.
*/
void SystemTB::calculateInverseReciprocalMatrix(){

	arma::mat coefs = arma::zeros(ndim, ndim);
	coefs = (Gbasis.t()) * Gbasis;
	arma::mat inverse;
	try{
		inverse = arma::inv(coefs);	
	}
	catch(std::runtime_error& e){
		std::cout << "Unable to compute inverse reciprocal coefficients" << std::endl;
		throw;
	}
	this->inverseReciprocalMatrix_ = inverse;

}

/**
 * Method to determine a kpoint equivalent to another within the BZ mesh.
 * @details Given a k-point, this function finds its reciprocal coordinates to determine if 
 * it is outside of the BZ mesh. If it is, shift it by a reciprocal lattice vector so it falls
 * on a point of the mesh. Intended to be used only with mesh of full BZ (i.e. Monkhorst-Pack).
 * @param kpoint k-point to be mapped back.
 * @param ncell Number of points used in the original BZ mesh, usually equivalent to the number of cells.
 * @return Index (row) of the equivalent kpoint from the BZ mesh matrix.
 */ 
int SystemTB::findEquivalentPointBZ(const arma::colvec& kpoint, int ncell){

	if(inverseReciprocalMatrix.empty()){
		calculateInverseReciprocalMatrix();
	}
	ncell = ncell * factor_;
	arma::rowvec independentTerm = (kpoint.t())*Gbasis;
	arma::rowvec coefs = independentTerm*inverseReciprocalMatrix*(2*ncell);
	coefs = (ncell % 2 == 1) ? coefs - 1 : coefs; 
	coefs = arma::round(coefs);

	for(uint i = 0; i < coefs.n_elem; i++){
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
	for(uint i = 0; i < coefs.n_elem; i++){
		index += coefs(i)*cells_array[i];
	}

	return index;

}

/**
 * Method to generate a list of cells (stored by columns) within a sphere of specified radius.
 * @param ncell List of cells along one axis.
 * @param radius Cutoff radius of the sphere.
 * @return arma::mat List of cells within the sphere, stored by columns.
*/
arma::mat SystemTB::truncateSupercell(const int ncell, const double radius){

	arma::mat combinations = generateCombinations(ncell, true);
	std::vector<arma::colvec> cells_vector;
	for (uint32_t i = 0; i < combinations.n_rows; i++){
		arma::colvec lattice_vector = arma::zeros<arma::colvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * Rbasis.col(j);
		}
		double distance = arma::norm(lattice_vector);
		if (distance < radius + 1E-5){
			cells_vector.push_back(lattice_vector);
		}
	}

	uint32_t total_cells = cells_vector.size();
	arma::mat cells = arma::zeros(3, total_cells);
	for (uint32_t i = 0; i < total_cells; i++){
		cells.col(i) = cells_vector[i];
	}

	return cells;

}

/**
 * Method to generate a list of reciprocal cells (stored by columns) contained within a sphere of specified radius.
 * @param ncell List of cells in one direction.
 * @param radius Radius of the cutoff sphere.
 * @return arma:mat List of reciprocal cells in cartesian coordinates, stored by columns.
*/
arma::mat SystemTB::truncateReciprocalSupercell(const int ncell, const double radius){

	arma::mat combinations = generateCombinations(ncell, true);
	std::vector<arma::colvec> cells_vector;
	for (uint32_t i = 0; i < combinations.n_rows; i++){
		arma::colvec lattice_vector = arma::zeros<arma::colvec>(3);
		for (int j = 0; j < ndim; j++){
			lattice_vector += combinations.row(i)(j) * Gbasis.col(j);
		}
		double distance = arma::norm(lattice_vector);
		if (distance < radius + 1E-5){
			cells_vector.push_back(lattice_vector);
		}
	}

	uint32_t total_cells = cells_vector.size();
	arma::mat cells = arma::zeros(3, total_cells);
	for (uint32_t i = 0; i < total_cells; i++){
		cells.col(i) = cells_vector[i];
	}

	return cells;

}

/**
 * Method to obtain the lattice parameters of 2D lattices.
 * @details For simplicity, it takes a as the norm of the first Bravais vector
 * and c as the height of the 2D crystal, taking as reference the hexagonal lattice.
 * @return void.
*/
void SystemTB::extractLatticeParameters(){

	try{
		if (motif.is_empty() || Rbasis.is_empty()){
			throw "Error: Can not obtain lattice parameters (no Bravais lattice or motif)";
		}
	}
	catch (std::string e){
			std::cerr << e;
	}
	this->a_ = arma::norm(Rbasis.col(0));

	double reference_height = motif.col(0)(2);
	double c = 0;
	for (arma::uword i = 0; i < motif.n_cols; i++){
		double diff = abs(motif.col(i)(2) - reference_height);
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
 * Routine to rotate a position by 2pi/3, either on real space
 * or on reciprocal space to enforce C3 symmetry.
 * @param position Vector to rotate.
 * @return arma::rowvec Rotated vector.
 */
arma::rowvec SystemTB::rotateC3(const arma::rowvec& position){

	double theta = 2*PI/3;
	arma::mat C3rotation = {{cos(theta), -sin(theta), 0},
							{sin(theta),  cos(theta), 0},
							{         0,		   0, 1}};
	
	arma::vec rotated_position = arma::inv(C3rotation)*(position.t());

	return rotated_position.t();

}

}