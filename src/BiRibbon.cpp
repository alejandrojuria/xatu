#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "xatu/BiRibbon.hpp"

using namespace arma;
using namespace std::chrono;

namespace xatu {

/**
 * Bismuth zigzag nanoribbon constructor.
 * @param N Number of dimers in the unit cells.
 * @param zeeman_axis String to specify axis to apply magnetic term ('x', 'y' or 'z'). 
 */
BiRibbon::BiRibbon(int N, std::string zeeman_axis){

	this->N = N;
	this->zeeman_axis = zeeman_axis;
	
	initializeConstants();
	calculateReciprocalLattice();
	createMotif();
	initializeBlockMatrices();
	prepareHamiltonian();
	
	cout << "Correctly initiallized Ribbon object" << endl;
};

/**
 * Bismuth ribbon destructor. 
 */
BiRibbon::~BiRibbon(){
	cout << "Destroying Ribbon object..." << endl;
};

/**
 * Routine to initialize the parameters of the model. To be called from the constructor.
 * @return void 
 */
void BiRibbon::initializeConstants(){
	//// ------------ Global variable initialization ------------

    // System class attributes
	this->ndim_      = 1;
    this->orbitalsPerSpecies_  = arma::urowvec{8};
    this->filling_   = 5;

	//// Lattice parameters
	this->a = 4.5332;
	this->c = 1.585;

	//// Tight-binding parameters
	// On-site energies
	this->Es = -10.906;
	this->Ep = -0.486;

	// Interaction amplitudes
	this->Vsss = -0.608;
	this->Vsps = 1.320;
	this->Vpps = 1.854;
	this->Vppp = -0.600;

	// Spin-orbit coupling 
	this->lambda = 1.5;

	// Zeeman term (infinitesimal, only for spin splitting)
	this->zeeman = 1E-7;

	// Infinitesimal on-site energy to split edges
	this->onsiteEdge = 0.0;

	//// Lattice vectors
	// Bulk Bravais basis
	a1 = { sqrt(3) / 2, 1.0 / 2, 0.0 };
	a1 *= a;
	a2 = { sqrt(3) / 2, -1.0 / 2, 0.0 };
	a2 *= a;
    // Ribbon Bravais lattice
	this->bravaisLattice_ = arma::mat{ a1 - a2 };

	// Motif
	tau = { a / sqrt(3), 0, -c};

	// First neighbours
	n1 = a1 - tau;
	n2 = a2 - tau;
	n3 = tau;

	//// High symmetry points of IBZ
	Gamma = { 0, 0, 0 };
	K = { 1.0 / sqrt(3), 1.0 / 3, 0 };
	K = 2 * PI * K / a;
	M = { 1.0 / sqrt(3), 0, 0 };
	M = 2 * PI * M / a;
};

/** 
 * Setter for Zeeman term amplitude.
 * @param zeeman Zeeman amplitude.
 * @return void
 */
void BiRibbon::setZeeman(double zeeman){
	this->zeeman = zeeman;
};

/** 
 * Routine to calculate the positions of the atoms inside the unit cell of the ribbon.
 * @return void
 */
void BiRibbon::createMotif(){
	arma::mat motif = arma::zeros(2*(N + 1), 3);
	motif.row(0) = arma::rowvec({0,0,0});
	motif.row(1) = n1;
	motif.row(2) = a1;
	motif.row(3) = a1 + n2;

	for(int n = 2; n <= N; n++){
		if(n%2 == 0){
			motif.row(2*n) = (a1 + a2)*n/2;
			motif.row(2*n + 1) = (a1 + a2)*n/2 + n1;
		}
		else{
			motif.row(2*n) = (a1 + a2)*(n-1)/2 + a1;
			motif.row(2*n + 1) = (a1 + a2)*(n-1)/2 + a1 + n2;
		};
	};

	// Add species (zeros)
	motif.insert_cols(3, 1);

    // Initiallize remaining System attributes
    this->motif_      = motif;
    this->natoms_     = motif.n_rows;
    int basisdim = 0;
    for(int i = 0; i < natoms; i++){
        int species = this->motif.row(i)(3);
        basisdim += orbitalsPerSpecies(species);
    }
    this->basisdim_   = basisdim;
    this->filling_ = (int)natoms*filling;
	this->fermiLevel_  = fermiLevel - 1; // Overwrite filling to match fermiLevel
};


//// Matrix routines for hamiltonian initialization

/** 
 * Kronecker product to incorporate spin to a matrix by doubling it.
 * @param matrix Matrix without spin in its basis (polarized).
 * @return Spinful matrix.
 */
mat BiRibbon::matrixWithSpin(const mat& matrix) {
	mat id2 = arma::eye(2, 2);
	mat M = arma::zeros(8, 8);

	// New implementation: Basis is made of two spin sectors 
	// for the orbitals (up, down)
	M = arma::kron(matrix, arma::eye(2, 2));

	return M;
};

/** 
 * Create tight-binding matrix with the hoppings for system with atoms with one s orbital and three p orbitals, 
 * based on Slater-Koster approximation. 
 * @param n Displacemente vector connecting the two atoms involved.
 * @return Hopping matrix.
 */
mat BiRibbon::tightbindingMatrix(const rowvec& n) {
	
	double vNorm = arma::norm(n);
	arma::vec ndir = { n(0)/vNorm, n(1)/vNorm, n(2)/vNorm };
	mat M = arma::zeros(4, 4);

	M(0, 0) = Vsss;
	for (int i = 0; i < 3; i++) {
		M(0, i + 1) = ndir(i) * Vsps;
		M(i + 1, 0) = -M(0, i + 1);
		M(i + 1, i + 1) = ndir(i) * ndir(i) * Vpps + (1 - ndir(i) * ndir(i)) * Vppp;
		for (int j = i + 1; j < 3; j++) {
			M(i + 1, j + 1) = ndir(i) * ndir(j) * (Vpps - Vppp);
			M(j + 1, i + 1) = ndir(i) * ndir(j) * (Vpps - Vppp);
		};
	};


	return matrixWithSpin(M);
};

/**
 * Initialize tight-binding block matrices and spin-orbit coupling of the hamiltonian for Bi bilayers
 * @return void
 */
void BiRibbon::initializeBlockMatrices() {

	std::complex<double> imagNum(0, 1);

	M0 = arma::zeros(4, 4);
	M0(0, 0) = Es;
	M0.submat(1, 1, 3, 3) = Ep * arma::eye(3, 3);
	M0 = matrixWithSpin(M0);

	M1 = tightbindingMatrix(n3);		
	M2p = tightbindingMatrix(n1);
	M2m = tightbindingMatrix(n2);

	// This has to be changed to match new ordering (up, down)
	Mso = arma::zeros<cx_mat>(8, 8);
	Mso(2, 4) = -imagNum;
	Mso(2, 7) = 1;
	Mso(3, 5) = imagNum;
	Mso(3, 6) = -1;
	Mso(4, 7) = -imagNum;
	Mso(5, 6) = -imagNum;

	Mso = Mso + Mso.t();
	Mso *= lambda / (3.0);

	arma::cx_mat Mzeeman = arma::zeros<cx_mat>(8,8);
	if(zeeman_axis == "z"){
		Mzeeman.diag() = arma::cx_vec({1., -1., 1., -1., 1., -1., 1., -1.,})*zeeman;
	}
	// Wrong implementation (for old basis ordering)
	else if(zeeman_axis == "x"){
		arma::cx_mat Sx(2,2);
		Sx(0,1) = 1;
		Sx(1,0) = 1;
		Mzeeman.submat(0,0, 1,1) = Sx;
		Mzeeman.submat(2,2, 7,7) = arma::kron(Sx, arma::eye<cx_mat>(3,3));
 		Mzeeman *= zeeman;
	}
	else if(zeeman_axis == "y"){
		std::complex<double> i(0,1);
		arma::cx_mat Sy(2,2);
		Sy(0,1) = -i;
		Sy(1,0) = i;
		Mzeeman.submat(0,0, 1,1) = Sy;
		Mzeeman.submat(2,2, 7,7) = arma::kron(Sy, arma::eye<cx_mat>(3,3));
		Mzeeman *= zeeman;
	}
	else{
		cout << "Incorrect Zeeman axis given, defaulting Zeeman term to zero" << endl;
	};

	this->M0 = M0;
	this->M1 = M1;
	this->M2p = M2p;
	this->M2m = M2m;
	this->Mso = Mso;
	this->Mzeeman = Mzeeman;
};

/** Create hamiltonian matrix of the semi-infinite tight binding system (Bi ribbon).
 * NB: Requires previous initialization of block matrices; updates previously declared bloch hamiltonian matrices.
 * @return void
 */
void BiRibbon::prepareHamiltonian() {
	if (N < 2) {
		std::cout << "Invalid value for N (Expected N >= 2)" << std::endl;
		return;
	};

	int hDim = 2 * (N + 1) * 8;
	arma::cx_mat H0 = arma::zeros<cx_mat>(hDim, hDim);
	arma::cx_mat Ha = arma::zeros<cx_mat>(hDim, hDim);

	mat H0block1 = arma::kron(arma::eye(4, 4), M0 / 2);
	mat H0block2 = arma::kron(arma::eye(4, 4), M0 / 2);

	mat HaBlock1 = arma::zeros(4 * 8, 4 * 8);
	mat HaBlock2 = arma::zeros(4 * 8, 4 * 8);

	H0block1.submat(0, 1 * 8, 8 - 1, 2 * 8 - 1) = M2p;
	H0block1.submat(8, 2 * 8, 2 * 8 - 1, 3 * 8 - 1) = M1;
	H0block1.submat(2 * 8, 3 * 8, 3 * 8 - 1, 4 * 8 - 1) = M2m;

	H0block2.submat(0, 1 * 8, 8 - 1, 2 * 8 - 1) = M2m;
	H0block2.submat(8, 2 * 8, 2 * 8 - 1, 3 * 8 - 1) = M1;
	H0block2.submat(2 * 8, 3 * 8, 3 * 8 - 1, 4 * 8 - 1) = M2p;

	H0block1 = H0block1 + H0block1.t();
	H0block2 = H0block2 + H0block2.t();

	HaBlock1.submat(8, 0, 2 * 8 - 1, 8 - 1) = M2m.t();
	HaBlock1.submat(2 * 8, 3 * 8, 3 * 8 - 1, 4 * 8 - 1) = M2p;

	HaBlock2.submat(0, 8, 8 - 1, 2 * 8 - 1) = M2p;
	HaBlock2.submat(3 * 8, 2 * 8, 4 * 8 - 1, 3 * 8 - 1) = M2m.t();

	int i = 4 * 8; // N=1
	int j;
	H0.submat(0, 0, i - 1, i - 1) = conv_to<cx_mat>::from(H0block1);
	Ha.submat(0, 0, i - 1, i - 1) = conv_to<cx_mat>::from(HaBlock1);
	for (int m = 0; m < N - 1; m++) {
		i -= 2 * 8;
		j = i + 4 * 8;
		if (m % 2 == 0) {
			H0.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(H0block2);
			Ha.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(HaBlock2);
		}
		else {
			H0.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(H0block1);
			Ha.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(HaBlock1);
		};
		i = j;
	};

	arma::cx_mat Hsoc = kron(arma::eye(2 * (N + 1), 2 * (N + 1)), Mso);
	arma::cx_mat Hzeeman = arma::kron(arma::eye<cx_mat>(2 * (N + 1), 2 * (N + 1)), Mzeeman);

	for(int i = 0; i < 8; i++){
		H0(0 + i, 0 + i) += onsiteEdge;
		H0(8*1 + i, 8*1 + i) += onsiteEdge;
		H0(2*(N+1)*8 - 8 + i, 2*(N+1)*8 - 8 + i) -= onsiteEdge;
		H0(2*(N+1)*8 - 2*8 + i, 2*(N+1)*8 - 2*8 + i) -= onsiteEdge;
	};

	this->ncells_ = 3;
	this->unitCellList_ = arma::zeros(ncells, 3);
	this->unitCellList_.row(1) = bravaisLattice.row(0);
	this->unitCellList_.row(2) = -bravaisLattice.row(0);
    this->hamiltonianMatrices_ = arma::zeros<arma::cx_cube>(basisdim, basisdim, ncells);
    this->hamiltonianMatrices_.slice(0) = H0 + Hsoc + Hzeeman;
    this->hamiltonianMatrices_.slice(1) = Ha;
	this->hamiltonianMatrices_.slice(2) = Ha.t();

};


/**
 * Routine to apply the inversion operator P over a given state.
 * @param state State over which we want to act.
 * @return P|state>
 */
cx_mat BiRibbon::inversionOperator(const cx_vec& state){

	int dimTB = 2*(N+1)*8;
	cx_mat P = arma::zeros<cx_mat>(dimTB/8, dimTB/8);
	for(int i = 0; i < dimTB/8; i++){
		P(dimTB/8 - 1 - i, i) = 1.;
	};
	cx_mat spinOperator = arma::zeros<cx_mat>(8,8);
	spinOperator(0, 1) = -1;
	spinOperator(1, 0) = -1;
	spinOperator.submat(2,5, 4,7) = -arma::eye<cx_mat>(3,3);
	spinOperator.submat(5,2, 7,4) = -arma::eye<cx_mat>(3,3);
	P = arma::kron(P, spinOperator);

	return P*state;
};


/**
 * Method to modify the onsite energies according to a electric field perpendicular to the periodic direction.
 * @param amplitude Amplitude of applied electric field.
 * @return void
 */
void BiRibbon::applyElectricField(double amplitude){

	arma::cx_mat onsiteField = arma::eye<arma::cx_mat>(8, 8)*amplitude;
	for(int i = 0; i < natoms; i++){
		hamiltonianMatrices_.slice(0).submat(i*8, i*8, (i + 1)*8 - 1, (i + 1)*8 - 1) += onsiteField*motif.row(i)(0);
	}
};

/**
 * Method to introduce additional onsite energy at the edges to break the edge state degeneracy 
 * and identify them. 
 * @param energy Value of onsite edge energies.
 * @return void
 */
void BiRibbon::offsetEdges(double energy){

	// arma::uvec indices = {(unsigned)natoms - 2, (unsigned)natoms - 1};
	arma::uvec indices = {0, 1};
	arma::cx_mat onsiteField = arma::eye<arma::cx_mat>(8, 8)*energy;
	for (unsigned int i : indices){
		hamiltonianMatrices_.slice(0).submat(i*8, i*8, (i + 1)*8 - 1, (i + 1)*8 - 1) += onsiteField;
	}
};

/**
 * Method to introduce a substrate that modifies the onsite energies of the atoms of the lower sublattice.
 * @param energy Value of onsite sublattice energy.
 * @return void
 */
void BiRibbon::addSubstrate(double energy){

	arma::cx_mat onsiteField = arma::eye<arma::cx_mat>(8, 8)*energy;
	for(unsigned int i = 0; i < natoms/2; i++){
		hamiltonianMatrices_.slice(0).submat(2*i*8, 2*i*8, (2*i + 1)*8 - 1, (2*i + 1)*8 - 1) += onsiteField;
	}
}


}