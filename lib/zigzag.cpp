#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "Zigzag.hpp"

using namespace arma;
using namespace std::chrono;

Zigzag::Zigzag(int N, std::string zeeman_axis){

	this->N = N;
	this->zeeman_axis = zeeman_axis;
	initializeConstants();
	createMotif();
	initializeBlockMatrices();
	prepareHamiltonian();
	cout << "Correctly initiallized Zigzag object" << endl;

};

Zigzag::~Zigzag(){
	cout << "Destroying Zigzag object..." << endl;
};

void Zigzag::initializeConstants(){
	//// ------------ Global variable initialization ------------
	//// Lattice parameters
	a = 4.5332;
	c = 1.585;

	//// Tight-binding parameters
	// On-site energies
	Es = -10.906;
	Ep = -0.486;

	// Interaction amplitudes
	Vsss = -0.608;
	Vsps = 1.320;
	Vpps = 1.854;
	Vppp = -0.600;

	// Spin-orbit coupling 
	lambda = 1.5;

	// Zeeman term (infinitesimal, only for spin splitting)
	zeeman = 0.0;

	// Infinitesimal on-site energy to split edges
	onsiteEdge = 1E-7;

	//// Lattice vectors
	// Bravais basis
	a1 = { sqrt(3) / 2, 1.0 / 2, 0.0 };
	a1 *= a;


	a2 = { sqrt(3) / 2, -1.0 / 2, 0.0 };
	a2 *= a;

	// Motif
	tau = { a / sqrt(3), 0, -c };

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

/* Setter to access Zeeman term, which is private */
void Zigzag::setZeeman(double Zeeman){
	this->zeeman = Zeeman;
};

/* Routine to calculate the positions of the atoms inside the unit cell
Input: int N (cells in the finite direction). Output: mat motiv */
void Zigzag::createMotif(){
	motif = zeros(3, 2*(N + 1));
	motif.col(0) = arma::vec({0,0,0});
	motif.col(1) = n1;
	motif.col(2) = a1;
	motif.col(3) = a1 + n2;

	for(int n = 2; n <= N; n++){
		if(n%2 == 0){
			motif.col(2*n) = (a1 + a2)*n/2;
			motif.col(2*n + 1) = (a1 + a2)*n/2 + n1;
		}
		else{
			motif.col(2*n) = (a1 + a2)*(n-1)/2 + a1;
			motif.col(2*n + 1) = (a1 + a2)*(n-1)/2 + motif.col(3);
		};
	};
	motif = motif.t();
};


//// Matrix routines for hamiltonian initialization

/* Kronecker product to incorporate spin to a matrix (spinless interactions).
   Input: 4x4 matrix. Output: 8x8 matrix */
mat Zigzag::matrixWithSpin(const mat& matrix) {
	mat id2 = arma::eye(2, 2);
	mat M = arma::zeros(8, 8);
	// Previous implementation; a bit wierd
	M.submat( 0,0, 1,1 ) = id2*matrix(0,0);
	M.submat( 2,0, 7,1 ) = arma::kron(id2, matrix.submat( 1,0, 3,0 ));
	M.submat( 0,2, 1,7 ) = arma::kron(id2, matrix.submat( 0,1, 0,3 ));
	M.submat( 2,2, 7,7 ) = arma::kron(id2, matrix.submat( 1,1, 3,3 ));

	// New implementation: Basis is made of two spin sectors 
	// for the orbitals (up, down)
	//M.submat(0,0, 3,3) = matrix;
	//M.submat(4,4, 7,7) = matrix;

	return M;
};

/* Create tight-binding matrix for system with one s orbital and three p orbitals, based on Slater-Koster
   approximation. Input: 3x1 vector. Output: 8x8 matrix */
mat Zigzag::tightbindingMatrix(const vec& n) {
	double vNorm = norm(n);
	vec ndir = { n(0)/vNorm, n(1)/vNorm, n(2)/vNorm };
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

/* Initialize tight-binding block matrices and spin-orbit coupling of the hamiltonian for Bi bilayers
   Void function since we want multiple return (to update previously declared matrices) */
void Zigzag::initializeBlockMatrices() {

	std::complex<double> imagNum(0, 1);

	M0 = arma::zeros(4, 4);
	M0(0, 0) = Es;
	M0.submat(1, 1, 3, 3) = Ep * eye(3, 3);
	M0 = matrixWithSpin(M0);

	M1 = tightbindingMatrix(n3);		
	M2p = tightbindingMatrix(n1);
	M2m = tightbindingMatrix(n2);

	Mso = arma::zeros<cx_mat>(8, 8);
	Mso(2, 3) = -imagNum;
	Mso(3, 7) = -imagNum;
	Mso(2, 7) = 1;
	Mso(4, 5) = -1;
	Mso(4, 6) = imagNum;
	Mso(5, 6) = imagNum;

	Mso = Mso + Mso.t();
	Mso *= lambda / (3.0);

	arma::cx_mat Mzeeman = arma::zeros<cx_mat>(8,8);
	if(zeeman_axis == "z"){
		Mzeeman.diag() = arma::cx_vec({1., -1., 1., 1., 1., -1., -1., -1.,})*zeeman;
	}
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

/* Create hamiltonian matrix of the semi-infinite tight binding system (Bi ribbon).
   NB: Requires previous initialization of block matrices !!!
   Expected input: integer N (size of the system along the finite direction) Output: None (void function); updates
   previously declared bloch hamiltonian matrices */
void Zigzag::prepareHamiltonian() {
	if (N < 2) {
		std::cout << "Invalid value for N (Expected N >= 2)" << std::endl;
		return;
	};

	int hDim = 2 * (N + 1) * 8;
	H0 = arma::zeros<cx_mat>(hDim, hDim);
	Ha = arma::zeros<cx_mat>(hDim, hDim);

	mat H0block1 = arma::kron(eye(4, 4), M0 / 2);
	mat H0block2 = arma::kron(eye(4, 4), M0 / 2);

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

	Hsoc = kron(arma::eye(2 * (N + 1), 2 * (N + 1)), Mso);
	Hzeeman = kron(arma::eye<cx_mat>(2 * (N + 1), 2 * (N + 1)), Mzeeman);

	for(int i = 0; i < 8; i++){
		H0(0 + i, 0 + i) += onsiteEdge;
		H0(8*1 + i, 8*1 + i) += onsiteEdge;
		H0(2*(N+1)*8 - 8 + i, 2*(N+1)*8 - 8 + i) -= onsiteEdge;
		H0(2*(N+1)*8 - 2*8 + i, 2*(N+1)*8 - 2*8 + i) -= onsiteEdge;
	};

	this->H0 = H0;
	this->Ha = Ha;
	this->Hsoc = Hsoc;
	this->Hzeeman = Hzeeman;
};

/* Initialize Bloch hamiltonian for posterior diagonalization. Input: double k (wave number), H0, Ha 2*8*(N+1)x2*8*(N+1) matrices,
   and Hsoc 2*8*(N+1)x2*8*(N+1) complex matrix */
cx_mat Zigzag::hamiltonian(double k) {

	cx_mat h = arma::zeros<cx_mat>(H0.n_rows, H0.n_cols);
	std::complex<double> i(0, 1);
	h = H0 + Hzeeman +
	    Ha * std::exp(-i * k * a) + Ha.t() * std::exp(i * k * a) + Hsoc;

	return h;
};

/*------------------ Utilities/Additional observables ------------------*/

/* Write a text file to save calculation results and posterios graphical representation */
void Zigzag::writeEigenvaluesToFile(FILE* file, const vec& eigenval, double k) {

	fprintf(file, "%11.8f\t", k);
	for (unsigned int it = 0; it < eigenval.n_elem; it++) {
		fprintf(file, "%11.8f\t", eigenval(it));
	}
	fprintf(file, "\n");

	return;
};

/* Routine to apply the inversion operator P over a eigenstate */
cx_mat Zigzag::inversionOperator(const cx_vec& eigenvector){

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

	return P*eigenvector;
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double Zigzag::expectedSpinZValue(const arma::cx_vec& eigvec){

	int dimTB = 2*(N+1)*8;
	cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
	cx_vec spinVector = arma::kron(arma::ones(dimTB/8, 1), spinEigvalues);

	cx_vec spinEigvec = eigvec % spinVector;

	return real(arma::cdot(eigvec, spinEigvec));
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double Zigzag::expectedSpinYValue(const arma::cx_vec& eigvec){

	int dimTB = 2*(N+1)*8;
	std::complex<double> i(0,1);
	cx_mat operatorSy = arma::zeros<cx_mat>(8, 8);
	operatorSy(0,1) = -i;
	operatorSy(1,0) = i;
	operatorSy.submat(2,5, 4,7) = -i*arma::eye<cx_mat>(3,3);
	operatorSy.submat(5,2, 7,4) = i*arma::eye<cx_mat>(3,3);

	cx_mat syMatrix = 0.5*arma::kron(arma::eye<cx_mat>(dimTB/8, dimTB/8), operatorSy);

	cx_vec spinEigvec = syMatrix*eigvec;
	
	return abs(arma::cdot(eigvec, spinEigvec));
};

/* Routine to calculate the expected value of the spin projection
Sz of a TB eigenstate. */
double Zigzag::expectedSpinXValue(const arma::cx_vec& eigvec){

	int dimTB = 2*(N+1)*8;
	cx_mat operatorSx = arma::zeros<cx_mat>(8, 8);
	operatorSx(0,1) = 1;
	operatorSx(1,0) = 1;
	operatorSx.submat(2,5, 4,7) = arma::eye<cx_mat>(3,3);
	operatorSx.submat(5,2, 7,4) = arma::eye<cx_mat>(3,3);

	cx_mat sxMatrix = 0.5*arma::kron(arma::eye<cx_mat>(dimTB/8, dimTB/8), operatorSx);

	cx_vec spinEigvec = sxMatrix*eigvec;
	
	return abs(arma::cdot(eigvec, spinEigvec));
};


/* Definition of non-interacting retarded Green function */
std::complex<double> Zigzag::rGreenF(double energy, double delta, double eigEn){

	std::complex<double> i(0,1);
	return 1./((energy + i*delta) - eigEn);
};

/* Routine to calcule the density of states at a given energy,
associated to a given set of eigenvalues (e.g. bulk or edge).
NB: It is NOT normalized. */
double Zigzag::densityOfStates(double energy, double delta, const mat& energies){

		double dos = 0;
		for(int i = 0; i < (int)energies.n_rows; i++){
			for(int j = 0; j < (int)energies.n_cols; j++){
				double eigEn = energies(i,j);
				dos += -PI*imag(rGreenF(energy, delta, eigEn));
			};
		};
		dos /= energies.n_cols*a; // Divide by number of k's and length a

		return dos;
}

/* Routine to calculate and write the density of states associated 
to a given set of eigenenergies. Density of states is normalized to 1
(integral of DOS over energies equals 1)
Input: mat energies (eigenvalues), double delta (convergence parameter),
FILE* dosfile output file
Output: void (write results to output file) */
void Zigzag::writeDensityOfStates(const mat& energies, double delta, FILE* dosfile){

	double minE = energies.min();
	double maxE = energies.max();
	int nE = 2000; // Number of points in energy mesh

	vec energyMesh = arma::linspace(minE - 0.5, maxE + 0.5, nE);

	// Main loop
	double totalDos = 0;
	double dos = 0;
	// First loop over energies to normalize
	for(int n = 0; n < (int)energyMesh.n_elem; n++){
		double energy = energyMesh(n);
		double deltaE = energyMesh(1) - energyMesh(0);

		dos = densityOfStates(energy, delta, energies);
		totalDos += dos*deltaE;
	};
	// Main loop to write normalized DOS
	for(int n = 0; n < (int)energyMesh.n_elem; n++){
		double dos = 0;
		double energy = energyMesh(n);

		dos = densityOfStates(energy, delta, energies);
		dos = dos/totalDos; // Normallise
		fprintf(dosfile, "%lf\t%lf\n", energy, dos);
	};
	return;
};
