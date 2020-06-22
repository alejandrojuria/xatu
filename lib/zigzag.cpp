#include <iostream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "zigzag.hpp"

#define PI 3.14159265359

using namespace arma;
using namespace std::chrono;

// Shared variable definition
double a, c;
double Es, Ep, Vsss, Vsps, Vpps, Vppp;
double lambda;
vec a1, a2, tau;
vec n1, n2, n3;
vec Gamma, K, M;
mat M0, M1, M2p, M2m;
cx_mat Mso;
cx_mat H0, Ha, Hsoc;

auto start = high_resolution_clock::now();

void initializeConstants() {
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

//// Matrix routines for hamiltonian initialization

/* Kronecker product to incorporate spin to a matrix (spinless interactions).
   Input: 4x4 matrix. Output: 8x8 matrix */
mat matrixWithSpin(const mat& matrix) {
	mat id2 = eye(2, 2);
	mat M = arma::zeros(8, 8);
	M.submat( 0,0, 1,1 ) = id2*matrix(0,0);
	M.submat( 2,0, 7,1 ) = kron(id2, matrix.submat( 1,0, 3,0 ));
	M.submat( 0,2, 1,7 ) = kron(id2, matrix.submat( 0,1, 0,3 ));
	M.submat( 2,2, 7,7 ) = kron(id2, matrix.submat( 1,1, 3,3 ));

	return M;
};

/* Create tight-binding matrix for system with one s orbital and three p orbitals, based on Slater-Koster
   approximation. Input: 3x1 vector. Output: 8x8 matrix */
mat tightbindingMatrix(const vec& n) {
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
void initializeBlockMatrices() {

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
	Mso *= lambda / 3.0;

	return;
};

/* Create hamiltonian matrix of the semi-infinite tight binding system (Bi ribbon).
   Expected input: integer N (size of the system along the finite direction) Output: None (void function); updates
   previously declared bloch hamiltonian matrices */
void prepareHamiltonian(int N) {
	if (N < 2) {
		std::cout << "Invalid value for N (Expected N >= 2)" << std::endl;
		return;
	};
	int hDim = 2 * (N + 1) * 8;
	H0 = arma::zeros<cx_mat>(hDim, hDim);
	Ha = arma::zeros<cx_mat>(hDim, hDim);

	mat H0block1 = kron(eye(4, 4), M0 / 2);
	mat H0block2 = kron(eye(4, 4), M0 / 2);

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
	H0.submat(0, 0, i - 1, i - 1) = conv_to<cx_mat>::from(H0block1);
	Ha.submat(0, 0, i - 1, i - 1) = conv_to<cx_mat>::from(HaBlock1);
	for (int m = 1; m < N; m++) {
		i -= 2 * 8;
		int j = i + 4 * 8;
		if (m % 2 == 0) {
			H0.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(H0block1);
			Ha.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(HaBlock1);
		}
		else {
			H0.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(H0block2);
			Ha.submat(i, i, j - 1, j - 1) = conv_to<cx_mat>::from(HaBlock2);
		};
		i = j;
	};

	Hsoc = kron(eye(2 * (N + 1), 2 * (N + 1)), Mso);

	return;
};

/* Initialize Bloch hamiltonian for posterior diagonalization. Input: double k (wave number), H0, Ha 2*8*(N+1)x2*8*(N+1) matrices,
   and Hsoc 2*8*(N+1)x2*8*(N+1) complex matrix */
cx_mat hamiltonian(double k, const cx_mat& H0, const cx_mat& Ha, const cx_mat& Hsoc) {

	cx_mat h = arma::zeros<cx_mat>(H0.n_rows, H0.n_cols);
	std::complex<double> i(0, 1);
	h = H0 + Ha * std::exp(-i * k * a) + Ha.t() * std::exp(i * k * a) + Hsoc;

	return h;
};

/* Write a text file to save calculation results and posterios graphical representation */
void writeEigenvaluesToFile(FILE* file, const vec& eigenval, double k) {

	fprintf(file, "%11.8f\t", k);
	for (unsigned int it = 0; it < eigenval.n_elem; it++) {
		fprintf(file, "%11.8f\t", eigenval(it));
	}
	fprintf(file, "\n");
};


int main() {

	initializeConstants();
	// Number of cells in the finite direction
	int N = 15;

	initializeBlockMatrices();
	prepareHamiltonian(N);

	int nk = 100;
	vec kpoints = linspace(-PI / a, PI / a, nk);

	std::string filename = "test_file";
	FILE* textfile = fopen(filename.c_str(), "w");

	for (int i = 0; i < nk; i++) {
		cx_mat h = hamiltonian(kpoints[i], H0, Ha, Hsoc);
		cx_mat eigenvec;
		vec eigenval;

		eig_sym(eigenval, eigenvec, h);
		writeEigenvaluesToFile(textfile, eigenval, kpoints[i]);
	}

	fclose(textfile);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};
