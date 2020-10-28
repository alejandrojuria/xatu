#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "zigzag.hpp"
#include "libinvariant_sparse.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

// Shared variable declaration
double a, c;
double Es, Ep, Vsss, Vsps, Vpps, Vppp;
double lambda, zeeman;
vec a1, a2, tau;
vec n1, n2, n3;
vec Gamma, K, M;
mat M0, M1, M2p, M2m, Mzeeman;
cx_mat Mso;
cx_mat H0, Ha, Hsoc, Hzeeman;
cx_mat spinOperator;
sp_cx_mat PpTB_matrix, PmTB_matrix, fullTB_matrix;

int main(){

    auto start = high_resolution_clock::now();

    initializeConstants();
    // lambda = 0;
	// Number of cells in the finite direction
	int N = 15;
	int nOrbitals = 8;
    int dimTB = 2*(N+1)*nOrbitals;
	int Ncell = 1;
    int nk = 2*Ncell + 1;
	vec kpoints = arma::linspace(-PI / a, PI / a, nk);
	bool plotSpinMix = false;
	int dimFull = dimTB*nk;

	initializeBlockMatrices();
	prepareHamiltonian(N);

	// Calculation to plot spin mixing in occupied states
	if(plotSpinMix == true){
		std::string filename = "spin_eigenvalues";
		FILE* textfile = fopen(filename.c_str(), "w");

		cx_cube eigvecKStack(dimTB, dimTB, nk);
		cx_mat eigenvec;
		vec eigenval;
		for (int i = 0; i < nk; i++) {
			cx_mat h = hamiltonian(kpoints(i), H0, Ha, Hsoc);

			arma::eig_sym(eigenval, eigvecKStack.slice(i), h);

			spinOperator = PspinMatrix(eigvecKStack.slice(i));
			arma::eig_sym(eigenval, eigenvec, spinOperator);
			writeEigenvaluesToFile(textfile, eigenval, kpoints[i]);
		};
		fclose(textfile);
	}; 

	calculateFullTBMatrices(kpoints, N, nOrbitals);

	mat motif = createMotiv(N);
	sp_cx_mat eX = positionMatrix(0, motif, nOrbitals, Ncell, N);
	cout << "exp(X) computed " << endl;
	sp_cx_mat eY = positionMatrix(1, motif, nOrbitals, Ncell, N);
	cout << "exp(Y) computed " << endl;

	sp_cx_mat sectorPlus = projectorMatrix(0, fullTB_matrix, 
												PpTB_matrix);
	cout << "Projector + computed" << endl;

	sp_cx_mat sectorMinus = projectorMatrix(1, fullTB_matrix, 
												PmTB_matrix);
	cout << "Projector - computed" << endl;

	double bottPlus = sectorBottIndex(eX, eY, sectorPlus);
	cout << "Bott index + computed" << endl;
	double bottMinus = sectorBottIndex(eX, eY, sectorMinus);
	cout << "Bott index - computed" << endl;

	double spinBott = bottPlus - bottMinus;

	cout << "Ncell = " << nk << endl;
	cout << "Spin Bott index is: " << spinBott << endl;

	// Orthonormality check
    // cout << arma::cdot(eigenvec.col(0), eigenvec.col(0)) << endl;

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};