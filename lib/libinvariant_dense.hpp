#pragma once
#include <armadillo>
#include <complex>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

extern arma::cx_mat PpTB_matrix, PmTB_matrix, fullTB_matrix;

arma::cx_mat PspinMatrix(const arma::cx_mat& eigvectors);
void calculateFullTBMatrices(const arma::vec& kpoints, int N, int nOrb);
arma::cx_mat positionMatrix(int posVar, const arma::mat& motif, 
                               int nOrbitals, int Ncell, int N);
arma::cx_mat projectorMatrix(int sign, const arma::cx_mat& tbStack,
                                const arma::cx_mat& Psector, int Ncell);
double sectorBottIndex(const arma::cx_mat& eX, const arma::cx_mat& eY, 
                       const arma::cx_mat& Psector);