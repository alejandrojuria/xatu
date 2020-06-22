#pragma once
#include <armadillo>
#include <string>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

// Shared variable declaration
extern double a, c;
extern double Es, Ep, Vsss, Vsps, Vpps, Vppp;
extern double lambda, zeeman, onsiteEdge;
extern arma::vec a1, a2, tau;
extern arma::vec n1, n2, n3;
extern arma::vec Gamma, K, M;
extern arma::mat M0, M1, M2p, M2m, Mzeeman;
extern arma::cx_mat Mso;
extern arma::cx_mat H0, Ha, Hsoc, Hzeeman;

/* Global variable initialization */
void initializeConstants();

/* Matrix routines for hamiltonian initialization */
arma::mat matrixWithSpin(const arma::mat& matrix);
arma::mat tightbindingMatrix(const arma::vec& n);
arma::mat createMotiv(int);
void initializeBlockMatrices();
void prepareHamiltonian(int N);
arma::cx_mat hamiltonian(double k, const arma::cx_mat& H0, const arma::cx_mat& Ha, const arma::cx_mat& Hsoc);

/* Some utilities/extra information */
void writeEigenvaluesToFile(FILE* file, const arma::vec& eigenval, double k);
arma::cx_mat inversionOperator(const arma::cx_vec&, int);

/* Expected value of spin components */
double expectedSpinZValue(const arma::cx_vec&, int);
double expectedSpinYValue(const arma::cx_vec&, int);
double expectedSpinXValue(const arma::cx_vec&, int);

/* Routines for DoS calculation */
std::complex<double> rGreenF(double, double, double);
double densityOfStates(double, double, const arma::mat&);
void writeDensityOfStates(const arma::mat&, double, FILE*);