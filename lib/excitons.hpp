#pragma once
#include <armadillo>
#include <complex>
#include <omp.h>
#include <stdlib.h>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

extern arma::cx_mat HBS;
extern arma::mat HK;
extern arma::mat states;

// Routines to initialize BSE elements

void STVH0(double, double*);
arma::cx_mat eigenstatesH0(double);
double potential(double);
std::complex<double> realfourierTrans(double, int Ncell);
std::complex<double> fourierTrans(double, int Ncell);
// Implementation of approximated direct and exchange terms
std::complex<double> tDirect(std::complex<double>,
                             const arma::cx_vec&, 
                             const arma::cx_vec&,
                             const arma::cx_vec&, 
                             const arma::cx_vec&, int);
std::complex<double> tExchange(std::complex<double>, 
                               const arma::cx_vec&, 
                               const arma::cx_vec&,
                               const arma::cx_vec&, 
                               const arma::cx_vec&, int);

arma::mat initializePotentialMatrix(int, const arma::mat&);
// Implementation of exact direct and exchange terms
std::complex<double> exactInteractionTerm(const arma::cx_vec&, 
                                     const arma::cx_vec&,
                                     const arma::cx_vec&, 
                                     const arma::cx_vec&, 
                                     const arma::vec&,
                                     const arma::mat&,
                                     const arma::mat&,
                                     int);

std::complex<double> selfEnergy(const arma::cx_cube&, 
                                const arma::rowvec&,
                                const arma::mat&,
                                const arma::vec& kpoints,
                                int, int);

// Routines for BSE matrix initialization

arma::mat createBasis(int, double, const arma::vec&, int, int nEdgeStates = 0);
arma::mat createSOCBasis(int, double, const arma::vec&, int);
void fixBandCrossing(arma::vec&, arma::cx_mat&, int);
int determineKIndex(double k, const arma::vec& kpoints);
arma::cx_cube atomicGCoefs(const arma::cx_cube&, const arma::mat&, const arma::vec&, int);
void BShamiltonian(int, int, const arma::mat&, const arma::vec&, std::string ordering = "kpoints");
arma::vec computeEnergies(const arma::cx_vec&, const arma::cx_mat&, const arma::mat&);

// Exciton information
arma::cx_vec spinX(const arma::cx_vec&, const arma::mat&, const arma::vec&, int);

// Routines to compute Fermi Golden Rule
arma::cx_vec ehPairCoefs(double, const arma::vec&, const arma::mat&);
double fermiGoldenRule(const arma::cx_vec&, double, int, double, int, int, const arma::vec&);