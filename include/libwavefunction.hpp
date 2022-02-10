#pragma once
#include <armadillo>
#include <complex>

#include "System.hpp"
#include "GExciton.hpp"

arma::cx_mat RScoefficientCalc(GExciton&, const arma::cx_vec&, 
                               int);
double atomCoefficientSquared(int, const arma::rowvec&, const arma::rowvec&, 
                              const arma::cx_mat&, GExciton&);
double fourierTransformExciton(const arma::cx_vec&, GExciton&, const arma::rowvec&, const arma::rowvec&);

double realSpaceWavefunction(GExciton&, const arma::cx_vec&, int, int, const arma::rowvec& eCell, const arma::rowvec& hCell);
std::complex<double> densityMatrixElement(GExciton&, int, int, int, int);

std::complex<double> densityMatrix(GExciton&, const arma::cx_vec&, int, int);
std::complex<double> densityMatrixK(int, GExciton&, const arma::cx_vec&, int, int);