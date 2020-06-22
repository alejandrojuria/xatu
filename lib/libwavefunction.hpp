#pragma once
#include <armadillo>
#include <complex>

double sWaveFunction(const arma::vec&, const arma::vec&);
arma::cx_mat collapseTBcoefs(const arma::cx_mat&);
arma::cx_mat RScoefficientCalc(const arma::cx_vec&, 
                               const arma::vec&, int, int, int);
double electronicDensity(const arma::vec&, const arma::cx_mat&, int, 
                         int, const arma::vec&);
double atomCoefficientSquared(int, int, int, const arma::cx_mat&, const arma::vec&);