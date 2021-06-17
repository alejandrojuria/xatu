#pragma once
#include <armadillo>
#include <complex>

#include "Zigzag.hpp"
#include "Exciton.hpp"

double sWaveFunction(const arma::vec&, const arma::vec&);
arma::cx_mat collapseTBcoefs(const arma::cx_mat&);
arma::cx_mat RScoefficientCalc(Exciton, const arma::cx_vec&, 
                               int);
double electronicDensity(const arma::vec&, const arma::cx_mat&, int, 
                         int, const arma::vec&);
double atomCoefficientSquared(int, int, int, const arma::cx_mat&, const arma::vec&);