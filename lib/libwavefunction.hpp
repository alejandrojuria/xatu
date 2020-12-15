#pragma once
#include <armadillo>
#include <complex>

#include "Zigzag.hpp"
#include "Exciton.hpp"

arma::cx_mat RScoefficientCalc(Exciton&, const arma::cx_vec&, 
                               int);
double atomCoefficientSquared(int, int, int, const arma::cx_mat&, 
                              Exciton&);