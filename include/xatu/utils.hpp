#pragma once
#include <armadillo>
#include "xatu/Result.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

/* Output */
void writeVectorToFile(arma::vec, FILE*);
void writeVectorToFile(arma::rowvec, FILE*);
void writeVectorsToFile(const arma::mat&, FILE*, std::string mode = "row");
std::vector<std::vector<double>>  detectDegeneracies(const arma::vec&, int, int);
void printEnergies(Result, int n = 8, int precision = 6);

/* Input */
arma::vec readVectorFromFile(std::string);


/* Routines for DoS calculation */
std::complex<double> rGreenF(double, double, double);
double densityOfStates(double, double, const arma::mat&);
void writeDensityOfStates(const arma::mat&, double, FILE*);

}