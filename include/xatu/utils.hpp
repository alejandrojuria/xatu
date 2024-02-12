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
void printHeader();

/* Input */
arma::vec readVectorFromFile(std::string);


/* Routines for DoS calculation */
std::complex<double> rGreenF(double, double, double);
double densityOfStates(double, double, const arma::mat&);
void writeDensityOfStates(const arma::mat&, double, FILE*);

/* Routines for tests */
template<typename T>
double array2hash(const T& array){
    T copyArray(array);
    for (unsigned int i = 0; i < array.n_rows; i++){
        for(unsigned int j = 0; j < array.n_cols; j++){
            copyArray(i, j) *= i*j + j + 1;
        }   
    }

    auto realArray = arma::abs(array);
    double hash    = (arma::accu(realArray)+ 
                      arma::accu(array != 0) + 
                      arma::accu(arma::real(copyArray)*1.5) + 
                      arma::accu(arma::imag(copyArray)*1.5) + 
                      arma::accu(arma::imag(copyArray) % arma::real(copyArray)))/array.n_elem;
                    
    return hash;
};


}