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
void printHeader();

/* Input */
arma::vec readVectorFromFile(std::string);

/* Routines for DoS calculation */
std::complex<double> rGreenF(double, double, double);
double densityOfStates(double, double, const arma::mat&);
void writeDensityOfStates(const arma::mat&, double, FILE*);

/* Matrix properties */
bool checkIfTriangular(const arma::cx_mat&);

/** Routine to pretty print the eigenenergies from the BSE calculation. Computes number of degenerate states 
* corresponding to each energy level. 
* @param results Unique pointer to Result object
* @param n Number of energies 
* @param precision Number of decimals (sets degeneracy threshold)
**/
template<typename T>
void printEnergies(const std::unique_ptr<T>& results, int n = 8, int precision = 6){

    // Print header
    printf("+---------------+-----------------------------+-----------------------------+\n");
    printf("|       N       |          Eigval (eV)        |          Degeneracy         |\n");
    printf("+---------------+-----------------------------+-----------------------------+\n");

    std::vector<std::vector<double>> pairs = detectDegeneracies(results->eigval, n, precision);
    int it = 1;

    for(auto pair : pairs){
        double energy  = pair[0];
        int degeneracy = (int)pair[1];

        printf("|%15d|%29.*lf|%29d|\n", it, precision, energy, degeneracy);
        printf("+---------------+-----------------------------+-----------------------------+\n");

        it++;
    }
}

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