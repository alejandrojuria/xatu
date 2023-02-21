#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <string>
#include <vector>

#include <xatu.hpp>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

int main(int argc, char* argv[]){

    std::cout << "Testing exciton spectrum in spinful hBN nk=20, direct diagonalization... " << std::flush;
    std::cout.setstate(std::ios_base::failbit);

    int nbands = 2;
    int nrmbands = 0;
    int ncell = 20;
    int nstates = 20;
    arma::rowvec parameters = {1., 1., 10.};
    std::string modelfile = "../models/hBN_spinful.txt";    
    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);
    xatu::Exciton bulkExciton = xatu::Exciton(config, ncell, nbands, nrmbands, parameters);
    arma::cout << "Orbitals: " << bulkExciton.orbitals << arma::endl;
    bulkExciton.setMode("realspace");

    bulkExciton.brillouinZoneMesh(ncell);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize("diag", nstates);

    std::cout.clear();
    bool testPassed = true;
    auto energies = xatu::detectDegeneracies(results.eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 8}, {6.074062, 4}, {6.164494, 4}, {6.164585, 4}};
    for(int i = 0; i < energies.size(); i++){
        if(abs(energies[i][0] - expectedEnergies[i][0]) > 1E-5){
            std::cout << "Incorrect eigval computed. " << std::flush; 
            testPassed = false;
            break;
        }
        else if(abs(energies[i][1] - expectedEnergies[i][1]) > 1E-5){
            std::cout << "Incorrect degeneracy. " << std::flush;
            testPassed = false; 
            break;
        }
    }

    if (testPassed){
        std::cout << "\033[1;32mPassed\033[0m" << std::endl;
        return 0;
    }
    else{
        std::cout << "\033[1;31mFailed\033[0m" << std::endl;
        return 1;
    }
};