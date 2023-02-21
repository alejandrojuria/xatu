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

    std::cout << "Testing exciton spectrum in hBN nk=40, davidson method... " << std::flush;
    std::cout.setstate(std::ios_base::failbit);

    int nbands = 1;
    int nrmbands = 0;
    int ncell = 40;
    int nstates = 8;
    arma::rowvec parameters = {1., 1., 10.};
    std::string modelfile = "../models/hBN.txt";    
    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);
    xatu::Exciton bulkExciton = xatu::Exciton(config, ncell, nbands, nrmbands, parameters);
    arma::cout << "Orbitals: " << bulkExciton.orbitals << arma::endl;
    bulkExciton.setMode("realspace");

    bulkExciton.brillouinZoneMesh(ncell);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize("davidson", nstates);

    std::cout.clear();
    bool testPassed = true;
    auto energies = xatu::detectDegeneracies(results.eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335687, 2}, {6.073800, 1}, {6.164057, 2}, {6.172253, 1}, {6.351066, 2}};
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