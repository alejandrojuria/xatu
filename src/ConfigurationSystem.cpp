#include "xatu/ConfigurationSystem.hpp"

namespace xatu {

/**
 * Method to parse several three-dimensional vectors into a matrix, where they are stored by columns.
 * @details This method is intented to be used with the Bravais vectors. 
 * @param vectors Array with the strings encoding the vectors to be parsed.
 * @return Matrix with the parsed vectors as rows.
 */
arma::mat ConfigurationSystem::parseVectors(const std::vector<std::string>& vectors){
    arma::mat colVectors = arma::zeros(3, vectors.size());
    for (uint i = 0; i < vectors.size(); i++){
        std::string line = vectors[i];
        std::vector<double> colVect = parseLine<double>(line);
        if (colVect.size() != 3){ 
            throw std::logic_error("Bravais vectors must have three components"); }

        colVectors.col(i) = arma::colvec(colVect);
    }

    return colVectors;
}

/**
 * Auxiliary method to print the contents of the configuration struct.
 * @details Useful for debugging. 
 */
void ConfigurationSystem::printConfiguration(const bool printH, const bool printS) {

    arma::cout << "------------------------------ Common configuration attributes ------------------------------" << arma::endl;

    arma::cout << "Lattice dimensionality: " << arma::endl; 
    arma::cout << ndim << arma::endl;

    arma::cout << "Basis of direct lattice vectors: " << arma::endl;
    arma::cout << Rbasis << arma::endl;

    arma::cout << "Atomic motif: " << arma::endl;
    arma::cout << motif << arma::endl;

    arma::cout << "Number of filled bands: " << arma::endl;
    arma::cout << filling << arma::endl;

    arma::cout << "List of direct lattice vectors for which the Hamiltonian is given: " << arma::endl;
    arma::cout << Rlist << arma::endl;

    arma::cout << "Number of unit cells considered for the Hamiltonian: " << arma::endl;
    arma::cout << ncells << arma::endl;

    arma::cout << "Number of orbitals per species: " << arma::endl;
    arma::cout << orbitalsPerSpecies << arma::endl;

    if(printH){
        arma::cout << "Hamiltonian matrices H(R): " << arma::endl;
        arma::cout << hamiltonianMatrices << arma::endl;
    }

    if(printS){
        arma::cout << "Overlap matrices S(R): " << arma::endl;
        arma::cout << overlapMatrices << arma::endl;
    }

    arma::cout << "---------------------------------------------------------------------------------------------" << arma::endl;

}


}