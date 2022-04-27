#include <armadillo>
#include <fstream>

#include "utils.hpp"

void writeVectorsToFile(const arma::mat& vectors, FILE* textfile){
    for(unsigned int i = 0; i < vectors.n_rows; i++){
        arma::rowvec vector = vectors.row(i);
        fprintf(textfile, "%10.7lf\t%10.7lf\t%10.7lf\n", vector(0), vector(1), vector(2));
    }
}

arma::vec readVectorFromFile(std::string filename){
    std::ifstream file(filename.c_str());
    std::string line;
    std::vector<double> vector;
    double value;
    while(std::getline(file, line)){
        std::istringstream iss(line);
        while(iss >> value){
            vector.push_back(value);
        };
    };

    arma::vec coefs(vector);
    return coefs;
};
