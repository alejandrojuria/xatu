#include <complex>
#include "SystemConfiguration.hpp"

SystemConfiguration::SystemConfiguration(){
    throw std::invalid_argument("SystemConfiguration must be called with one argument (filename)");
};

SystemConfiguration::SystemConfiguration(std::string filename) : ConfigurationBase(filename){
    expectedArguments = {"dimension", "bravaislattice", "motif", "norbitals", "bravaisvectors", "hamiltonian"};
    parseContent();
    checkContentCoherence();
}

void SystemConfiguration::parseContent(){
    if (contents.empty()){
        throw std::logic_error("File contents must be extracted first");
    }
    for(int i = 0; i < foundArguments.size(); i++){
        std::string arg = foundArguments[i];
        auto content = contents[arg];

        if (arg == "dimension"){
            if(content.size() != 1){throw std::logic_error("Expected only one line for 'Dimension'");}
            std::string line = content[0];
            systemInfo.ndim = parseScalar<int>(line);
        }
        else if(arg == "bravaislattice"){ 
            systemInfo.bravaisLattice = parseVectors(content);
        }
        else if (arg == "motif")
        {
            systemInfo.motif = parseMotif(content);
        }
        else if (arg == "norbitals"){
            systemInfo.norbitals = parseOrbitals(content);
        }
        else if (arg == "bravaisvectors") {
            systemInfo.bravaisVectors = parseVectors(content);
        }
        else if (arg == "hamiltonian")
        {
            systemInfo.hamiltonian = parseMatrices(content);
        }
        else if (arg == "overlap") {
            systemInfo.overlap = parseMatrices(content);
        }
        else
        {
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
}

arma::mat SystemConfiguration::parseVectors(std::vector<std::string>& vectors){
    arma::mat bravaisLattice = arma::zeros(vectors.size(), 3);
    for (int i = 0; i < vectors.size(); i++){
        std::string line = vectors[i];
        std::vector<double> latticeVector = parseLine<double>(line);
        if (latticeVector.size() != 3){ 
            throw std::logic_error("Bravais vectors must have three components"); }

        bravaisLattice.row(i) = arma::rowvec(latticeVector);
    }

    return bravaisLattice;
}

arma::mat SystemConfiguration::parseMotif(std::vector<std::string>& content){
    arma::mat motif = arma::zeros(content.size(), 4);
    double x, y, z;
    std::string species;
    int index = 0;
    std::vector<double> atom;

    for (int i = 0; i < content.size(); i++){
        std::string line = content[i];
        std::istringstream iss(line);
        if ((iss >> x >> y >> z).fail()){
            throw std::invalid_argument("Motif must be of shape (x,y,z,'species')");
        }
        if ((iss >> species).fail()){
            if (systemInfo.atomToIndex.size() == 0 && index == 0){
                species = "Default";
            }
            else { throw std::logic_error("Error: Must specify all species in motif"); }
        }
        auto it = systemInfo.atomToIndex.find(species);
        if (it == systemInfo.atomToIndex.end()){
            systemInfo.atomToIndex[species] = index;
        }
        atom = { x, y, z, (float)index };

        motif.row(i) = arma::rowvec(atom);
    }

    return motif;
}

arma::urowvec SystemConfiguration::parseOrbitals(std::vector<std::string>& content) {
    if (content.size() != 1) {
        throw std::invalid_argument("Error: Orbital information must be one line only");
    }
    std::vector<int> orbitalVec = parseLine<int>(content[0]);
    arma::urowvec orbitals = arma::zeros<arma::urowvec>(orbitalVec.size());
    for (int i = 0; i < orbitalVec.size(); i++){
        orbitals(i) = orbitalVec[i];
    }

    return orbitals;
}

arma::cx_cube SystemConfiguration::parseMatrices(std::vector<std::string>& content) {
    std::vector<arma::cx_mat> matrixVector;
    double re, im;
    std::string sign, imagNumber, vartype = "real";
    int ndim = 0;

    // Scan first matrix to get information
    // ndim, triangular, real or complex values
    std::string line = content[0];
    if (line.find('i') != std::string::npos || line.find('j') != std::string::npos) {
        vartype = "complex";
    }
    for (auto i = content.begin(); i != content.end(); i++) {
        line = *i;
        if (line.find('&') != std::string::npos) {
            break;
        }
        else {
            ndim++;
        }
    }

    arma::cx_mat matrix = arma::zeros<arma::cx_mat>(ndim, ndim);
    arma::cx_rowvec row = arma::zeros<arma::cx_rowvec>(ndim);
    std::string strValue;

    // Start appending matrices
    for (int i = 0; i < content.size(); i++) {
        line = content[i];
        if (line.find('&') != std::string::npos) {
            matrixVector.push_back(matrix);
            matrix.zeros();
        }
        else{
            std::istringstream iss(line);
            int j = 0;
            while (iss >> strValue) {
                std::istringstream valuestream(strValue);
                if (vartype == "real") {
                    valuestream >> re;
                    im = 0.0;
                }
                else {
                    valuestream >> re >> sign >> im >> imagNumber;
                }
                std::complex<double> value(re, im);
                row(j) = value;
                j++;
            }
            matrix.row(i) = arma::cx_rowvec(row);
        }
    }
    if (!matrix.is_zero()) {
        matrixVector.push_back(matrix);
    }

    arma::cx_cube matrices(ndim, ndim, matrixVector.size());
    for (int i = 0; i < matrixVector.size(); i++) {
        matrices.slice(i) = matrixVector[i];
    }

    return matrices;
}

void SystemConfiguration::checkContentCoherence() {
    if (systemInfo.ndim != systemInfo.bravaisLattice.n_rows) {
        throw std::invalid_argument("Error: Dimensions must match number of Bravais basis");
    }
    if (systemInfo.hamiltonian.n_slices != systemInfo.bravaisVectors.n_rows) {
        throw std::invalid_argument("Error: Number of H matrices must match number of Bravais vectors");
    }
    if (!systemInfo.overlap.is_empty() && (systemInfo.overlap.n_slices != systemInfo.hamiltonian.n_slices)) {
        throw std::invalid_argument("Error: Number of overlap matrices must match number of H matrices");
    }
    if (systemInfo.norbitals.size() != systemInfo.atomToIndex.size()) {
        throw std::invalid_argument("Error: Number of different species must match be consistent in motif and orbitals");
    }
}