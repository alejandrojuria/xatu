#include <complex>
#include "SystemConfiguration.hpp"

SystemConfiguration::SystemConfiguration(){
    throw std::invalid_argument("Error: SystemConfiguration must be invoked with one argument (filename)");
};

SystemConfiguration::SystemConfiguration(std::string filename) : ConfigurationBase(filename){
    expectedArguments = {"dimension", "bravaislattice", "motif", "norbitals", "filling", "bravaisvectors", "hamiltonian"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

void SystemConfiguration::parseContent(){
    extractArguments();
    extractRawContent();

    if (contents.empty()){
        throw std::logic_error("File contents must be extracted first");
    }
    for(const auto& arg : foundArguments){
        auto content = contents[arg];

        if (arg == "dimension"){
            if(content.size() != 1){throw std::logic_error("Expected only one line for 'Dimension'");}
            systemInfo.ndim = parseScalar<int>(content[0]);
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
        else if (arg == "filling"){
            if(content.size() != 1){
                throw std::logic_error("Expected only one line in 'filling' field");
            }
            systemInfo.filling = parseFraction(content[0]);
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
            species = "Default";
        }
        auto it = systemInfo.atomToIndex.find(species);
        if (it == systemInfo.atomToIndex.end()){
            systemInfo.atomToIndex[species] = index;
            index++;
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
    char placeholder; // To read uninteresting characters
    std::string sign, imagNumber, vartype = "real";
    bool isTriangular = false;
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
    int row_index = 0;

    // Start appending matrices
    for (int i = 0; i < content.size(); i++) {
        line = content[i];
        if (line.find('&') != std::string::npos) {
            matrixVector.push_back(matrix);
            matrix.zeros();
            row_index = 0;
        }
        else{
            std::istringstream iss(line);
            int j = 0;
            row.zeros();
            while (iss >> strValue) {
                std::istringstream valuestream(strValue);
                if (vartype == "real") {
                    valuestream >> re;
                    im = 0.0;
                    std::cout << "this" << std::endl;
                }
                else {
                    valuestream >> placeholder >> re >> im >> placeholder;
                }
                std::complex<double> value(re, im);
                row(j) = value;
                j++;
            }
            matrix.row(row_index) = arma::cx_rowvec(row);
            row_index++;
        }
    }
    /*if (!matrix.is_zero()) {
        // Check if matrix is triangular
        if (matrix.row(0)(matrix.n_cols - 1) != std::conj(matrix.row(matrix.n_rows - 1)(0))){
            isTriangular = true;

        }
        matrixVector.push_back(matrix);
    }*/

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