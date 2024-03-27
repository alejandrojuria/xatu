#include <complex>
#include <algorithm>
#include "xatu/SystemConfiguration.hpp"

namespace xatu {

/**
 * Default constructor.
 * @details Note that it can be called, but it will throw an error upon execution
 * of ConfigurationBase default construtor.
 */
SystemConfiguration::SystemConfiguration(){};

/**
 * File constructor. 
 * @details SystemConfiguration should always be initialized with this constructor.
 * Upon call, the configuration file is fully parsed.
 * @param filename Name of file containing the information about the system of interest.
 */
SystemConfiguration::SystemConfiguration(std::string filename) : ConfigurationBase{filename}{
    expectedArguments = {"dimension", "bravaislattice", "motif", "norbitals", "filling", "bravaisvectors", "hamiltonian"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

/**
 * Method to parse the content from the configuration file.
 * @details First, all the argument and the semi-structured contents are extracted.
 * From this, the content of each argument is parsed according to its expected shape,
 * and is finally stored in the configuration struct.
 */
void SystemConfiguration::parseContent(){
    extractArguments();
    extractRawContent();

    if (contents.empty()){
        throw std::logic_error("File contents must be extracted first");
    }
    for(const auto& arg : foundArguments){
        auto content = contents[arg];

        if (arg == "dimension"){
            if (content.size() != 1){
                throw std::logic_error("Expected only one line for 'Dimension'");
            }
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
            systemInfo.filling = parseScalar<double>(content[0]);
            double intpart;
            if (std::modf(systemInfo.filling, &intpart) != 0){
                throw std::invalid_argument("Filling must be a positive integer.");
            };
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

/**
 * Method to parse several three-dimensional vectors.
 * @details This method is intented to be used with the Bravais vectors. 
 * @param vectors Array with the strings encoding the vectors to be parsed.
 * @return Matrix with the parsed vectors as rows.
 */
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

/**
 * Method to parse the motif of the system.
 * @details This method expects to parse four-dimensional arrays corresponding to the
 * atomic coordinates plus the chemical species.
 * @param content Array storing the string with the motif information.
 * @return Motif matrix. 
 */
arma::mat SystemConfiguration::parseMotif(std::vector<std::string>& content){
    arma::mat motif = arma::zeros(content.size(), 4);
    double x, y, z, species;
    int index = 0;
    std::vector<double> atom;

    for (int i = 0; i < content.size(); i++){
        std::string line = content[i];
        std::istringstream iss(line);
        if ((iss >> x >> y >> z >> species).fail()){
            throw std::invalid_argument("Motif must be of shape (x,y,z,'species')");
        }

        atom = { x, y, z, species };
        motif.row(i) = arma::rowvec(atom);
    }

    return motif;
}

/**
 * Method to parse the orbitals.
 * @details Based upon parseLine template method.
 * @param content Array with the string to be parsed.
 * @return Vector with the orbitals per chemical species. 
 */
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

/**
 * Method to parse matrices in the configuration file.
 * @details This method can be used for both the Fock and the overlap matrices.
 * It expects dense or complete matrices, which can also be real or complex. 
 * @param content Array with the matrices to be parsed.
 * @return Cube storing the matrices.
 */
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
                }
                else {
                    valuestream >> re;
                    iss >> strValue;
                    std::istringstream valuestream(strValue);
                    valuestream >> im;
                }
                std::complex<double> value(re, im);
                row(j) = value;
                j++;
            }
            matrix.row(row_index) = arma::cx_rowvec(row);
            row_index++;
        }
    }

    arma::cx_cube matrices(ndim, ndim, matrixVector.size());
    for (int i = 0; i < matrixVector.size(); i++) {
        matrices.slice(i) = matrixVector[i];
    }

    return matrices;
}

/**
 * Method to check whether the parsed contents of the configuration file make sense.
 * @details This method makes basic checks such as making sure that the dimensionality of
 * the system if consistent, or that the number of Fock matrices matches that of Bravais vectors provided. 
 */
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

    std::vector<int> species_vec;
    for(unsigned int atomIndex = 0; atomIndex < systemInfo.motif.n_rows; atomIndex++){
        int species = systemInfo.motif.row(atomIndex)(3);
        if ( std::find(species_vec.begin(), species_vec.end(), species) == species_vec.end() ){
            species_vec.push_back(species);
        }
    }
    int nspecies = species_vec.size();

    if (systemInfo.norbitals.size() != nspecies) {
        throw std::invalid_argument("Error: Number of different species must match be consistent in motif and orbitals");
    }
}

/**
 * Auxiliary method to print the contents of the configuration struct.
 * @details Useful for debugging. 
 */
void SystemConfiguration::printConfiguration(std::ostream& stream) const {

    stream << "Dimension: " << systemInfo.ndim << "\n" << std::endl;

    stream << "Bravais lattice: " << std::endl;
    stream << systemInfo.bravaisLattice << "\n" << std::endl;

    stream << "Motif: " << std::endl;
    stream << systemInfo.motif << "\n" << std::endl;

    stream << "Orbitals: " << std::endl;
    stream << systemInfo.norbitals << "\n" << std::endl;

    stream << "Filling: " << std::endl;
    stream << systemInfo.filling << "\n" << std::endl;

    stream << "Bravais vectors (connected unit cells): " << std::endl;
    stream << systemInfo.bravaisVectors << "\n" << std::endl;

    stream << "Hamiltonian matrices: " << std::endl;
    stream << systemInfo.hamiltonian << "\n" << std::endl;

    if(!systemInfo.overlap.is_empty()){
        stream << "Overlap matrices: " << std::endl;
        stream << systemInfo.overlap << "\n" << std::endl;
    }
}

/**
 * Overload of the stream operator to print the contents of the configuration struct. 
 */
std::ostream& operator<<(std::ostream& stream, const SystemConfiguration& config){
    config.printConfiguration(stream);
}


}