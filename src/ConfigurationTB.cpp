#include "xatu/ConfigurationTB.hpp"

namespace xatu {

/**
 * File constructor. 
 * @details ConfigurationTB should always be initialized with this constructor.
 * Upon call, the configuration file is fully parsed.
 * @param filename Name of file containing the information about the system of interest.
 */
ConfigurationTB::ConfigurationTB(const std::string& filename) : ConfigurationBase{filename}{
    expectedArguments = {"dimension", "bravaislattice", "motif", "norbitals", "filling", "bravaisvectors", "hamiltonian"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

/**
 * Method to parse the content from the configuration file.
 * @details First, all the argument and the semi-structured contents are extracted.
 * From this, the content of each argument is parsed according to its expected shape,
 * and is finally stored in the corresponding attributes of ConfigurationSystem.
 */
void ConfigurationTB::parseContent(){
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
            this->ndim_ = parseScalar<int>(content[0]);
        }
        else if(arg == "bravaislattice"){ 
            this->Rbasis_ = parseVectors(content);
        }
        else if (arg == "motif")
        {
            this->motif_ = parseMotif(content);
        }
        else if (arg == "norbitals"){
            this->orbitalsPerSpecies_ = parseOrbitals(content);
        }
        else if (arg == "filling"){
            if(content.size() != 1){
                throw std::logic_error("Expected only one line in 'filling' field");
            }
            this->filling_ = parseScalar<double>(content[0]);
            double intpart;
            if (std::modf(filling, &intpart) != 0){
                throw std::invalid_argument("Filling must be a positive integer.");
            };
        }
        else if (arg == "bravaisvectors") {
            this->Rlist_ = parseVectors(content);
        }
        else if (arg == "hamiltonian")
        {
            this->hamiltonianMatrices_ = parseMatrices(content);
            this->ncells_ = hamiltonianMatrices.n_slices;
        }
        else if (arg == "overlap") {
            this->overlapMatrices_ = parseMatrices(content);
        }
        else
        {
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
}

/**
 * Method to parse the motif of the system.
 * @details This method expects to parse four-dimensional arrays corresponding to the
 * atomic coordinates plus the chemical species.
 * @param content Array storing the string with the motif information.
 * @return arma::mat Motif matrix, with each column representing an atom in the unit cell. 
 */
arma::mat ConfigurationTB::parseMotif(std::vector<std::string>& content){
    arma::mat motif = arma::zeros(4, content.size());
    double x, y, z, species;
    std::vector<double> atom;

    for (uint i = 0; i < content.size(); i++){
        std::string line = content[i];
        std::istringstream iss(line);
        if ((iss >> x >> y >> z >> species).fail()){
            throw std::invalid_argument("Motif must be of shape (x,y,z,'species')");
        }

        atom = { x, y, z, species };
        motif.col(i) = arma::colvec(atom);
    }

    return motif;
}

/**
 * Method to parse the orbitals.
 * @details Based upon parseLine template method.
 * @param content Array with the string to be parsed.
 * @return Vector with the orbitals per chemical species. 
 */
arma::urowvec ConfigurationTB::parseOrbitals(std::vector<std::string>& content) {
    if (content.size() != 1) {
        throw std::invalid_argument("Error: Orbital information must be one line only");
    }
    std::vector<int> orbitalVec = parseLine<int>(content[0]);
    arma::urowvec orbitalsPerSpecies = arma::zeros<arma::urowvec>(orbitalVec.size());
    for (uint i = 0; i < orbitalVec.size(); i++){
        orbitalsPerSpecies(i) = static_cast<uint>(orbitalVec[i]);
    }

    return orbitalsPerSpecies;
}

/**
 * Method to parse matrices in the configuration file.
 * @details This method can be used for both the Fock and the overlap matrices.
 * It expects dense or complete matrices, which can also be real or complex. 
 * @param content Array with the matrices to be parsed.
 * @return Cube storing the matrices.
 */
arma::cx_cube ConfigurationTB::parseMatrices(std::vector<std::string>& content) {
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
    int row_index = 0;

    // Start appending matrices
    for (uint i = 0; i < content.size(); i++) {
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
    for (uint i = 0; i < matrixVector.size(); i++) {
        matrices.slice(i) = matrixVector[i];
    }

    return matrices;
}

/**
 * Method to check whether the parsed contents of the configuration file make sense.
 * @details This method makes basic checks such as making sure that the dimensionality of
 * the system if consistent, or that the number of Fock matrices matches that of Bravais vectors provided. 
 * @return void.
 */
void ConfigurationTB::checkContentCoherence() {
    if(static_cast<uint>(ndim) != Rbasis.n_cols) {
        throw std::invalid_argument("Error: Dimensions must match number of Bravais basis");
    }
    if (hamiltonianMatrices.n_slices != Rlist.n_cols) {
        throw std::invalid_argument("Error: Number of H matrices must match number of Bravais vectors");
    }
    if (!overlapMatrices.is_empty() && (overlapMatrices.n_slices != hamiltonianMatrices.n_slices)) {
        throw std::invalid_argument("Error: Number of overlap matrices must match number of H matrices");
    }

    std::vector<int> species_vec;
    for(unsigned int atomIndex = 0; atomIndex < motif.n_cols; atomIndex++){
        int species = motif.col(atomIndex)(3);
        if ( std::find(species_vec.begin(), species_vec.end(), species) == species_vec.end() ){
            species_vec.push_back(species);
        }
    }
    int nspecies = species_vec.size();

    if (orbitalsPerSpecies.size() != static_cast<uint>(nspecies)) {
        throw std::invalid_argument("Error: Number of different species must match be consistent in motif and orbitals");
    }
}

}