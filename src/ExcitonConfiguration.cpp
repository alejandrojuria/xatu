#include "xatu/ExcitonConfiguration.hpp"
#include <fstream>
#include <iostream>

namespace xatu {

/**
 * Default constructor.
 * @details Throws an error if called. Class must be constructed providing a filename with the configuration. 
 */
ExcitonConfiguration::ExcitonConfiguration(){
    throw std::invalid_argument("Error: ExcitonConfiguration must be invoked with one argument (filename)");
};

/**
 * File constructor. 
 * @details ExcitonConfiguration must be initialized always with this constructor.
 * Upon call, the exciton configuration file is fully parsed and the information extracted.
 * @param filename Name of file with the exciton configuration.
 */
ExcitonConfiguration::ExcitonConfiguration(std::string filename) : ConfigurationBase(filename){
    this->expectedArguments = {"ncells", "dielectric"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

/**
 * Method to parse the exciton configuration from its file.
 * @details This method extracts all information from the configuration file and
 * stores it with the adequate format in the information struct. 
 */
void ExcitonConfiguration::parseContent(){
    extractArguments();
    extractRawContent();

    if (contents.empty()){
        throw std::logic_error("File contents must be extracted first");
    }
    for(const auto& arg : foundArguments){
        auto content = contents[arg];

        if (content.size() == 0){
            continue;
        }
        else if(content.size() != 1){
            throw std::logic_error("Expected only one line per field");
        }

        if(arg == "label"){
            excitonInfo.label = standarizeLine(content[0]);
        }
        else if(arg == "ncells"){
            excitonInfo.ncell = parseScalar<int>(content[0]);
        }
        else if(arg == "submesh"){
            excitonInfo.submeshFactor = parseScalar<int>(content[0]);
        }
        else if(arg == "shift"){
            std::vector<double> shift = parseLine<double>(content[0]);
            excitonInfo.shift = arma::rowvec(shift);
        }
        else if(arg == "bands"){
            excitonInfo.nbands = parseScalar<int>(content[0]);
        }
        else if(arg == "bandlist"){
            // uint does not work with arma::urowvec, so we use typedef arma::uword
            std::vector<arma::s64> bands = parseLine<arma::s64>(content[0]);
            excitonInfo.bands = arma::ivec(bands);
        }
        else if(arg == "totalmomentum"){
            std::vector<double> Q = parseLine<double>(content[0]);
            excitonInfo.Q = arma::rowvec(Q);
        }
        else if(arg == "cutoff"){
            excitonInfo.cutoff = parseScalar<double>(content[0]);
        }
        else if(arg == "dielectric"){
            std::vector<double> eps = parseLine<double>(content[0]);
            excitonInfo.eps = arma::vec(eps);
        }
        else if(arg == "reciprocal"){
            excitonInfo.mode = "reciprocalspace";
            excitonInfo.nReciprocalVectors = parseScalar<int>(content[0]);
        }
        else if(arg == "exchange"){
            std::string str = content[0];
            str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
            std::transform(str.begin(), str.end(), str.begin(), ::tolower);
            if ((str != "true") && (str != "false")){
                throw std::invalid_argument("Exchange option must be set to 'true' or 'false'.");
            }
            if (str == "true"){
                excitonInfo.exchange = true;
            }
        }
        else if(arg == "scissor"){
            excitonInfo.scissor = parseScalar<double>(content[0]);
        }
        else{    
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
};

/**
 * Method to check whether the information extracted from the configuration file is
 * consistent and well-defined. 
 */
void ExcitonConfiguration::checkContentCoherence(){
    if(excitonInfo.Q.n_elem != 3){
        throw std::logic_error("Q must be a 3d vector");
    };
    if(excitonInfo.ncell <= 0){
        throw std::logic_error("ncell must be a positive number");
    };
    if(excitonInfo.bands.empty() && excitonInfo.nbands == 0){
        throw std::logic_error("bands must be specified");
    };
    if(excitonInfo.eps.empty()){
        throw std::logic_error("eps must be specified");
    };
    if(excitonInfo.nbands == 0 && excitonInfo.bands.empty()){
        throw std::invalid_argument("Must specify 'nbands' or 'bandlist' parameters");
    }
};

}