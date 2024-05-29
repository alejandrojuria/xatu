#include "xatu/ConfigurationExciton.hpp"

namespace xatu {

/**
 * File constructor for ConfigurationExciton. It extracts the relevant information from the exciton file.
 * @param exciton_file Name of file with the exciton configuration.
 */
ConfigurationExciton::ConfigurationExciton(const std::string& exciton_file) : ConfigurationBase(exciton_file){
    this->expectedArguments = {"ncells", "dielectric"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

/**
 * Method to parse the exciton configuration from its file.
 * @details This method extracts all information from the configuration file and
 * stores it with the adequate format in the information struct. 
 * @return void.
 */
void ConfigurationExciton::parseContent(){
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
            excitonInfo.shift = arma::colvec(shift);
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
            excitonInfo.Q = arma::colvec(Q);
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
            std::string str = parseWord(content[0]);
            if ((str != "true") && (str != "false")){
                throw std::invalid_argument("Exchange option must be set to 'true' or 'false'.");
            }
            if (str == "true"){
                excitonInfo.exchange = true;
            }
        }
        else if(arg == "exchange.potential"){
            std::string str = parseWord(content[0]);
            std::cout << arg << " " << str << std::endl;
            excitonInfo.exchangePotential = str;
        }
        else if(arg == "potential"){
            std::string str = parseWord(content[0]);
            excitonInfo.potential = str;
        }
        else if(arg == "scissor"){
            excitonInfo.scissor = parseScalar<double>(content[0]);
        }
        else if(arg == "regularization"){
            excitonInfo.regularization = parseScalar<double>(content[0]);
        }
        else{    
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
}

/**
 * Method to check whether the information extracted from the configuration file is
 * consistent and well-defined. 
 * @return void.
 */
void ConfigurationExciton::checkContentCoherence(){
    if(excitonInfo.Q.n_elem != 3){
        throw std::logic_error("'Q' must be a 3d vector");
    };
    if(excitonInfo.ncell <= 0){
        throw std::logic_error("'ncell' must be a positive number");
    };
    if(excitonInfo.bands.empty() && excitonInfo.nbands == 0){
        throw std::logic_error("'bands' must be specified");
    };
    if(excitonInfo.eps.empty()){
        throw std::logic_error("'dielectric' must be specified");
    };
    if(excitonInfo.nbands == 0 && excitonInfo.bands.empty()){
        throw std::invalid_argument("Must specify 'nbands' or 'bandlist' parameters");
    };

    bool potentialFound = false;
    bool exchangePotentialFound = false;
    for (auto potential : supportedPotentials){
        if(excitonInfo.potential == potential){
            potentialFound = true;
        }
        if(excitonInfo.exchange && excitonInfo.exchangePotential == potential){
            exchangePotentialFound = true;
        }
    }
    if (!potentialFound){
        throw std::invalid_argument("Specified 'potential' not supported. Use 'keldysh' or 'coulomb'");
    }
    if (excitonInfo.exchange && !exchangePotentialFound){
        throw std::invalid_argument("Specified 'exchange.potential' not supported. Use 'keldysh' or 'coulomb'");
    }
    if (excitonInfo.mode != "realspace" && excitonInfo.mode != "reciprocalspace"){
        throw std::invalid_argument("Invalid mode. Use 'realspace' or 'reciprocalspace'");
    }
}

}