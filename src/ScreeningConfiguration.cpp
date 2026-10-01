#include <fstream>
#include <iostream>
#include "xatu/ScreeningConfiguration.hpp"

namespace xatu {

/**
 * Default constructor.
 * @details Throws an error if called. Class must be constructed providing a filename with the configuration. 
 */
ScreeningConfiguration::ScreeningConfiguration(){
    throw std::invalid_argument("Error: ScreeningConfiguration must be invoked with one argument (filename)");
};

/**
 * File constructor. 
 * @details ScreeningConfiguration must be initialized always with this constructor.
 * Upon call, the exciton configuration file is fully parsed and the information extracted.
 * @param filename Name of file with the exciton configuration.
 */
ScreeningConfiguration::ScreeningConfiguration(std::string filename) : ConfigurationBase(filename){
    this->expectedArguments = {"function","valence.bands","conduction.bands","ncell_aux","spin","gcutoff"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

/**
 * Method to parse the exciton configuration from its file.
 * @details This method extracts all information from the configuration file and
 * stores it with the adequate format in the information struct. 
 */
void ScreeningConfiguration::parseContent(){
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
        else if(arg == "valence.bands"){
            screeningInfo.nvbands = parseScalar<int>(content[0]);
        }
        else if(arg == "conduction.bands"){
            screeningInfo.ncbands = parseScalar<int>(content[0]);
        }
        else if(arg == "vectors"){
            std::vector<int> gs = parseLine<int>(content[0]);
            screeningInfo.Gs(0) = gs[0];
            screeningInfo.Gs(1) = gs[1];
        }
        else if(arg == "momentum"){
            std::vector<double> momentum = parseLine<double>(content[0]);
            screeningInfo.q = arma::rowvec(momentum);
        }
        else if(arg == "function"){
            screeningInfo.function = parseWord(content[0]);
        }
        else if(arg == "gcutoff"){
            screeningInfo.Gcutoff = parseScalar<double>(content[0]);
        }
        else if(arg == "spin"){
            std::string spin_string = parseWord(content[0]);
            if (spin_string == "false") {
                screeningInfo.spin = false;
            } else if (spin_string == "true") {
                screeningInfo.spin = true;
            } else {
                std::cout << "Option for spin not recognized, skipping block..." << std::endl;
            }
        } else if(arg == "ncell_aux"){
            screeningInfo.ncell_aux = parseScalar<int>(content[0]);        
        } else if(arg == "thickness"){
            screeningInfo.d = parseScalar<double>(content[0]);        
        } else if(arg == "isotropic"){
            std::string isotropic_string = parseWord(content[0]);
            
            if (isotropic_string == "false") {
                screeningInfo.isotropic = false;
            } else if (isotropic_string == "true") {
                screeningInfo.isotropic = true;
            } else {
                std::cout << "Option for isotropic not recognized, skipping block..." << std::endl;
            }
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
void ScreeningConfiguration::checkContentCoherence(){
    if(screeningInfo.nvbands <= 0){
        throw std::logic_error("'nvbands' must be a positive number");
    };
    if(screeningInfo.ncbands <= 0){
        throw std::logic_error("'ncbands' must be a positive number");
    }
    if(screeningInfo.function != "dielectric" && screeningInfo.function != "polarizability" && screeningInfo.function != "inversedielectric" && screeningInfo.function != "exciton"){
        throw std::logic_error("'function' must be 'dielectric', 'polarizability', 'inversedielectric' or 'exciton'");
    }
    if(screeningInfo.Gs(0) < 0 || screeningInfo.Gs(1) < 0){
        throw std::invalid_argument("The indeces of the vectors must be zero or a positive integer.");
    }
    if(screeningInfo.ncell_aux < 1){
        throw std::invalid_argument("The number of points in the coarser BZ mesh in each direction, ncell_aux, has to be positive.");
    }
    if(screeningInfo.d < 0){
        throw std::invalid_argument("The thickness of the material must not be negative.");
    }
};

}