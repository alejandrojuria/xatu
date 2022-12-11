#include "xatu/ExcitonConfiguration.hpp"
#include <fstream>
#include <iostream>

namespace xatu {

ExcitonConfiguration::ExcitonConfiguration(){
    throw std::invalid_argument("Error: ExcitonConfiguration must be invoked with one argument (filename)");
};

ExcitonConfiguration::ExcitonConfiguration(std::string filename) : ConfigurationBase(filename){
    this->expectedArguments = {"ncells", "bands", "dielectric"};
    parseContent();
    checkArguments();
    checkContentCoherence();
}

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

        if(arg == "ncells"){
            excitonInfo.ncell = parseScalar<int>(content[0]);
        }
        else if(arg == "submesh"){
            excitonInfo.nbands = parseScalar<int>(content[0]);
        }
        else if(arg == "shiftmesh"){
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
        }
        else if(arg == "exchange"){
            excitonInfo.exchange = true;
        }
        else{    
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
};


void ExcitonConfiguration::checkContentCoherence(){
    if(excitonInfo.Q.n_elem != 3){
        throw std::logic_error("Q must be a 3d vector");
    };
    if(excitonInfo.ncell < 0){
        throw std::logic_error("ncell must be a positive number");
    };
    if(excitonInfo.bands.empty() && excitonInfo.nbands == 0){
        throw std::logic_error("bands must be specified");
    };
    if(excitonInfo.eps.empty()){
        throw std::logic_error("eps must be specified");
    };
};

}