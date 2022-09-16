#include <armadillo>
#include "CrystalDFTConfiguration.hpp"
#include <fstream>

CrystalDFTConfiguration::CrystalDFTConfiguration(std::string file) : ConfigurationBase(file) {
    parseContent();
}

void CrystalDFTConfiguration::parseContent(){
    // Parse Crystal output file

    std::string line;
    while(std::getline(m_file, line)){
        // Bravais lattice
        if (line.find("DIRECT LATTICE") != std::string::npos){
            parseBravaisLattice();
        }

        // Motif
        else if(line.find("N. OF ATOMS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ATOMS PER CELL");
            int strsize = strlen("N. OF ATOMS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> natoms;
        }
        else if (line.find("ATOM N.AT.") != std::string::npos){
            parseAtoms();
        }
    }
}

void CrystalDFTConfiguration::parseBravaisLattice(){
    std::string line;
    std::vector<std::string> vectors;
    for(int i = 0; i < 3; i++){
        std::getline(m_file, line);
        vectors.push_back(line);     
    }
    printVector(vectors);
    systemInfo.bravaisLattice = parseVectors(vectors);
}

void CrystalDFTConfiguration::parseAtoms(){
    std::string line;
}