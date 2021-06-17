#include "SystemConfiguration.hpp"

SystemConfiguration::SystemConfiguration(){
    throw std::invalid_argument("SystemConfiguration must be called with one argument (filename)");
};

SystemConfiguration::SystemConfiguration(std::string filename) : ConfigurationBase(filename){
    expectedArguments = {"dimension", "bravaislattice", "motif", "norbitals", "hamiltonian"};
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
            systemInfo.bravaisLattice = parseBravaisLattice(content);
        }
        else if (arg == "motif")
        {
            systemInfo.motif = parseMotif(content);
        }
        else if (arg == "norbitals"){
            systemInfo.norbitals = parseOrbitals(content);
        }
        else if (arg == "hamiltonian")
        {
            systemInfo.hamiltonian = parseHamiltonian(content);
        }
        else
        {
            std::cout << "Unexpected argument: " << arg << ", skipping block..." << std::endl;
        }
    }
}

arma::mat SystemConfiguration::parseBravaisLattice(std::vector<std::string>& vectors){
    arma::mat bravaisLattice = arma::zeros(vectors.size(), 3);
    for (int i = 0; i < vectors.size(); i++){
        std::string line = vectors[i];
        std::vector<float> latticeVector = parseLine<float>(strvector);
        if (latticeVector.size() != 3){ 
            throw std::logic_error("Bravais lattice vectors must have three components"); }

        bravaisLattice.row(i) = arma::rowvec(latticeVector);
    }

    return bravaisLattice;
}

arma::mat SystemConfiguration::parseMotif(std::vector<std::string>& content){
    arma::mat motif = arma::zeros(content.size(), 4);
    float x, y, z;
    std::string species;
    int index = 0;
    stdd::vector<float> atom;

    for (int i = 0; i < content.size(); i++){
        std::string line = content[i];
        std::istringstream iss(line);
        if ((iss >> x >> y >> z).fail()){
            throw std::invalid_argument("Motif must be composed of (x,y,z;'species')");
        }
        if ((iss >> species).fail()){
            if (systemInfo.atomToIndex.size() == 0){
                systemInfo.atomToIndex["Default"] = index;
            }
        }
        else{
            auto it = systemInfo.atomToIndex.find(species);
            if (it == systemInfo.atomToIndex.end()){
                systemInfo.atomToIndex[species] = index;
            }
            if (systemInfo.atomToIndex.find("Default") == systemInfo.atomToIndex.end())
            


            }
        }


        catch(){
            throw std::logic_error("Motif must be composed of (x,y,z;'species')");
        }

        std::vector<float> latticeVector = parseLine<float>(line)

        // Fill vector if no atomic species was given
        if ((latticeVector.size() != 3) || (latticeVector.size() != 4))
        if (latticeVector.size() == 3){
            latticeVector.push_back(0);
        }
        else if (latticeVector.size() != 4){
            throw std::logic_error("Motif positions must be (x, y, z, species_index)"); }

        motif.row(i) = arma::rowvec(latticeVector);
    }

    return motif;
}