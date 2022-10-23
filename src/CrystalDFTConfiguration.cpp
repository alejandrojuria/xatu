#include <armadillo>
#include "CrystalDFTConfiguration.hpp"
#include "SystemConfiguration.hpp"
#include <fstream>
#include <algorithm>

CrystalDFTConfiguration::CrystalDFTConfiguration(std::string file, int ncells) : ConfigurationBase(file) {
    parseContent(ncells);
    mapContent();
}

void CrystalDFTConfiguration::parseContent(int ncells, double threshold){
    // Parse Crystal output file

    std::string line;
    while(std::getline(m_file, line)){
        // Bravais lattice
        if (line.find("DIRECT LATTICE VECTOR COMPONENTS") != std::string::npos){
            parseBravaisLattice(threshold);
            arma::cout << "Bravais lattice: " << arma::endl;
            arma::cout << bravaisLattice << arma::endl;
        }

        // Motif
        else if(line.find("N. OF ATOMS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ATOMS PER CELL");
            int strsize = strlen("N. OF ATOMS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> natoms;
            arma::cout << "N. atoms per cell: " << natoms << arma::endl;
        }

        // N. shells (total)
        else if(line.find("NUMBER OF SHELLS") != std::string::npos){
            int pos = line.find("NUMBER OF SHELLS");
            int strsize = strlen("NUMBER OF SHELLS");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> nshells;
            arma::cout << "N. shells: " << nshells << arma::endl;
        }

        // N. orbitals (total)
        else if(line.find("NUMBER OF AO") != std::string::npos){
            int pos = line.find("NUMBER OF AO");
            int strsize = strlen("NUMBER OF AO");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> norbitals;
            arma::cout << "N. orbitals: " << norbitals << arma::endl;
        }

        // N. electrons
        else if(line.find("N. OF ELECTRONS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ELECTRONS PER CELL");
            int strsize = strlen("N. OF ELECTRONS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> valenceElectrons;

            arma::cout << "N. electrons per cell: " << valenceElectrons << arma::endl;
        }

        // N. electrons
        else if(line.find("CORE ELECTRONS PER CELL") != std::string::npos){
            int pos = line.find("CORE ELECTRONS PER CELL");
            int strsize = strlen("CORE ELECTRONS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> coreElectrons;

            arma::cout << "N. core electrons per cell: " << coreElectrons << arma::endl;
        }

        // 
        else if ((line.find("ATOM") != std::string::npos) && line.find("SHELL") != std::string::npos){
            if(natoms == 0){
                throw std::logic_error("Must parse first number of atoms");
            }
            parseAtoms();
            arma::cout << "Motif: " << arma::endl;
            arma::cout << motif << arma::endl;
            printVector(shellsPerAtom);
            arma::cout << "N. species: " << nspecies << arma::endl;
        }

        // Parse atomic basis info
        else if(line.find("LOCAL ATOMIC FUNCTIONS BASIS SET") != std::string::npos){
            parseAtomicBasis();
            printVector(orbitalsPerAtom);
            for(auto const& [key, val]: gaussianCoefficients){
                arma::cout << val << arma::endl;
            }
            for(auto const& [key, val]: shellTypesPerAtom){
                for(int i = 0; i < val.size(); i++){
                    arma::cout << key << arma::endl;
                    arma::cout << val[i] << arma::endl;
                }
            }
        }

        else if(line.find("OVERLAP MATRIX") != std::string::npos){
            int pos = line.find("OVERLAP MATRIX - CELL N.");
            int strsize = strlen("OVERLAP MATRIX - CELL N.");
            line = line.substr(pos + strsize, line.length());

            int cellIndex, x, y, z;
            std::string parenthesis;

            std::istringstream iss(line);
            iss >> cellIndex >> parenthesis >> x >> y >> z;
            std::vector coefCombinations = {x, y, z};

            if(cellIndex <= ncells){
                arma::rowvec cell = arma::rowvec(3);
                for (int i = 0; i < ndim; i++){
                    cell += bravaisLattice.row(i)*coefCombinations[i];
                }
                this->bravaisVectors = arma::join_vert(bravaisVectors, cell);

                arma::cx_mat overlapMatrix = parseMatrix();

                this->overlapMatrices = arma::join_slices(this->overlapMatrices, overlapMatrix);
            }
        }

        else if(line.find("FOCK MATRIX") != std::string::npos){
            int pos = line.find("FOCK MATRIX - CELL N.");
            int strsize = strlen("FOCK MATRIX - CELL N.");
            line = line.substr(pos + strsize, line.length());

            int cellIndex;
            std::string parenthesis;
            std::istringstream iss(line);
            iss >> cellIndex;

            if(cellIndex <= ncells){
                arma::cx_mat fockMatrix = parseMatrix();

                this->fockMatrices = arma::join_slices(this->fockMatrices, fockMatrix);
            }
        }
    }    
}

void CrystalDFTConfiguration::parseBravaisLattice(double threshold){
    std::string line;
    std::vector<std::string> vectors;
    for(int i = 0; i < 3; i++){
        std::getline(m_file, line);
        vectors.push_back(line);     
    }
    this->bravaisLattice = parseVectors(vectors);
    extractDimension(threshold);
}

void CrystalDFTConfiguration::extractDimension(double threshold){
    for(unsigned int i = 0; i < bravaisLattice.n_rows; i++){
        double norm = arma::norm(bravaisLattice.row(i));
        if (norm > threshold){
            bravaisLattice.shed_row(i);
        }
    }

    ndim = bravaisLattice.n_rows;
}

void CrystalDFTConfiguration::parseAtoms(){
    std::string line;
    int index, natom, nshells, nspecies = 0;
    double x, y, z;
    std::string chemical_species;
    std::vector<int> shellsPerAtom;
    std::vector<double> atom;
    std::map<std::string, int> chemical_species_to_index;
    std::vector<std::string> species;
    arma::mat motif = arma::zeros(natoms, 4);

    std::getline(m_file, line); // Parse asterisks
    for(int i = 0; i < natoms; i++){
        std::getline(m_file, line);
        std::istringstream iss(line);
        iss >> index >> natom >> chemical_species >> nshells >> x >> y >> z;

        if(std::find(species.begin(), species.end(), chemical_species) == species.end()){
            species.push_back(chemical_species);
            shellsPerAtom.push_back(nshells);
            chemical_species_to_index[chemical_species] = chemical_species_to_index.size();
            nspecies++;
        }

        int index = chemical_species_to_index[chemical_species];
        atom = {x, y, z, (double)index};
        motif.row(i) = arma::rowvec(atom);
    }

    this->motif = motif;
    this->shellsPerAtom = shellsPerAtom;
    this->nspecies = nspecies;
}

void CrystalDFTConfiguration::parseAtomicBasis(){
    std::string line;
    int norbitals, totalOrbitals = 0;
    std::string shellType;
    std::vector<std::string> shellTypes;
    double exponent, sCoef, pCoef, dCoef;
    std::vector<double> coefs;
    arma::cube gaussianCoefficients;

    std::getline(m_file, line); // Parse asterisks
    std::getline(m_file, line); // Parse header
    std::getline(m_file, line); // Parse asterisks
    
    for(int atomIndex = 0; atomIndex < this->natoms; atomIndex++){
        std::getline(m_file, line); // Skip line with positions
        gaussianCoefficients = arma::zeros(4, 4, shellsPerAtom[atomIndex]);
        shellTypes.clear();
        for(int shellIndex = 0; shellIndex < shellsPerAtom[atomIndex]; shellIndex++){
            std::getline(m_file, line);
            std::istringstream iss(line);
        
            iss >> norbitals >> shellType;
            if(shellType == "-"){
                iss >> norbitals >> shellType;
            }

            shellTypes.push_back(shellType);

            for(int i = 0; i < 4; i++){
                // Suppose that there are only four coefs per shell
                std::getline(m_file, line);
                std::istringstream iss(line);
                iss >> exponent >> sCoef >> pCoef >> dCoef;
                coefs = {exponent, sCoef, pCoef, dCoef};
                
                gaussianCoefficients.slice(shellIndex).row(i) = arma::rowvec(coefs);
            }

            this->gaussianCoefficients[atomIndex] = gaussianCoefficients;
        }
        this->orbitalsPerAtom.push_back(norbitals - totalOrbitals);
        totalOrbitals += norbitals;

        this->shellTypesPerAtom[atomIndex] = shellTypes;
    }

}

arma::cx_mat CrystalDFTConfiguration::parseMatrix(){
    std::string line;
    arma::cx_mat matrix = arma::zeros<arma::cx_mat>(norbitals, norbitals);
    bool firstNonEmptyLineFound = false;
    arma::cx_rowvec matrixRow;
    double coef;
    int i, j = 0, index = 0;

    while(std::getline(m_file, line)){
        if (line.empty() && matrix.is_zero()){
            continue;
        }
        else if (!line.empty() && !firstNonEmptyLineFound){
            firstNonEmptyLineFound = true;
            continue;
        }
        else if (line.empty() && firstNonEmptyLineFound){
            return matrix;
        }
        std::istringstream iss(line);
        matrixRow = arma::zeros<arma::cx_rowvec>(norbitals);

        i = 0;
        iss >> index;
        while(iss >> coef){
            matrixRow(i) = coef;
            i++;
        }

        matrix.row(j) = matrixRow;
        j++;
    }
    return matrix;
}

void CrystalDFTConfiguration::mapContent(){
    systemInfo.bravaisLattice = bravaisLattice;
    systemInfo.motif          = motif;
    systemInfo.filling        = (valenceElectrons + coreElectrons)/2;
    systemInfo.hamiltonian    = fockMatrices;
    systemInfo.overlap        = overlapMatrices;
    systemInfo.ndim           = ndim;
    systemInfo.bravaisVectors = bravaisVectors;

    arma::urowvec norbitals = arma::zeros<arma::urowvec>(orbitalsPerAtom.size());
    for (int i = 0; i < orbitalsPerAtom.size(); i++){
        norbitals(i) = orbitalsPerAtom[i];
    }
    systemInfo.norbitals      = norbitals;
}