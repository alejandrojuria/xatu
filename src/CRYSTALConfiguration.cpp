#include "xatu/CRYSTALConfiguration.hpp"

#define SOC_STRING "to_be_defined_for_crystal23"
#define MAGNETIC_STRING "UNRESTRICTED OPEN SHELL"

namespace xatu {

/**
 * File constructor for CRYSTALConfiguration. It extracts the relevant information
 * from the file and stores it in an adequate format.
 * @details This class is intended to be used with .outp files from the CRYSTAL code.
 * Since orbitals in CRYSTAL extend over several unit cells, the Fock matrices that define the
 * Hamiltonian also cover several unit cells. Therefore, one can specify how many unit cells to take
 * for the actual exciton calculation. 
 * @param file Name of the .outp from the CRYSTAL SCF calculation.
 * @param ncells Number of unit cells for which the Hamiltonian and overlap matrices are read from file.
 */
CRYSTALConfiguration::CRYSTALConfiguration(std::string file, int ncells) : ConfigurationBase(file) {
    parseContent(ncells);
    mapContent();
}

/**
 * Method to extract all the content from the file.
 * @details The lattice dimensionality is determined according to CRYSTAL's criterion for the 
 * direct lattice vectors: in 2D, R3 == [0 0 500]; and in 1D, in addition R2 = [0 0 500]. The 
 * hypothetical case of R3 and/or R2 == [0 0 0] is also assumed to lower the dimensionality.
 * @param ncells Number of unit cells to be parsed.
 * @return void
 */
void CRYSTALConfiguration::parseContent(int ncells){
    // Parse Crystal output file

    int countr = 0;
    std::string line;
    while(std::getline(m_file, line)){
        // Bravais lattice
        if (line.find("DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)") != std::string::npos){
            parseBravaisLattice();
            countr++;
        }

        // Motif
        else if(line.find("N. OF ATOMS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ATOMS PER CELL");
            int strsize = strlen("N. OF ATOMS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> natoms;
        }

        // N. shells (total)
        else if(line.find("NUMBER OF SHELLS") != std::string::npos){
            int pos = line.find("NUMBER OF SHELLS");
            int strsize = strlen("NUMBER OF SHELLS");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> nsh;
        }

        // N. orbitals (total)
        else if(line.find("NUMBER OF AO") != std::string::npos){
            int pos = line.find("NUMBER OF AO");
            int strsize = strlen("NUMBER OF AO");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> norbitals;            
        }

        // N. electrons
        else if(line.find("N. OF ELECTRONS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ELECTRONS PER CELL");
            int strsize = strlen("N. OF ELECTRONS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> totalElectrons;
        }

        // N. core electrons
        else if(line.find("CORE ELECTRONS PER CELL") != std::string::npos){
            int pos = line.find("CORE ELECTRONS PER CELL");
            int strsize = strlen("CORE ELECTRONS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> coreElectrons;
        }

        // 
        else if ((line.find("ATOM") != std::string::npos) && line.find("SHELL") != std::string::npos){
            if(natoms == 0){
                throw std::logic_error("Must parse first number of atoms");
            }
            parseAtoms();
        }

        // Parse atomic basis info
        else if(line.find("LOCAL ATOMIC FUNCTIONS BASIS SET") != std::string::npos){
            parseAtomicBasis();
        }

        else if(line.find("OVERLAP MATRIX") != std::string::npos){
            int pos = line.find("OVERLAP MATRIX - CELL N.");
            int strsize = strlen("OVERLAP MATRIX - CELL N.");
            line = line.substr(pos + strsize, line.length());

            int cellIndex, x, y, z;
            std::string parenthesis;

            std::istringstream iss(line);
            iss >> cellIndex >> parenthesis >> x >> y >> z;
            std::vector<int> coefCombinations = {x, y, z};

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

        // Toggle SOC flag if calculation is done with spin-orbit coupling
        else if(line.find(SOC_STRING) != std::string::npos){
            this->SOC_FLAG = true;
        }

        else if(line.find(MAGNETIC_STRING) != std::string::npos){
            this->MAGNETIC_FLAG = true;
        }

        if ((line.find("BETA") != std::string::npos) && (line.find("ELECTRONS") != std::string::npos)){
            alpha_electrons = false;
        }

        if(line.find("FOCK MATRIX") != std::string::npos){
            int pos = line.find("FOCK MATRIX - CELL N.");
            int strsize = strlen("FOCK MATRIX - CELL N.");
            line = line.substr(pos + strsize, line.length());

            int cellIndex, x, y, z;
            std::string parenthesis;
            std::istringstream iss(line);
            iss >> cellIndex >> parenthesis >> x >> y >> z;

            if(cellIndex <= ncells){
                arma::cx_mat fockMatrix = parseMatrix();
                if(MAGNETIC_FLAG){
                    if(alpha_electrons){
                        this->alphaMatrices = arma::join_slices(this->alphaMatrices, fockMatrix);
                    }
                    else{
                        this->betaMatrices = arma::join_slices(this->betaMatrices, fockMatrix);
                    }
                }
                else if(SOC_FLAG){
                    // Do something here
                    continue;
                }
                else{
                    this->fockMatrices = arma::join_slices(this->fockMatrices, fockMatrix);
                }
            }
        }
    }    

    if (countr==0){
        ndim = 0;
        throw std::logic_error("Finite systems are not yet implemented");
    } 
}

/**
 * Method to parse and format the Bravais basis vectors from the file. 
 * @return void
 */
void CRYSTALConfiguration::parseBravaisLattice(){
    std::string line;
    std::vector<std::string> vectors;
    for(int i = 0; i < 3; i++){
        std::getline(m_file, line);
        vectors.push_back(line);     
    }
    this->bravaisLattice = parseVectors(vectors);
    extractDimension();
}

/**
 * Method to obtain the dimension (1D,2D,3D) of the system.
 * @return void 
 */
void CRYSTALConfiguration::extractDimension(){
    arma::rowvec R0 = arma::zeros<arma::rowvec>(3);
    arma::rowvec R2 = {0.0, 500.0, 0.0};
    arma::rowvec R3 = {0.0, 0.0, 500.0}; 
                 
    if ( approx_equal(bravaisLattice.row(2),R3,"absdiff",0.1) || 
         approx_equal(bravaisLattice.row(2),R0,"absdiff",0.1) ){
        bravaisLattice.shed_row(2);
        if ( approx_equal(bravaisLattice.row(1),R2,"absdiff",0.1) || 
             approx_equal(bravaisLattice.row(1),R0,"absdiff",0.1) ){
            bravaisLattice.shed_row(1); 
        }
    }
    ndim = bravaisLattice.n_rows;
}

/**
 * Method to extract the motif, the chemical species and the number of shells per species.
 * @return void 
 */
void CRYSTALConfiguration::parseAtoms(){
    std::string line;
    int index, atomic_number, nsh, nspecies = 0;
    double x, y, z;
    std::string chemical_species;
    std::vector<int> shellsPerSpecies;
    std::vector<double> atom;
    std::map<std::string, int> chemical_species_to_index;
    std::vector<std::string> species;
    arma::mat motif = arma::zeros(natoms, 4);

    std::getline(m_file, line); // Parse asterisks
    for(int i = 0; i < natoms; i++){
        std::getline(m_file, line);
        std::istringstream iss(line);
        iss >> index >> atomic_number >> chemical_species >> nsh >> x >> y >> z;

        if(std::find(species.begin(), species.end(), chemical_species) == species.end()){
            species.push_back(chemical_species);
            shellsPerSpecies.push_back(nsh);
            chemical_species_to_index[chemical_species] = chemical_species_to_index.size();
            nspecies++;
        }

        int index = chemical_species_to_index[chemical_species];
        atom = {x, y, z, (double)index};
        motif.row(i) = arma::rowvec(atom);
    }
    this->motif = motif;
    this->shellsPerSpecies = shellsPerSpecies;
    this->nspecies = nspecies;
}

/**
 * Method to extract the details of the basis used in the CRYSTAL calculation.
 * @details This method extracts all the orbitals per chemical species and
 * the corresponding coefficients of the gaussian expansion.
 * @return void 
 */
void CRYSTALConfiguration::parseAtomicBasis(){
    std::string line, chemical_species;
    int norbitals, natom, totalOrbitals = 0, nspecies = 0;
    std::string shellType;
    // std::vector<std::string> shellTypes;
    double exponent, sCoef, pCoef, dCoef;
    std::vector<double> coefs;
    cube_vector gaussianCoefficients;
    std::vector<std::string> species;

    std::getline(m_file, line); // Parse asterisks
    std::getline(m_file, line); // Parse header
    std::getline(m_file, line); // Parse asterisks
    
    for(int atomIndex = 0; atomIndex < this->natoms; atomIndex++){

        // First parse chemical species and add to list if not present.
        std::getline(m_file, line);
        std::istringstream iss(line);
        iss >> natom >> chemical_species; 

        if(std::find(species.begin(), species.end(), chemical_species) == species.end()){
            species.push_back(chemical_species);
        }
        else{
            totalOrbitals += orbitalsPerSpecies[motif.row(atomIndex)(3)]; // Track total num. orbitals
            continue; // Skip already parsed
        }

        gaussianCoefficients.clear();
        // shellTypes.clear();

        for(int shellIndex = 0; shellIndex < shellsPerSpecies[nspecies]; shellIndex++){

            std::getline(m_file, line);
            std::istringstream iss(line);
            
            iss >> norbitals >> shellType;
            if(shellType == "-"){
                iss >> norbitals >> shellType;
            }

            // shellTypes.push_back(shellType);

            std::vector<std::vector<double>> coefList;
            long int previousLine = m_file.tellg(); // Store beginning of next line
            while (std::getline(m_file, line)){
                std::vector<double> tokenized_line;
                split(line, tokenized_line);
                if (tokenized_line.size() != 4){
                    m_file.seekg(previousLine); // Restore line
                    break;
                }
                
                std::istringstream iss(line);
                iss >> exponent >> sCoef >> pCoef >> dCoef;
                coefs = {exponent, sCoef, pCoef, dCoef};

                coefList.push_back(coefs);

                previousLine = m_file.tellg(); // Store beginning of next line
            }
            gaussianCoefficients.push_back(coefList);

            this->gaussianCoefficients[nspecies] = gaussianCoefficients;
        }
        this->orbitalsPerSpecies.push_back(norbitals - totalOrbitals);
        totalOrbitals = norbitals;

        // this->shellTypesPerSpecies[nspecies] = shellTypes;
        nspecies++;
    }

}

/**
 * Method to parse the Fock and overlap matrices from the input file.
 * @return void
 */
arma::cx_mat CRYSTALConfiguration::parseMatrix(){
    std::string line;
    arma::cx_mat matrix = arma::zeros<arma::cx_mat>(norbitals, norbitals);
    bool firstNonEmptyLineFound = false;
    arma::cx_rowvec matrixRow;
    double coef;
    std::string colIndicesStr;
    std::vector<int> colIndices;
    int i, index, colIndex = 0, rowIndex = 0;

    while(std::getline(m_file, line)){
        if (line.empty()){
            colIndices.clear();
            std::getline(m_file, colIndicesStr); // After blank line get indices
            std::istringstream iss(colIndicesStr);
            while(iss >> index){
                colIndices.push_back(index);
            }
            std::getline(m_file, line); // Get next line
            if (line.empty()){
                continue;
            }
        }
        
        i = 0;
        std::istringstream iss(line);
        iss >> rowIndex;
        while(iss >> coef){
            colIndex = colIndices[i];
            matrix(rowIndex - 1, colIndex - 1) = coef;
            i++;
        }

        if (rowIndex == norbitals && colIndex == norbitals){
            return matrix;
        }
    }
}

/**
 * Method to write all the extracted information into a struct.
 * @return void 
 */
void CRYSTALConfiguration::mapContent(bool debug){

    systemInfo.ndim           = ndim;
    systemInfo.bravaisLattice = bravaisLattice;
    systemInfo.motif          = motif;
    systemInfo.filling        = totalElectrons/2;
    systemInfo.bravaisVectors = bravaisVectors;
    systemInfo.overlap        = overlapMatrices;

    arma::urowvec norbitals = arma::zeros<arma::urowvec>(orbitalsPerSpecies.size());
    for (int i = 0; i < orbitalsPerSpecies.size(); i++){
        norbitals(i) = orbitalsPerSpecies[i];
    }
    systemInfo.norbitals      = norbitals;

    if (SOC_FLAG) {
        systemInfo.filling   *= 2; 
        systemInfo.norbitals *= 2;
        // Pending Hamiltonian initialization; for future CRYSTAL releases
    }
    else if(MAGNETIC_FLAG){
        systemInfo.filling   *= 2;
        systemInfo.norbitals *= 2;

        arma::mat spinUpBlock = {{1, 0}, {0, 0}};
        arma::mat spinDownBlock = {{0, 0}, {0, 1}};

        arma::cx_cube newOverlapMatrices;
        for(unsigned int i = 0; i < alphaMatrices.n_slices; i++){
            arma::cx_mat totalFockMatrix = arma::kron(alphaMatrices.slice(i), spinUpBlock) + arma::kron(betaMatrices.slice(i), spinDownBlock);
            this->fockMatrices = arma::join_slices(this->fockMatrices, totalFockMatrix);

            arma::cx_mat totalOverlapMatrix = arma::kron(this->overlapMatrices.slice(i), arma::eye(2, 2));
            newOverlapMatrices = arma::join_slices(newOverlapMatrices, totalOverlapMatrix);
        }
        systemInfo.overlap    = newOverlapMatrices;
    }
    systemInfo.hamiltonian    = fockMatrices;

    if (debug){
        // Print contents
        std::cout << "Dim: " << std::endl;
        std::cout << systemInfo.ndim << "\n" << std::endl;

        std::cout << "Bravais lattice: " << std::endl;
        std::cout << bravaisLattice << "\n" << std::endl;

        std::cout << "Motif: " << std::endl;
        std::cout << motif << "\n" << std::endl;

        std::cout << "Orbitals: " << std::endl;
        std::cout << norbitals << "\n" << std::endl;

        std::cout << "Filling: " << systemInfo.filling << "\n" << std::endl;

        std::cout << "Hamiltonian: " << std::endl;
        std::cout << systemInfo.hamiltonian << "\n" << std::endl;

        std::cout << "Unit cells: " << std::endl;
        std::cout << systemInfo.bravaisVectors << "\n" << std::endl;

        std::cout << "Overlap: " << std::endl;
        std::cout << systemInfo.overlap << "\n" << std::endl;
    }
}

/**
 * Auxiliary routine to split a numeric string and get the corresponding length.
 * @param txt Reference to string to be splitted.
 * @param strs Reference to vector which will store the splitted string.
 * @return Length of splitted string.
 */
size_t split(const std::string &txt, std::vector<double> &strs)
{
    std::istringstream iss(txt);
    double token;
    size_t size = 0;
    while(iss >> token){
        strs.push_back(token);
        size++;
    }

    return size;
}
    
}


