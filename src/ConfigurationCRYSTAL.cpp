#include "xatu/ConfigurationCRYSTAL.hpp"

#define SOC_STRING "to_be_defined_for_crystal23"
#define MAGNETIC_STRING "UNRESTRICTED OPEN SHELL"

namespace xatu {

/**
 * File constructor for ConfigurationCRYSTAL. It extracts the relevant information from CRYSTAL's .outp file.
 * @details This class is intended to be used with .outp files from the CRYSTAL code.
 * Since orbitals in CRYSTAL extend over several unit cells, the Fock matrices that define the
 * Hamiltonian also cover several unit cells. Therefore, one can specify how many unit cells to take
 * for the actual exciton calculation. 
 * @param outp_file Name of the .outp from the CRYSTAL post-scf calculation.
 * @param ncells Number of unit cells for which the Hamiltonian and overlap matrices are read from file.
 * @param isGTFmode True (false) indicates GTF (TB, resp.) mode. It determines whether hamiltonianMatrices or alpha/betaMatrices are 
 *        used in the case of magnetism without SOC.
 */
ConfigurationCRYSTAL::ConfigurationCRYSTAL(const std::string& outp_file, const int ncells, const bool isGTFmode) : ConfigurationBase{outp_file} {

    this->ncells_ = ncells;
    parseContent();
    Rlistsfun();

    if(!isGTFmode && MAGNETIC_FLAG && !SOC_FLAG){ // In the TB mode, build H(R) and S(R) matrices in spin space (in GTF mode, spin matrices are treated separately)
        arma::mat spinUpBlock = {{1, 0}, {0, 0}};
        arma::mat spinDownBlock = {{0, 0}, {0, 1}};
        arma::cx_cube newOverlapMatrices;
        for(int i = 0; i < ncells; i++){
            arma::cx_mat totalFockMatrix = arma::kron(this->alphaMatrices.slice(i), spinUpBlock) + arma::kron(this->betaMatrices.slice(i), spinDownBlock);
            this->hamiltonianMatrices_ = arma::join_slices(this->hamiltonianMatrices_, totalFockMatrix);
            arma::cx_mat totalOverlapMatrix = arma::kron(this->overlapMatrices.slice(i), arma::eye(2, 2));
            newOverlapMatrices = arma::join_slices(newOverlapMatrices, totalOverlapMatrix);
        }
        this->overlapMatrices_ = newOverlapMatrices;
    }

}

/**
 * Method to extract all the content from CRYSTAL's .outp file.
 * @details The lattice dimensionality is determined according to CRYSTAL's criterion for the 
 * direct lattice vectors: in 2D, R3 == [0 0 500] (Angstrom); and in 1D, in addition R2 = [0 0 500]. The 
 * hypothetical case of R3 and/or R2 == [0 0 0] is also assumed to lower the dimensionality.
 * @param ncells Number of unit cells to be parsed.
 * @return void.
 */
void ConfigurationCRYSTAL::parseContent(){

    std::vector<int> shellsPerSpecies;
    int totalElectrons;
    bool alpha_electrons = true;
    int countr = 0;
    std::string line;
    while(std::getline(m_file, line)){

        // Basis of Bravais vectors
        if (line.find("DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)") != std::string::npos){
            parseRbasis();
            countr++;
        }

        // Atomic motif
        else if(line.find("N. OF ATOMS PER CELL") != std::string::npos){
            int pos = line.find("N. OF ATOMS PER CELL");
            int strsize = strlen("N. OF ATOMS PER CELL");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> natoms_;
        }

        // N. shells (in total among all atoms in the unit cell)
        // else if(line.find("NUMBER OF SHELLS") != std::string::npos){
        //     int pos = line.find("NUMBER OF SHELLS");
        //     int strsize = strlen("NUMBER OF SHELLS");
        //     line = line.substr(pos + strsize, line.length());
        //     std::istringstream iss(line);
        //     iss >> nsh;
        // }

        // N. orbitals (total)
        else if(line.find("NUMBER OF AO") != std::string::npos){
            int pos = line.find("NUMBER OF AO");
            int strsize = strlen("NUMBER OF AO");
            line = line.substr(pos + strsize, line.length());
            std::istringstream iss(line);
            iss >> norbitals_;            
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
        // else if(line.find("CORE ELECTRONS PER CELL") != std::string::npos){
        //     int pos = line.find("CORE ELECTRONS PER CELL");
        //     int strsize = strlen("CORE ELECTRONS PER CELL");
        //     line = line.substr(pos + strsize, line.length());
        //     std::istringstream iss(line);
        //     iss >> coreElectrons;
        // }

        // 
        else if ((line.find("ATOM") != std::string::npos) && line.find("SHELL") != std::string::npos){
            if(natoms == 0){
                throw std::logic_error("Must parse number of atoms first");
            }
            shellsPerSpecies = parseAtoms();
        }

        // Parse atomic basis info
        else if(line.find("LOCAL ATOMIC FUNCTIONS BASIS SET") != std::string::npos){
            parseAtomicBasis(shellsPerSpecies);
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
                arma::colvec cell = arma::zeros<arma::colvec>(3);
                for (int i = 0; i < ndim; i++){
                    cell += Rbasis.col(i)*coefCombinations[i];
                }
                this->Rlist_ = arma::join_horiz(Rlist, cell);

                arma::cx_mat overlapMatrix = parseMatrix();
                this->overlapMatrices_ = arma::join_slices(this->overlapMatrices_, overlapMatrix);
            }
        }

        // Toggle SOC flag if calculation is done with spin-orbit coupling
        else if(line.find(SOC_STRING) != std::string::npos){
            this->SOC_FLAG_ = true;
        }

        else if(line.find(MAGNETIC_STRING) != std::string::npos){
            this->MAGNETIC_FLAG_ = true;
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
                if(MAGNETIC_FLAG && !SOC_FLAG){
                    if(alpha_electrons){
                        this->alphaMatrices_ = arma::join_slices(this->alphaMatrices_, fockMatrix);
                    }
                    else{
                        this->betaMatrices_ = arma::join_slices(this->betaMatrices_, fockMatrix);
                    }
                }
                else if(SOC_FLAG){
                    // Do something here (to be implemented when CRYSTAL prints H(R) matrices with SOC) 
                    continue;
                }
                else{ // No magnetism and no SOC
                    this->hamiltonianMatrices_ = arma::join_slices(this->hamiltonianMatrices_, fockMatrix);
                }
            }
        }
    }

    if (SOC_FLAG){
        this->filling_ = totalElectrons; 
        orbitalsPerSpecies_ *= 2;
        // Do something here (to be implemented when CRYSTAL prints H(R) matrices with SOC)
    }
    else if(MAGNETIC_FLAG){
        this->filling_ = totalElectrons;
        orbitalsPerSpecies_ *= 2;
    }
    else { // No magnetism and no SOC
        if(totalElectrons % 2 == 1){
            throw std::logic_error("CRYSTALConfiguration error: the system is a metal (odd number of electrons per cell and spin-independent Hamiltonian)");
        }
        else {
            this->filling_ = totalElectrons/2;
        }
    }

    if (countr==0){
        ndim_ = 0;
        throw std::logic_error("Finite systems are not yet implemented");
    } 
}

/**
 * Method to parse the Bravais basis vectors {R_1,..,R_ndim} from the file and store them by columns. 
 * @return void.
 */
void ConfigurationCRYSTAL::parseRbasis(){
    std::string line;
    std::vector<std::string> vectors;
    for(int i = 0; i < 3; i++){
        std::getline(m_file, line);
        vectors.push_back(line);     
    }
    this->Rbasis_ = parseVectors(vectors);
    extractDimension();
}

/**
 * Method to obtain the dimensionality (1, 2 or 3) of the lattice.
 * @return void. 
 */
void ConfigurationCRYSTAL::extractDimension(){
    arma::colvec R0 = arma::zeros<arma::colvec>(3);
    arma::colvec R2 = {0.0, 500.0, 0.0};   // Assuming that Rbasis is in Angstrom
    arma::colvec R3 = {0.0, 0.0, 500.0};   // Assuming that Rbasis is in Angstrom
                 
    if(approx_equal(Rbasis.col(2),R3,"absdiff",0.1) || approx_equal(Rbasis.col(2),R0,"absdiff",0.1) ){
        Rbasis_.shed_col(2);
        if(approx_equal(Rbasis.col(1),R2,"absdiff",0.1) || approx_equal(Rbasis.col(1),R0,"absdiff",0.1) ){
            Rbasis_.shed_col(1); 
        }
    }
    this->ndim_ = Rbasis.n_cols;
}

/**
 * Method to extract the motif, the chemical species and the number of shells per species.
 * @return std::vector<int> shellsPerSpecies, to be used later internally within parseContent. 
 */
std::vector<int> ConfigurationCRYSTAL::parseAtoms(){
    std::string line;
    int index, atomic_number, nsh, nspecies = 0;
    double x, y, z;
    std::string chemical_species;
    std::vector<int> shellsPerSpecies, atomic_number_ordering;
    std::vector<double> atom;
    std::map<std::string, int> chemical_species_to_index;
    std::vector<std::string> species;
    arma::mat motif = arma::zeros(4, natoms);

    std::getline(m_file, line); // Parse asterisks
    for(int i = 0; i < natoms; i++){
        std::getline(m_file, line);
        std::istringstream iss(line);
        iss >> index >> atomic_number >> chemical_species >> nsh >> x >> y >> z;

        if(std::find(species.begin(), species.end(), chemical_species) == species.end()){
            species.push_back(chemical_species);
            shellsPerSpecies.push_back(nsh);
            chemical_species_to_index[chemical_species] = chemical_species_to_index.size();
            atomic_number_ordering.push_back(atomic_number);
            nspecies++;
        }

        int index = chemical_species_to_index[chemical_species];
        atom = {x, y, z, (double)index};
        motif.col(i) = arma::colvec(atom);
    }
    this->motif_ = motif;
    this->atomic_number_ordering_ = atomic_number_ordering;
    this->nspecies_ = nspecies;
    return shellsPerSpecies;
}

/**
 * Method to parse the Fock and overlap matrices from the .outp file.
 * Note: apparently CRYSTAL prints **** instead of the corresponding row and column numbers when those are >= 10000, which doesn't disrupt this method.
 * @return arma::cx_mat.
 */
arma::cx_mat ConfigurationCRYSTAL::parseMatrix(){
    std::string colval, line;
    arma::cx_mat matrix = arma::zeros<arma::cx_mat>(norbitals, norbitals);
    double  coef;
    int     subCol = 0;  // Index for the column in the current block
    int32_t subRow = 0;  // Index for the row in the current block
    int32_t blockIndex = -10;  // The number of the current block (starting at 0) times 10
    while(std::getline(m_file, line)){
        if (line.empty()){ // Start of a new sub-block (including the first one)
            subRow = 0;
            blockIndex += 10;
            std::getline(m_file, line); // Get next line, which contains column indices or ****
            std::getline(m_file, line); // The next line is blank except for the last sub-block. If so, skip to the bottom of the while loop
            if (line.empty()){continue;}   
        }

        subCol = 0;
        std::istringstream iss(line);
        iss >> colval;
        while(iss >> coef){ // Read matrix elements in the current line (or row)
            matrix(blockIndex + subRow, blockIndex + subCol) = coef;
            subCol++;
        }
        subRow++;

        if(static_cast<uint32_t>(blockIndex + subCol) == norbitals){
            return matrix;
        }
    }
    throw std::logic_error("ERROR ConfigurationCRYSTAL::parseMatrix. Unable to read Overlap or Hamiltonian matrices, check the .outp file.");

}

/**
 * Method to extract the details of the basis used in the CRYSTAL calculation. Only relevant in TB mode.
 * @details This method extracts all the orbitals per chemical species and
 * the corresponding coefficients of the gaussian expansion.
 * @return void.
 */
void ConfigurationCRYSTAL::parseAtomicBasis(const std::vector<int>& shellsPerSpecies){
    std::string line, chemical_species;
    int natom, nspecies = 0;
    int norb; 
    int cumOrbitals = 0;
    std::string shellType;
    // double exponent, sCoef, pCoef, dCoef;
    // std::vector<double> coefs;
    // cube_vector gaussianCoefficients;
    std::vector<std::string> species;
    std::vector<int> orbitalsPerSpecies;

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
            cumOrbitals += orbitalsPerSpecies[motif.col(atomIndex)(3)]; // Track cumulative num. orbitals
            continue; // Skip already parsed
        }

        // gaussianCoefficients.clear();

        for(int shellIndex = 0; shellIndex < shellsPerSpecies[nspecies]; shellIndex++){

            std::getline(m_file, line);
            std::istringstream iss(line);
            
            iss >> norb >> shellType;
            if(shellType == "-"){
                iss >> norb >> shellType;
            }

            // std::vector<std::vector<double>> coefList;
            long int previousLine = m_file.tellg(); // Store beginning of next line
            while (std::getline(m_file, line)){
                std::vector<double> tokenized_line;
                split(line, tokenized_line);
                if (tokenized_line.size() != 4){
                    m_file.seekg(previousLine); // Restore line
                    break;
                }
                
                // std::istringstream iss(line);
                // iss >> exponent >> sCoef >> pCoef >> dCoef;
                // coefs = {exponent, sCoef, pCoef, dCoef};

                // coefList.push_back(coefs);

                previousLine = m_file.tellg(); // Store beginning of next line
            }
            // gaussianCoefficients.push_back(coefList);

            // this->gaussianCoefficients[nspecies] = gaussianCoefficients;
        }
        orbitalsPerSpecies.push_back(norb - cumOrbitals);
        cumOrbitals = norb;

        nspecies++;
    }
    this->orbitalsPerSpecies_ = arma::conv_to<arma::urowvec>::from(orbitalsPerSpecies);

}

/**
 * Auxiliary routine to split a numeric string and get the corresponding length.
 * @param txt Reference to string to be splitted.
 * @param strs Reference to vector which will store the splitted string.
 * @return Length of splitted string.
 */
size_t ConfigurationCRYSTAL::split(const std::string &txt, std::vector<double> &strs)
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
    
/**
 * Method to remove the unpaired Bravais vectors from Rlist and the Hamiltonian and overlap matrices. This ensures
 * that the k-dependent matrices are hermitian, irrespective of the supplied list of direct lattice vectors.
 * @return void
 */
void ConfigurationCRYSTAL::Rlistsfun(){

    int ncells_copy = ncells;
    for(int RindOpp = 0; RindOpp < ncells_copy; RindOpp++){ //Vector for which the opposite is searched
        int countr = 0;
        for(int Rind = 0; Rind < ncells_copy; Rind++){ //Looking for the opposite of RindOpp among all vectors
            arma::colvec Rsum = Rlist.col(Rind) + Rlist.col(RindOpp);
            if( Rsum.is_zero(0.001) ){
                countr++;
            }
        }
        if(countr == 0){
            std::cerr << "WARNING! Bravais vector number " << (RindOpp+1) <<
                " has no opposite in the list of vectors, and has been discarded to avoid non-hermiticity." << std::endl;

            Rlist_.shed_col(RindOpp);
            overlapMatrices_.shed_slice(RindOpp);
            if(!MAGNETIC_FLAG || SOC_FLAG){
                hamiltonianMatrices_.shed_slice(RindOpp);
            } 
            else{ //magnetism without SOC => 2 separate spin block matrices
                alphaMatrices_.shed_slice(RindOpp);
                betaMatrices_.shed_slice(RindOpp);
            }
            ncells_ -= 1;
        }
    }
}

}


