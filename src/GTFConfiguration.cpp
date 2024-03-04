#include <algorithm>
#include "xatu/GTFConfiguration.hpp"


namespace xatu {

/**
 * File constructor for GTFConfiguration. It extracts the relevant information (exponents and contraction 
 * coefficients per shell per atomic species) from the basis sets file and stores it in an adequate format.
 * @details This class is intended to be used with the basis sets file in addition to the .outp file from the 
 * CRYSTAL code (see CRYSTALConfiguration.cpp for details on the latter). The former must be given in CRYSTAL 
 * format and contain both the basis employed in the self-consistent DFT calculation (under a line with the 
 * string: SCF BASIS) in addition to a larger auxiliary basis employed for fitting the density (under a line 
 * with the string: AUXILIARY BASIS). The initial occupancies have no impact. Pseudo-potentials are ignored, 
 * but it is advisable to be consistent with the PP choices in the self-consistency.
 * @param bases_file Name of the file containing both Gaussian basis sets.
 * @param file Name of the .outp from the CRYSTAL SCF calculation.
 * @param ncells Number of unit cells for which the Hamiltonian and overlap matrices are read from file.
 */
GTFConfiguration::GTFConfiguration(std::string bases_file, std::string file, int ncells) : CRYSTALConfiguration{file, ncells}, ConfigurationBase{file} {
    
    this->b_filename = bases_file;
   
    if(bases_file.empty()){
        throw std::invalid_argument("GTFConfiguration: bases_file must not be empty");
    }
    b_file.open(bases_file.c_str());
    if(!b_file.is_open()){
        throw std::invalid_argument("GTFConfiguration: bases file does not exist");
    }

   arma::mat Rlist = (bravaisVectors.cols(0,ndim-1)).t();
   this-> Rlist = Rlist;
   parseBases();
}

/**
 * Method to parse both the DFT basis used in the self-consistency, and the auxiliary basis
 * @return void
 */
void GTFConfiguration::parseBases(){ 
    int countr = 0;
    std::string line;
    while(std::getline(b_file, line)){
        if (line.find("SCF BASIS") != std::string::npos){
            parseBasis(true);
            countr++;
        }

        if (line.find("AUXILIARY BASIS") != std::string::npos){
            parseBasis(false);
            countr++;
        }  
    } 
    if (countr!=2){
        throw std::logic_error("The basis sets file is not valid. Make sure that the bases are preceded by the SCF BASIS and AUXILIARY BASIS lines, and that the ordering of this and the .outp file is correct");
    }
}

/**
 * Method to parse a given basis already identified in the file
 * @return void
 */
void GTFConfiguration::parseBasis(bool basis_id){
    cube_vector shells_all_species;
    std::vector<std::vector<int>> L_all_species, nG_all_species; 
    std::vector<int> nshells_all_species(nspecies,0);
    for (int s = 0; s < nspecies; s++){ //initialize the first dimension with empty entries, so that they can be reordered later
        shells_all_species.push_back({});
        L_all_species.push_back({});
        nG_all_species.push_back({});
    }
    int atomic_number, nshells; 
    for (int s = 0; s < nspecies; s++){
        std::vector<std::vector<double>> shells_in_species; //vector assigning the vector of contracted GTF to each shell
        std::vector<int> L_in_species; //list of L (ang. mom. quant. num.) corresponding to each shell, for a given atomic species
        std::vector<int> nG_in_species; // list of nG (number of contracted Gaussians) in each shell, for a given atomic species
        int nshells_SP = 0;  //counts the number of SP shells in a given species, to increase the nshells_in_all_species attribute later
        int bs0, bs1, bs2;
        float bs3, bs4;
        double exponent, contraction, contractionSP;
        std::string exponent_s, contraction_s, contractionSP_s;
        std::string line1;
        std::getline(b_file, line1);
        std::istringstream iss(line1);
        iss >> atomic_number >> nshells;

        int itr_order;
        std::vector<int>::iterator itr = find(atomic_number_ordering.begin(),atomic_number_ordering.end(),atomic_number);
        if(itr == atomic_number_ordering.end()){ //reorder the species within the bases file to match the ordering in the .outp
            std::vector<int>::iterator itr1 = find(atomic_number_ordering.begin(),atomic_number_ordering.end(),atomic_number + 200);
            std::vector<int>::iterator itr2 = find(atomic_number_ordering.begin(),atomic_number_ordering.end(),atomic_number % 200);
            if(itr1 == atomic_number_ordering.end() && itr2 == atomic_number_ordering.end()){ //check if it was a PP mismatch
                 throw std::logic_error("The atomic numbers in the basis sets file don't match those of the DFT calculation");
            } else if(itr1 != atomic_number_ordering.end()){
                 itr_order = itr1 - atomic_number_ordering.begin();
            } else {
                 itr_order = itr2 - atomic_number_ordering.begin();
            }
        } else {
            itr_order = itr - atomic_number_ordering.begin(); 
        }

        if (atomic_number > 200){          //skip pseudo-potential
            std::getline(b_file, line1);
            if (line1.find("INPUT") != std::string::npos){ 
                skip_PP(); 
            }
            else if (line1.find("INPSOC") != std::string::npos){
                std::getline(b_file, line1);
                skip_PP();
            }
        }

        for(int i = 0; i < nshells; i++){   //iterate over (unfolded) shells
            std::vector<double> ContGs_in_Shell, ContGs_in_Shell_SP;  //list(s) of the contracted GTF for a given shell: {alpha_1,d_1,..,alpha_n,d_n}
            std::getline(b_file, line1);
            std::istringstream iss(line1);
            iss >> bs0 >> bs1 >> bs2 >> bs3 >> bs4;

            if (bs0 != 0 || bs4!= 1.0 ) throw std::logic_error("ITYB and SCAL parameters must be 0 and 1, see CRYSTAL manual");

            if (bs1 == 1){      //unfold SP into separate S and P orbitals
                nshells_SP++;
                nG_in_species.push_back(bs2);
                nG_in_species.push_back(bs2);
                L_in_species.push_back(0);
                L_in_species.push_back(1);
                for(int j = 0; j < bs2; j++){   //iterate over contracted GTF within a shell
                std::getline(b_file, line1);
                std::istringstream iss(line1);
                iss >> exponent_s >> contraction_s >> contractionSP_s;

                exponent = replace_D_2_double(exponent_s);
                contraction = replace_D_2_double(contraction_s);
                contractionSP = replace_D_2_double(contractionSP_s);

                ContGs_in_Shell.insert(ContGs_in_Shell.end(), {exponent, contraction});
                ContGs_in_Shell_SP.insert(ContGs_in_Shell_SP.end(), {exponent, contractionSP});
                }
                shells_in_species.push_back(ContGs_in_Shell);
                shells_in_species.push_back(ContGs_in_Shell_SP);
            } 
            else {
                nG_in_species.push_back(bs2);
                L_in_species.push_back( std::max(bs1-1,0) );
                for(int j = 0; j < bs2; j++){   //iterate over contracted GTF within a shell
                    std::getline(b_file, line1);
                    std::istringstream iss(line1);
                    iss >> exponent_s >> contraction_s;
                    
                    exponent = replace_D_2_double(exponent_s);
                    contraction = replace_D_2_double(contraction_s);

                    ContGs_in_Shell.insert(ContGs_in_Shell.end(), {exponent, contraction});
                }
                shells_in_species.push_back(ContGs_in_Shell);
            }
        }
        
        nshells_all_species[itr_order]  = nshells + nshells_SP;
        nG_all_species[itr_order]       = nG_in_species;
        L_all_species[itr_order]        = L_in_species;
        shells_all_species[itr_order]   = shells_in_species;
    }
    if(basis_id){ //basis_id == true => SCF basis; basis_id == false => auxiliary basis
        this->nshells_all_species_SCF = nshells_all_species;
        this->nG_all_species_SCF      = nG_all_species;
        this->L_all_species_SCF       = L_all_species;
        this->shells_all_species_SCF  = shells_all_species;
    } else {
        this->nshells_all_species_AUX = nshells_all_species;
        this->nG_all_species_AUX      = nG_all_species;
        this->L_all_species_AUX       = L_all_species;
        this->shells_all_species_AUX  = shells_all_species;
    }
}

/**
 * Method to skip the pseudo-potentials when parsing the bases 
 * @return void
 */
void GTFConfiguration::skip_PP(){ 
    float pp0;
    int pp1, pp2, pp3, pp4, pp5, pp6;
    std::string line;
    std::getline(b_file, line);
    std::istringstream iss(line);
    iss >> pp0 >> pp1 >> pp2 >> pp3 >> pp4 >> pp5 >> pp6;
    for(int i = 0; i < pp1+pp2+pp3+pp4+pp5+pp6; i++) std::getline(b_file, line);
}

/**
 * Method to replace Fortran's scientific notation (D) with the standard (E). 
 * The string is then converted to double.
 * @return double
 */
double GTFConfiguration::replace_D_2_double(std::string S){ 
    std::replace(S.begin(), S.end(), 'D', 'E');
    return std::stod(S);
}

}