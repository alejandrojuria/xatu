#pragma once
#include <map>
#include "xatu/CRYSTALConfiguration.hpp"

namespace xatu {

class GTFConfiguration : public CRYSTALConfiguration {

    public:
        // (3 x ncells) matrix of Bravais vectors, stored by columns.
        arma::mat Rlist;
        // Map of ncells entries, the n-th (n=0,..,ncells-1) of which is the index of -R_{n} (minus the n-th Bravais vector in Rlist).
        std::map<int,int> RlistOpposites;
        // Vector with the indices of the unpaired Bravais vectors, i.e. those with their opposites not present in Rlist.
        std::vector<int> R_unpaired;
        // Vector assigning the vector of shells (itself a vector of contracted GTF) to each atomic species.
        // One for the DFT basis and one for the auxiliary basis.
        cube_vector shells_all_species_SCF, shells_all_species_AUX;
        // Vector assigning the vector of L (ang. mom. quant. num.) to each atomic species.
        std::vector<std::vector<int>> L_all_species_SCF, L_all_species_AUX;
        // Vector assigning the vector of nG (number of contracted Gaussians) to each atomic species.
        std::vector<std::vector<int>> nG_all_species_SCF, nG_all_species_AUX;
        // Vector storing the number shells for each atomic species.
        std::vector<int> nshells_all_species_SCF, nshells_all_species_AUX;
        // Vector storing the number orbitals for each atomic species.
        std::vector<int> norbs_all_species_SCF, norbs_all_species_AUX;
        // Name of the basis sets file.
        std::string bases_filename;
        // Redefinition of the Hamiltonian matrices as real if SOC is absent.
        // arma::cube fockMatrices, alphaMatrices, betaMatrices;

    protected:
        // Open basis sets file.
        std::ifstream b_file;

    public:
        GTFConfiguration(std::string bases_file, std::string outp_file, int ncells);
        ~GTFConfiguration(){};

    protected:  
        // Method to convert possible Fortran convention for scientific notation in the bases file to that of normal people :)
        double replace_D_2_double(std::string&);

    private:
        // Preliminary parsing method common to both basis sets.
        void parseBases();
        // Method to parse a given basis (bool = true => SCF, bool = false => AUX).
        void parseBasis(bool);
        // Method to skip possible pseudo-potentials in the basis sets file (they are irrelevant to the matrix elements).
        void skip_PP();
        // Method to build the Rlist and RlistOpposites attributes, relative to the selected Bravais vectors.
        void Rlistsfun(const int ncells);
        // Method to redefine the Hamiltonian matrices as real if SOC is absent.
        // void makeHreal();
        
};

}