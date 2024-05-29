#pragma once
#include "xatu/ConfigurationCRYSTAL.hpp"

namespace xatu {

/**
 * The ConfigurationGTF class is used to parse the file containing both the SCF and the AUXILIARY Gaussian basis sets. 
 * This class is intended to be used with the basis sets file in addition to the .outp file from the 
 * CRYSTAL code (see ConfigurationCRYSTAL.cpp for details on the latter). The former must be given in CRYSTAL 
 * format and contain both the basis employed in the self-consistent DFT calculation (under a line with the 
 * string: SCF BASIS) in addition to a larger auxiliary basis employed for fitting the density (under a line 
 * with the string: AUXILIARY BASIS). The initial occupancies have no impact. Pseudo-potentials are ignored, 
 * but it is advisable to be consistent with the PP choices in the self-consistency.
 * Exclusive to the Gaussian mode
 */
class ConfigurationGTF final : public ConfigurationBase {

    protected:
        // Vector assigning the vector of shells (itself a vector of contracted GTF) to each atomic species
        // One for the SCF basis and one for the AUX basis.
        cube_vector shells_all_species_SCF_, shells_all_species_AUX_;
        // Vector assigning the vector of l (angular momentum quantum number) to each atomic species
        std::vector<std::vector<int>> L_all_species_SCF_, L_all_species_AUX_;
        // Vector assigning the vector of nG (number of contracted Gaussians) to each atomic species
        std::vector<std::vector<int>> nG_all_species_SCF_, nG_all_species_AUX_;
        // Vector storing the number shells (or entries in each basis set) for each atomic species
        std::vector<int> nshells_all_species_SCF_, nshells_all_species_AUX_;
        // Vector storing the number orbitals (defined by (l,m) and Gaussian coefficients) for each atomic species
        std::vector<int> norbs_all_species_SCF_, norbs_all_species_AUX_;
    // Attributes copied from ConfigurationCRYSTAL
        // Number of chemical species in the unit cell
        int nspecies_;
        // Vector with the atomic number (+200 if pseudo-potential) of each chemical species, in the ordering displayed in the .outp file
        std::vector<int> atomic_number_ordering_;

    public:  // Const references to attributes (read-only)
        const cube_vector& shells_all_species_SCF = shells_all_species_SCF_;
        const cube_vector& shells_all_species_AUX = shells_all_species_AUX_;
        const std::vector<std::vector<int>>& L_all_species_SCF = L_all_species_SCF_;
        const std::vector<std::vector<int>>& L_all_species_AUX = L_all_species_AUX_;
        const std::vector<std::vector<int>>& nG_all_species_SCF = nG_all_species_SCF_;
        const std::vector<std::vector<int>>& nG_all_species_AUX = nG_all_species_AUX_;
        const std::vector<int>& nshells_all_species_SCF = nshells_all_species_SCF_;
        const std::vector<int>& nshells_all_species_AUX = nshells_all_species_AUX_;
        const std::vector<int>& norbs_all_species_SCF = norbs_all_species_SCF_;
        const std::vector<int>& norbs_all_species_AUX = norbs_all_species_AUX_;
        const int& nspecies = nspecies_;
        const std::vector<int>& atomic_number_ordering = atomic_number_ordering_;

    public:
        ConfigurationGTF(const int nspecies, const std::vector<int>& atomic_number_ordering, const std::string& bases_file);
        ConfigurationGTF(const ConfigurationCRYSTAL&, const std::string& bases_file);

    private:
        // Preliminary parsing method common to both basis sets
        void parseContent() override;
        // Method to parse a given basis (bool = true => SCF, bool = false => AUX)
        void parseBasis(const bool);
        // Method to skip possible pseudo-potentials in the basis sets file (they are irrelevant to the matrix elements)
        void skip_PP();
        // Method to convert possible Fortran convention for scientific notation in the bases file to C++ convention
        double replace_D_2_double(std::string&);
        
};

}