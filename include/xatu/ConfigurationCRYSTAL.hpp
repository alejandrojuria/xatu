#pragma once
#include "xatu/ConfigurationSystem.hpp"

namespace xatu {

/**
 * The ConfigurationCRYSTAL class is used to parse CRYSTAL's .outp files, which can be used both 
 * in the TB and Gaussian modes. This file is a post-scf processing which must contain H(R) and S(R).
 * This class will remove the unpaired Bravais vectors from Rlist and the Hamiltonian and overlap matrices.
 */
class ConfigurationCRYSTAL final : public ConfigurationSystem {

    protected:
        // Number of atoms in the unit cell
        int natoms_;
        // Number of chemical species in the unit cell
        int nspecies_;
        // Single-particle matrix dimension, which equals the number of orbitals (with fixed l,m in the SCF basis for Gaussian mode)
        uint32_t norbitals_;
        // Vector with the atomic number (+200 if pseudo-potential) of each chemical species, in the ordering displayed in the .outp file
        std::vector<int> atomic_number_ordering_;
        // Boolean which indicates whether the single-particle Hamiltonian includes SOC
        bool SOC_FLAG_ = false;
        // Boolean which indicates whether the single-particle Hamiltonian includes magnetic terms. Irrelevant under SOC
        bool MAGNETIC_FLAG_ = false;
        // Substitutes of hamiltonianMatrices for GTF mode under magnetism but no SOC, i.e they are initialized 
        // if(MAGNETIC_FLAG && !SOC_FLAG). Alpha (beta) refer to up-up (down-down, respectively) spin blocks
        arma::cx_cube alphaMatrices_, betaMatrices_;

    public:  // Const references to attributes (read-only)
        const int& natoms = natoms_;
        const int& nspecies = nspecies_;
        const uint32_t& norbitals = norbitals_;
        const std::vector<int>& atomic_number_ordering = atomic_number_ordering_;
        const bool& SOC_FLAG = SOC_FLAG_;
        const bool& MAGNETIC_FLAG = MAGNETIC_FLAG_;
        const arma::cx_cube& alphaMatrices = alphaMatrices_;
        const arma::cx_cube& betaMatrices = betaMatrices_;

    protected:
        ConfigurationCRYSTAL() = default;
    public:
        ConfigurationCRYSTAL(const std::string& outp_file, const int ncells, const bool isGTFmode);

    private:
        // Central method to parse the content from the configuration file
        void parseContent() override;
        void parseRbasis();
        void extractDimension();
        std::vector<int> parseAtoms();
        arma::cx_mat parseMatrix();
        // Method to remove the unpaired Bravais vectors from Rlist and the Hamiltonian and overlap matrices. This ensures
        // that the k-dependent matrices are hermitian, irrespective of the supplied list of direct lattice vectors.
        void Rlistsfun();
        // Relevant only in TB mode
        void parseAtomicBasis(const std::vector<int>&);
        size_t split(const std::string&, std::vector<double>&);
        
};

}
