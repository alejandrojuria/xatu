#pragma once
#include <armadillo>
#include "xatu/ConfigurationBase.hpp"

namespace xatu {

/**
 * The SystemConfiguration class is an extension of ConfigurationBase to parse 
 * specifically system configuration files.
 */
class SystemConfiguration : public virtual ConfigurationBase {

    struct configuration {
        /// Dimension of the system. Given by the number of bravais basis vectors.
        int ndim;
        /// Number of electrons of the system.
        double filling;
        /// Bravais basis vectors, ordered by rows.
        arma::mat bravaisLattice;
        /// Matrix containing the positions of the atoms of the motif by rows. 
        /// Each row has the format {x,y,z; species}, where the last elements is an index representing the chemical species.
        arma::mat motif;
        /// Fock matrices that form the Bloch Hamiltonian.
        arma::cx_cube hamiltonian;
        /// Overlap matrices in real-space (if any).
        arma::cx_cube overlap;
        /// Bravais vectors associated to the stored Fock matrices.
        /// Note that the Fock (overlap) matrices are stored in the same order as this vectors.
        arma::mat bravaisVectors;
        /// Vector storing the number of orbitals for each chemical species.
        arma::urowvec norbitals;
    };
        
    public:
        /// Parsed information of file already fully structured.
        configuration systemInfo;

    protected:
        SystemConfiguration();
    public:
        SystemConfiguration(std::string);
        ~SystemConfiguration(){};

        void printConfiguration(std::ostream& stream = std::cout) const;
        friend std::ostream& operator<< (std::ostream&, const SystemConfiguration&);
        
    protected:
        arma::mat parseVectors(std::vector<std::string>&);
        arma::cx_cube parseMatrices(std::vector<std::string>&);
        arma::mat parseMotif(std::vector<std::string>&);
        arma::urowvec parseOrbitals(std::vector<std::string>&);
        void parseContent() override;
        void checkContentCoherence();

};

}