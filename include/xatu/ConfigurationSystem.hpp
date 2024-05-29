#pragma once
#include <complex>
#include <armadillo>
#include "xatu/ConfigurationBase.hpp"

typedef std::vector<std::vector<std::vector<double>>> cube_vector;

namespace xatu {

/**
 * The ConfigurationSystem class is an extension of ConfigurationBase declaring common attributes
 * for both the TB and Gaussian modes
 */
class ConfigurationSystem : public virtual ConfigurationBase {

    protected:  
        // Lattice dimensionality (0,1,2,3)
        int ndim_;
        // Basis of Bravais vectors (R1,R2,R3), stored by columns: (3,ndim)
        arma::mat Rbasis_;
        // Matrix containing the positions of the atoms of the motif by columns: (4,natoms). Each row has the format {x,y,z; species},
        // where x,y,z are in Angstrom and the last element is an index (0,1,...) labelling the chemical species
        arma::mat motif_;
        // Number of filled bands, considering that each band holds 2 electrons at each k-point for spin-independent Hamiltonians
        int filling_;
        // List of direct lattice vectors {Ri} associated to the stored Hamiltonian matrices, stored by columns: (3,ncells).
        // The Hamiltonian (and Overlap, if reading from CRYSTAL) matrices are stored in the same order as these vectors
        arma::mat Rlist_;
        // Number of unit cells for which the Hamiltonian (and possibly overlap) matrices are read from file
        int ncells_;
        // Hamiltonian or Fock matrices stored as a cube, where each slice represents an element R in Rlist (in the same ordering) 
        // and contains the matrix H(R) 
        arma::cx_cube hamiltonianMatrices_;
        // Overlap stored as a cube, where each slice represents an element R in Rlist (in the same ordering) 
        // and contains the matrix S(R) 
        arma::cx_cube overlapMatrices_;
        // Vector storing the number of orbitals for each chemical species. This is only used in TB mode (including TB with CRYSTAL,
        // in which case each (l,m) pair constitutes a different orbital)
        arma::urowvec orbitalsPerSpecies_;
    
    public:  // Const references to attributes (read-only)
        const int& ndim = ndim_;
        const arma::mat& Rbasis = Rbasis_;
        const arma::mat& motif = motif_;
        const int& filling = filling_;
        const arma::mat& Rlist = Rlist_;
        const int& ncells = ncells_;
        const arma::cx_cube& hamiltonianMatrices = hamiltonianMatrices_;
        const arma::cx_cube& overlapMatrices = overlapMatrices_;
        const arma::urowvec& orbitalsPerSpecies = orbitalsPerSpecies_;
        
    protected:
        ConfigurationSystem() = default;

    public:
        // Method to parse several three-dimensional vectors into a matrix, where they are stored by columns
        arma::mat parseVectors(const std::vector<std::string>&);
        // Optional method to print all the attributes above. printH (printS) print all the H(R) (S(R), respectively) matrices
        void printConfiguration(const bool printH = false, const bool printS = false);

};

}