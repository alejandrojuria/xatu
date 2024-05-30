#pragma once
#include <stdlib.h>
#include "xatu/ConfigurationCRYSTAL.hpp"
#include "xatu/Lattice.hpp"

namespace xatu {

/**
 * The System class is an abstract class that contains the information, common to both TB and Gaussian modes, relative to the 
 * system where we want to compute the exciton spectrum. It is defined as a sub-class of Lattice
*/
class System : public Lattice {

    public:
        // Imaginary unit i
        static constexpr std::complex<double> imag {0., 1.};
        // Hamiltonian or Fock matrices stored as a cube, where each slice represents an element R in Rlist (in the same ordering) 
        // and contains the matrix H(R) 
        const arma::cx_cube* ptr_hamiltonianMatrices;
        // Substitutes of hamiltonianMatrices under magnetism but no SOC. Alpha (beta) refer to up-up (down-down, respectively) spin blocks
        const arma::cx_cube* ptr_alphaMatrices; 
        const arma::cx_cube* ptr_betaMatrices;
        // Overlap stored as a cube, where each slice represents an element R in Rlist (in the same ordering) 
        // and contains the matrix S(R) 
        const arma::cx_cube* ptr_overlapMatrices;
    
    protected:
        // String labelling the system
        std::string systemName;
        // Index of the highest occupied band, starting at 0. The system is assumed to be non-metallic
        int highestValenceBand_;
        // List of direct lattice vectors {Ri} associated to the stored Hamiltonian matrices, stored by columns: (3,ncells).
        // The Hamiltonian (and Overlap, if reading from CRYSTAL) matrices are stored in the same order as these vectors
        arma::mat Rlist_;
        // Number of unit cells for which the Hamiltonian (and possibly overlap) matrices are stored 
        int ncells_;
        // Single-particle matrix dimension, which equals the number of orbitals (with fixed l,m in the SCF basis for Gaussian mode)
        uint32_t norbitals_; 
        // List of k-points that define the basis, stored by (3D) columns
        arma::mat kpoints_;
        // Number of k-points in kpoints
        uint32_t nk_;
    
    public:  // Const references to attributes (read-only)
        const int& highestValenceBand = highestValenceBand_;
        const arma::mat& Rlist        = Rlist_;
        const int& ncells             = ncells_;
        const uint32_t& norbitals     = norbitals_;
        const arma::mat& kpoints = kpoints_;
        const uint32_t& nk = nk_;

    protected:
        System(); 
    public:
        System(const ConfigurationSystem&);
        System(const ConfigurationCRYSTAL&);
        System(const System&) = default; 

    public:
        void setSystemName(const std::string& systemName);

        /* Bloch Hamiltonian */
        virtual arma::cx_mat hamiltonian(const arma::colvec& k) const = 0;
        virtual arma::cx_mat overlap(const arma::colvec& k) const = 0;
        virtual void solveBands(const arma::colvec&, arma::vec&, arma::cx_mat&) const = 0;
        void solveBands(const std::string& kpointsfile) const;
          
};

}
