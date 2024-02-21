#pragma once
#include "xatu/Lattice.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

/**
 * The System class contains all information regarding the system where we want to compute
 * the exciton spectrum. It is defined as a sub-class of Crystal.
*/
class System : public Lattice{
    
    //// Attributes
    protected:
        int basisdim_, filling_, fermiLevel_;
        std::string systemName;

        arma::urowvec orbitals_;
        arma::cx_cube hamiltonianMatrices_;
        arma::cx_cube overlapMatrices_;

    // Const references to expose relevant attributes in a read-only way
    public:
        // Returns vector with orbitals per chemical species
        const arma::urowvec& orbitals = orbitals_;
        // Returns size of basis used/hamiltonian dimension
        const int& basisdim = basisdim_;
        // Returns value of system filling (total number of electrons)
        const int& filling = filling_;
        // Returns index of band corresponding to Fermi level
        const int& fermiLevel = fermiLevel_;
        // Fock matrices
        const arma::cx_cube& hamiltonianMatrices = hamiltonianMatrices_;
        // Overlap matrices
        const arma::cx_cube& overlapMatrices = overlapMatrices_;


    //// Methods
    public:
        /* Constructor and destructor */

        System();
        System(const System&);
        System(const SystemConfiguration&);  

        void setFilling(int);

        /* Bloch Hamiltonian */

        arma::cx_mat hamiltonian(arma::rowvec k, bool isTriangular = false) const;
        arma::cx_mat overlap(arma::rowvec k, bool isTriangular = false) const;
        void solveBands(arma::rowvec&, arma::vec&, arma::cx_mat&, bool triangular = false) const;
        void solveBands(std::string, bool triangular = false) const;

        /* Modifiers */
        void addZeeman(double);

        /* Expected value of spin components */
        
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);      

        arma::cx_vec velocity(const arma::rowvec, int, int) const;  
    
        void initializeSystemAttributes(const SystemConfiguration&);

    private:
        void orthogonalize(const arma::rowvec&, arma::cx_mat&, bool) const;
        void orthogonalize_hamiltonian(const arma::rowvec&, arma::cx_mat&, bool) const;
};

}
