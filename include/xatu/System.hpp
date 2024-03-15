#pragma once
#include "xatu/Lattice.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

/**
 * The System class is an abstract class that contains the minimum information relative to the system 
 * where we want to compute the exciton spectrum. It is defined as a sub-class of Lattice.
*/
class System : public Lattice {
    
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
        /* Constructors and destructor */
        System(const System&);
        System(const SystemConfiguration&);
        virtual ~System(){};

        void initializeSystemAttributes(const SystemConfiguration&);

        /* Setters */
        void setFilling(int);
        void setSystemName(std::string);

        /* Bloch Hamiltonian */
        virtual arma::cx_mat hamiltonian(arma::rowvec k) const = 0;
        virtual arma::cx_mat overlap(arma::rowvec k) const = 0;
        virtual void solveBands(arma::rowvec&, arma::vec&, arma::cx_mat&) const = 0;
        void solveBands(std::string) const;

        /* Modifiers */
        void addZeeman(double);    
};

}
