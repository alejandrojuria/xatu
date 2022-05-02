#pragma once
#include "Crystal.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class System : public Crystal{
    
    //// Attributes
    protected:
        int norbitals_, basisdim_, fermiLevel_;
        double filling_;

        arma::urowvec orbitals;
        arma::cx_cube hamiltonianMatrices;
        arma::cx_cube overlapMatrices;

    // Const references to expose relevant attributes in a read-only way
    public:
        // Returns number of orbitals per atom
        const int& norbitals = norbitals_;
        // Returns size of basis used/hamiltonian dimension
        const int& basisdim = basisdim_;
        // Returns value of system filling
        const double& filling = filling_;
        // Returns band number corresponding to the Fermi level
        const int& fermiLevel = fermiLevel_;


    //// Methods
    public:
        /* Constructor and destructor */
        System(){};
        System(std::string);  
        ~System();

        void setFilling(double);

        /* Bloch Hamiltonian */
        arma::cx_mat hamiltonian(arma::rowvec k, bool isTriangular = false);
        arma::cx_mat overlap(arma::rowvec k, bool isTriangular = false);
        

        /* Expected value of spin components */
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);        
    
    private:
        void initializeSystemAttributes(const SystemConfiguration&);
};
