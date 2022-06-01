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
        int basisdim_, fermiLevel_;
        double filling_;

        arma::urowvec orbitals_;
        arma::cx_cube hamiltonianMatrices;
        arma::cx_cube overlapMatrices;

    // Const references to expose relevant attributes in a read-only way
    public:
        // Returns vector with orbitals per chemical species
        const arma::urowvec& orbitals = orbitals_;
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
    
        void initializeSystemAttributes(const SystemConfiguration&);
};
