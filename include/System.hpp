#pragma once
#include "Crystal.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class System : public Crystal{
    
    //// Attributes
    private:
        int norbitals_, basisdim_;
        arma::urowvec orbitals;
        arma::cx_cube hamiltonianMatrices;
        arma::cx_cube overlapMatrices;

    // Const references to expose relevant attributes in a read-only way
    public:
        const int& norbitals = norbitals_;
        const int& basisdim = basisdim_;

    //// Methods
    public:
        /* Constructor and destructor */
        System(std::string);    
        ~System();


        /* Bloch Hamiltonian */
        arma::cx_mat hamiltonian(arma::rowvec k, bool isTriangular = false);
        arma::cx_mat overlap(arma::rowvec k, bool isTriangular = false);

        /* Expected value of spin components */
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);

        /* Routines for DoS calculation */
        double densityOfStates(double, double, const arma::mat&);
        void writeDensityOfStates(const arma::mat&, double, FILE*);

    protected:
        /* Routines for DoS calculation */
        std::complex<double> rGreenF(double, double, double);
    
    private:
        void initializeSystemAttributes(const SystemConfiguration&);
};
