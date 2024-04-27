#pragma once
#include "xatu/System.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

/**
 * The System class contains all information regarding the system where we want to compute
 * the exciton spectrum. It is defined as a sub-class of System.
*/
class SystemTB : public System {
    
    //// Attributes
    protected:
        bool isTriangular_ = false;
        bool isAU_ = false;

    // Const references to expose relevant attributes in a read-only way
    public:
        // Returns true if the system matrices (hamiltonian and overlap) are triangular
        const bool& isTriangular = isTriangular_;
        // Returns true if the system is in atomic units 
        const bool& isAU = isAU_;
    
    //// Methods
    /* Constructors and destructor */
    protected:
        SystemTB(){}; // Protected so that SystemTB cannot be initialized without parameters
    public:
        SystemTB(const SystemTB&);
        SystemTB(const SystemConfiguration&);  
        ~SystemTB(){};

        /* Setters */
        void setTriangular(bool);
        void setAU(bool);

        /* Bloch Hamiltonian */
        arma::cx_mat hamiltonian(arma::rowvec k) const;
        arma::cx_mat overlap(arma::rowvec k) const;
        using System::solveBands;
        void solveBands(arma::rowvec&, arma::vec&, arma::cx_mat&) const;

        /* Gauge */
        arma::cx_vec latticeToAtomicGauge(const arma::cx_vec&, const arma::rowvec&);
        arma::cx_vec atomicToLatticeGauge(const arma::cx_vec&, const arma::rowvec&);

        /* Observables */
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);      

        arma::cx_vec velocity(const arma::rowvec, int, int) const;  
    
    private:
        void orthogonalize_hamiltonian(const arma::rowvec&, arma::cx_mat&) const;
};

}
