#pragma once
#include <armadillo>
#include <string>
#include "System.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

class BiRibbon : public System{
    
    public:
        //// Attributes
        int N;
        double a, c;
        double Es, Ep, Vsss, Vsps, Vpps, Vppp;
        double lambda, zeeman, onsiteEdge;
        std::string zeeman_axis;
        arma::rowvec a1, a2, tau;
        arma::rowvec n1, n2, n3;
        arma::rowvec Gamma, K, M;
        arma::mat M0, M1, M2p, M2m;
        arma::cx_mat Mso, Mzeeman;
    

    //// Methods
    // Constructor and destructor
    public:
        BiRibbon(int N = 15, std::string zeeman_axis = "z");    
        ~BiRibbon();

    protected:
        /* Attribute initialization */
        void initializeConstants();

        /* Matrix routines for hamiltonian initialization */
        arma::mat matrixWithSpin(const arma::mat&);
        arma::mat tightbindingMatrix(const arma::rowvec&);
        void createMotif();
        void initializeBlockMatrices();
        void prepareHamiltonian();

    public:
        void setZeeman(double);
        void applyElectricField(double);
        void offsetEdges(double);

        /* Some utilities/extra information */
        arma::cx_mat inversionOperator(const arma::cx_vec&);

};

}
