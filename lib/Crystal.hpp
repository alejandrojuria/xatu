#pragma once
#include <armadillo>
#include <string>
#include <iostream>
#include "SystemConfiguration.hpp"


#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class Crystal {
    
    //// Attributes
    public:
        int ndim, natoms, ncells;
        double a, c;
        arma::mat bravaisLattice, motif, unitCellList;
        arma::mat reciprocalLattice;
        arma::vec kpoints;
        std::map<std::string, int> atomToIndex;

    //// Methods
    protected:
        Crystal(){}; // Protected so that Crystal can not be initialized (abstract)
    public:   
        ~Crystal(){};

        /* Mesh generation routines */
        arma::mat brillouinZoneMesh(int);
        arma::mat c3BzMesh(int);
        arma::mat wignerSeitzSupercell(int);
        arma::mat truncateSupercell(int, double);
        arma::mat generateCombinations(int n, int ndim);
        arma::mat generateCombinationsGamma(int n, int ndim);

        /* Crystal operations */
        arma::cx_mat inversionOperator(const arma::cx_vec&);
        arma::rowvec rotateC3(const arma::rowvec&);


    protected:
        void initializeCrystalAttributes(const SystemConfiguration&);
        void extractLatticeParameters();
        void calculateReciprocalLattice();
        bool isInsideWsCell(const arma::rowvec&, const arma::mat&, 
                               const arma::rowvec&);
};
