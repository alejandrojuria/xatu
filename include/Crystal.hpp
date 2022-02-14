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
    protected:
        int ndim_, natoms_, ncells;
        double a, c;
        arma::mat bravaisLattice_, motif_, unitCellList_;
        arma::mat reciprocalLattice_;
        std::map<std::string, int> atomToIndex;

    // Const references to attributes (read-only)
    public:
        const int& ndim = ndim_;
        const int& natoms = natoms_;
        const arma::mat& bravaisLattice = bravaisLattice_;
        const arma::mat& motif = motif_;
        const arma::mat& unitCellList = unitCellList_;
        const arma::mat& reciprocalLattice = reciprocalLattice_;

    //// Methods
    protected:
        Crystal(){}; // Protected so that Crystal can not be initialized (abstract)
    public:   
        ~Crystal(){};

        int getDimension();
        int getNumAtoms();
        arma::mat getBravaisLattice();
        arma::mat getMotif();
        arma::mat getUnitCellList();
        arma::mat getReciprocalLattice();

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
