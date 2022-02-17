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
        int ndim_, natoms_, nk_, ncells;
        double a, c;
        arma::mat bravaisLattice_, motif_, unitCellList_;
        arma::mat reciprocalLattice_, kpoints_;
        std::map<std::string, int> atomToIndex;

    // Const references to attributes (read-only)
    public:
        // Returns system dimension
        const int& ndim = ndim_;
        // Returns number of atoms in the motif
        const int& natoms = natoms_;
        // Returns Bravais lattice vectors of system
        const arma::mat& bravaisLattice = bravaisLattice_;
        // Returns position of atoms in motif
        const arma::mat& motif = motif_;
        // Returns list of Bravais vectors whose atom connect with the unit cell at the origin
        const arma::mat& unitCellList = unitCellList_;
        // Returns reciprocal lattice vectors
        const arma::mat& reciprocalLattice = reciprocalLattice_;
        // Returns kpoints
        const arma::mat& kpoints = kpoints_;
        // Number of k points
        const int& nk = nk_;

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
        void brillouinZoneMesh(int);
        void reducedBrillouinZoneMesh(int, int);
        void shiftBZ(const arma::rowvec&);
        void preserveC3();
        arma::mat wignerSeitzSupercell(int);
        arma::mat truncateSupercell(int, double);
        arma::mat generateCombinations(int n, int ndim, bool centered = false);

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
