#pragma once
#include <armadillo>
#include <string>
#include <iostream>

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
        arma::mat bravais_lattice, motif, unitCellList;
        arma::mat reciprocal_lattice;
        arma::vec kpoints;

    //// Methods
    public:
        /* Constructor and destructor */
        Crystal(std::string, std::string overlapFile = "", bool isReal = false);    
        ~Crystal();

        /* Mesh generation routines */
        arma::mat brillouin_zone_mesh(int);
        arma::mat C3_BZ_Mesh(int);
        arma::mat wigner_seitz_supercell(int);
        arma::mat truncate_supercell(int, double);
        arma::mat generate_combinations(int n, int ndim);
        arma::mat generate_combinations_gamma(int n, int ndim);

        /* Some utilities/extra information */
        arma::cx_mat inversionOperator(const arma::cx_vec&);
        arma::rowvec rotateC3(const arma::rowvec&);


    private:
        void readConfigurationFile(std::string, bool isReal = false);
        void readOverlapFile(std::string, bool isReal = false);
        void extractLatticeParameters();
        void calculateReciprocalLattice();
        bool isInsideWSCell(const arma::rowvec&, const arma::mat&, 
                               const arma::rowvec&);
};
