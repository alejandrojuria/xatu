#pragma once
#include <armadillo>
#include <string>
#include <iostream>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class System {
    
    //// Attributes
    public:
        int ndim, natoms, norbitals, ncells, basisdim;
        double a, c;
        arma::mat bravais_lattice, motif, unitCellList;
        arma::mat reciprocal_lattice;
        arma::cx_cube hamiltonianMatrices;
        arma::cx_cube overlapMatrices;
        arma::vec kpoints;

    //// Methods
    public:
        /* Constructor and destructor */
        System(std::string, std::string overlapFile = "", bool isReal = false);    
        ~System();

        /* Mesh generation routines */
        arma::mat brillouin_zone_mesh(int);
        arma::mat C3_BZ_Mesh(int);
        arma::mat wigner_seitz_supercell(int);
        arma::mat truncate_supercell(int, double);
        arma::mat generate_combinations(int n, int ndim);
        arma::mat generate_combinations_gamma(int n, int ndim);

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

        /* Some utilities/extra information */
        arma::cx_mat inversionOperator(const arma::cx_vec&);
        arma::rowvec rotateC3(const arma::rowvec&);

    protected:
        /* Routines for DoS calculation */
        std::complex<double> rGreenF(double, double, double);


    private:
        void readConfigurationFile(std::string, bool isReal = false);
        void readOverlapFile(std::string, bool isReal = false);
        void extractLatticeParameters();
        void calculate_reciprocal_lattice();
        bool is_inside_ws_cell(const arma::rowvec&, const arma::mat&, 
                               const arma::rowvec&);
};
