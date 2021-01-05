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
    
    public:
        //// Attributes
        int ndim, natoms, norbitals, ncells, basisdim;
        double a, c;
        arma::mat bravais_lattice, motif, unitCellList;
        arma::mat reciprocal_lattice;
        arma::cx_cube hamiltonianMatrices;
        arma::vec kpoints;

    //// Methods
    // Constructor and destructor
    public:
        System(std::string);    
        ~System();

        arma::mat brillouin_zone_mesh(int);
        arma::cx_mat hamiltonian(arma::rowvec k);

        /* Some utilities/extra information */
        arma::cx_mat inversionOperator(const arma::cx_vec&);

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
        arma::mat generate_combinations(int n, int ndim);


    private:
        void readConfigurationFile(std::string);
        void extractLatticeParameters();
        void calculate_reciprocal_lattice();
};
