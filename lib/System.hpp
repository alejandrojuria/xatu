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
        int ndim, nmotif, norbitals, ncells;
        arma::mat bravais_lattice, motif, unitCellList;
        arma::mat reciprocal_lattice;
        arma::cx_cube hamiltonianMatrices;
        arma::vec kpoints;

    //// Methods
    // Constructor and destructor
    public:
        System(std::string);    
        ~System();

    protected:

        /* Routines for DoS calculation */
        std::complex<double> rGreenF(double, double, double);

    public:
        arma::cx_mat hamiltonian(arma::vec k);

        /* Some utilities/extra information */
        arma::cx_mat inversionOperator(const arma::cx_vec&);

        /* Expected value of spin components */
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);

        /* Routines for DoS calculation */
        double densityOfStates(double, double, const arma::mat&);
        void writeDensityOfStates(const arma::mat&, double, FILE*);

    private:
        void readConfigurationFile(std::string);
        arma::mat calculate_reciprocal_lattice();
        arma::mat brillouin_zone_mesh(int);
        arma::mat generate_combinations(int, int);
};
