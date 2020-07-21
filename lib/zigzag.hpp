#pragma once
#include <armadillo>
#include <string>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class Zigzag{
    
    public:
        //// Attributes
        int N;
        double a, c;
        double Es, Ep, Vsss, Vsps, Vpps, Vppp;
        double lambda, zeeman, onsiteEdge;
        std::string zeeman_axis;
        arma::vec a1, a2, tau;
        arma::vec n1, n2, n3;
        arma::mat motif;
        arma::vec Gamma, K, M;
        arma::mat M0, M1, M2p, M2m;
        arma::cx_mat Mso, Mzeeman;
        arma::cx_mat H0, Ha, Hsoc, Hzeeman;
        arma::vec kpoints;

    //// Methods
    // Constructor and destructor
    public:
        Zigzag(int N = 15, std::string zeeman_axis = "z");    
        ~Zigzag();

    protected:
        /* Attribute initialization */
        void initializeConstants();

        /* Matrix routines for hamiltonian initialization */
        arma::mat matrixWithSpin(const arma::mat&);
        arma::mat tightbindingMatrix(const arma::vec&);
        void createMotif();
        void initializeBlockMatrices();
        void prepareHamiltonian();

        /* Routines for DoS calculation */
        std::complex<double> rGreenF(double, double, double);

    public:
        void setZeeman(double);
        arma::cx_mat hamiltonian(double k);

        /* Some utilities/extra information */
        void writeEigenvaluesToFile(FILE* file, const arma::vec& eigenval, double k);
        arma::cx_mat inversionOperator(const arma::cx_vec&);

        /* Expected value of spin components */
        double expectedSpinZValue(const arma::cx_vec&);
        double expectedSpinYValue(const arma::cx_vec&);
        double expectedSpinXValue(const arma::cx_vec&);

        /* Routines for DoS calculation */
        double densityOfStates(double, double, const arma::mat&);
        void writeDensityOfStates(const arma::mat&, double, FILE*);
};
