#pragma once
#include <armadillo>
#include <complex>
#include <omp.h>
#include <stdlib.h>
#include "Zigzag.hpp"
#include "System.hpp"


#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


// Routines to initialize BSE elements
class GExciton : System {

    //// Attributes
    private:
        arma::mat eigvalKStack, eigvalKQStack;
        arma::cx_vec ftStack;
        std::map<int, int> bandToIndex;
        
        double pairEnergy;

    public:
        arma::mat potentialMat;
        int Ncell, nk, nbands, nrmbands, fermiLevel, basisdim, excitonbasisdim;
        double filling;
        arma::vec bands;
        arma::vec valenceBands, conductionBands;
        arma::uvec bandList;
        arma::rowvec Q;
        arma::cx_mat HBS;
        arma::mat HK;
        arma::mat basisStates;
        arma::mat kpoints;
        
        arma::cx_cube eigvecKStack, eigvecKQStack;

    //// Methods
    // Constructor & Destructor
    public:
        // Overload constructor:
        // Specify number of bands participating (int)
        GExciton(std::string filename, int Ncell = 200, 
                 const arma::rowvec& Q = {0., 0., 0.}, int nbands = 2, int nrmbands = 0, 
                 double filling = 0.5);
        // Specify which bands participate (vector with band numbers)
        GExciton(std::string filename, int Ncell = 200, 
                 double Q = 0, arma::vec bands = {}, double filling = 0.5);
        ~GExciton();

    private:
        // Methods for BSE matrix initialization
        void STVH0(double, double*);
        double potential(double);
        std::complex<double> fourierTrans(arma::rowvec k);
        std::complex<double> tDirect(std::complex<double>,
                                     const arma::cx_vec&, 
                                     const arma::cx_vec&,
                                     const arma::cx_vec&, 
                                     const arma::cx_vec&);
        std::complex<double> tExchange(std::complex<double>, 
                                    const arma::cx_vec&, 
                                    const arma::cx_vec&,
                                    const arma::cx_vec&, 
                                    const arma::cx_vec&);
        std::complex<double> exactInteractionTerm(const arma::cx_vec&, 
                                 const arma::cx_vec&,
                                 const arma::cx_vec&, 
                                 const arma::cx_vec&,
                                 const arma::rowvec&);

        // Initializers
        void initializeResultsH0();
        void initializePotentialMatrix();

        // Utilities
        void generateBandDictionary();
        void createMesh();
        void fixBandCrossing(arma::vec&, arma::cx_mat&);
        int determineKIndex(double k);
        arma::cx_cube atomicGCoefs(const arma::cx_cube&);

        // Routines to compute Fermi Golden Rule
        arma::cx_mat fixDegeneracyIteration(const arma::cx_vec&, const arma::cx_vec&);
        arma::cx_vec ehPairCoefs(double, const arma::vec&, bool zone = true);

    public:
        // Basis creation methods
        void initializeBasis();
        arma::mat createBasis(const arma::vec&, const arma::vec&);
        arma::mat specifyBasisSubset(const arma::vec& bands);
        void createSOCBasis();

        // BSE initialization and energies
        void BShamiltonian(const arma::mat& basis = {});
        arma::vec computeEnergies(const arma::cx_vec&);

        // Observables
        arma::cx_vec spinX(const arma::cx_vec&);
        
        arma::cx_vec wavePacket(double, double);
        arma::cx_mat fixDegeneracy(const arma::cx_vec&, const arma::cx_vec&, int iterations = 5);
        double pairDensityOfStates(const arma::vec&, const arma::vec&, double, double);
        double fermiGoldenRule(const arma::cx_vec&, double);  
};


