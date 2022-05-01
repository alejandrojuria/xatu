#pragma once
#include <armadillo>
#include <complex>
#include <omp.h>
#include <stdlib.h>

#include "BiRibbon.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


// Routines to initialize BSE elements
class Exciton : public BiRibbon{

    //// Attributes
    private:
        int basisDimTB, nk;
        arma::vec specifyEdges;
        arma::uvec bandList;
        arma::mat eigvalKStack, eigvalKQStack;
        arma::cx_vec ftStack;
        std::map<int, int> bandToIndex;
        arma::uvec edgeV, edgeC;
        
        double pairEnergy;

    public:
        arma::mat potentialMat;
        int Ncell;
        double Q;
        int nBulkBands, nEdgeBands;
        arma::cx_mat HBS;
        arma::mat HK;
        arma::mat basisStates;
        arma::mat kpoints;
        arma::uvec valenceBands, conductionBands;
        arma::cx_cube eigvecKStack, eigvecKQStack;

    //// Methods
    // Constructor & Destructor
    public:
        Exciton(int N = 15, int Ncell = 200, double Q = 0, int nBulkBands = 2, int nEdgeBands = 0, arma::vec specifyEdges = {});
        ~Exciton();

    private:
        // Methods for BSE matrix initialization
        void STVH0(double, double*);
        double potential(double);
        std::complex<double> fourierTrans(double);
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
                                 const arma::vec&);

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
        arma::mat createBasis(int nBulkBands, int nEdgeBands);
        void initializeBasis();
        void createSOCBasis();

        // BSE initialization and energies
        void BShamiltonian();
        arma::vec computeEnergies(const arma::cx_vec&);

        // Observables
        arma::cx_vec spinX(const arma::cx_vec&);
        
        arma::cx_vec wavePacket(double, double);
        arma::cx_mat fixDegeneracy(const arma::cx_vec&, const arma::cx_vec&, int iterations = 5);
        double pairDensityOfStates(double, double);
        double fermiGoldenRule(const arma::cx_vec&, double);  
};


