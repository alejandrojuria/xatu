#pragma once
#include <armadillo>
#include <complex>
#include <omp.h>
#include <stdlib.h>
#include <memory>

#include "BiRibbon.hpp"
#include "System.hpp"
#include "ExcitonConfiguration.hpp"
#include "Result.hpp"
#include "forward_declaration.hpp"
#include "utils.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


class GExciton : public System {

    // ----------------------------------- Attributes -----------------------------------
    private:
        // Read-only parameters
        int ncell_, totalCells_, nbands_, nrmbands_, excitonbasisdim_;
        double eps_m_, eps_s_, r0_;
        arma::ivec bands_, valenceBands_, conductionBands_;
        arma::uvec bandList_;
        arma::imat basisStates_;
        arma::rowvec Q_;
        double cutoff_;
        arma::cx_mat HBS_;

        // Flags
        std::string gauge_ = "lattice";
        std::string mode_  = "realspace";

        // Internal attributes
        arma::mat eigvalKStack_, eigvalKQStack_;
        arma::cx_cube eigvecKStack_, eigvecKQStack_;
        arma::cx_mat ftStack;
        arma::cx_cube ftMotifStack;
        std::complex<double> ftX;
        arma::mat potentialMat;
        arma::mat HK_;
        int nReciprocalVectors_ = 1;
        
        double pairEnergy;
        

    public:
        // Number of unit cells along one axis
        const int& ncell = ncell_;
        // Total number of unit cells
        const int& totalCells = totalCells_;
        // Dimension of electron-hole pair basis used to build excitons
        const int& excitonbasisdim = excitonbasisdim_;
        // Number of bands participating in exciton formation, starting from the Fermi level
        const int& nbands = nbands_;
        // Remove bands starting from the Fermi level to build excited excitons directly
        const int& nrmbands = nrmbands_;
        // List of bands used to build the exciton
        const arma::ivec& bands = bands_;
        // List of bands used to build the exciton relative to the Fermi level
        const arma::uvec& bandList = bandList_;
        // 3d array with the center-of-mass momentum of the exciton
        const arma::rowvec& Q = Q_;
        // List of valence bands that form the exciton relative to the Fermi level
        const arma::ivec& valenceBands = valenceBands_;
        // List of conduction bands that form the exciton relative to the Fermi level
        const arma::ivec& conductionBands = conductionBands_;
        // Returns Bethe-Salpeter Hamiltonian
        const arma::cx_mat& HBS = HBS_;
        // Returns kinetic term of BSE
        const arma::mat& HK = HK_;
        // Returns dielectric constant of embedding medium
        const double& eps_m = eps_m_;
        // Returns dielectric constante of substrate
        const double& eps_s = eps_s_;
        // Returns effective screening length r0
        const double& r0 = r0_;
        // Returns cutoff for potential
        const double& cutoff = cutoff_;
        // Returns gauge for Bloch states
        const std::string gauge = gauge_;
        // Return type of interaction matrix elements
        const std::string& mode = mode_;
        // Return number of reciprocal lattice vectors to use in summations (mode="reciprocalspace")
        const int& nReciprocalVectors = nReciprocalVectors_;

        const arma::mat& eigvalKStack = eigvalKStack_;
        const arma::mat& eigvalKQStack = eigvalKQStack_;
        const arma::cx_cube& eigvecKStack = eigvecKStack_;
        const arma::cx_cube& eigvecKQStack = eigvecKQStack_;
        const arma::imat& basisStates = basisStates_;

        // BEWARE: This dictionary had to be exposed to be able to access it,
        // do not call outside of class methods.
        std::map<int, int> bandToIndex;

    // ----------------------------------- Methods -----------------------------------
    // Constructor & Destructor
    private:
        // Remove default constructor
        GExciton();
        // Private constructor to leverage to another method parameter initialization
        GExciton(int, const arma::ivec&, const arma::rowvec&, const arma::rowvec&);

    public:
        // Specify number of bands participating (int)
        GExciton(std::string filename, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Specify which bands participate (vector with band numbers)
        GExciton(std::string filename, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});
        
        // Use two files: the mandatory one for system config., and one for exciton config.
        GExciton(std::string systemfile, std::string excitonfile);

        // Initialize exciton passing directly a System object instead of a file using removed bands
        GExciton(System&, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Initialize exciton passing directly a System object instead of a file using bands vector
        GExciton(System&, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        ~GExciton();

        // Setters
        void setUnitCells(int);
        void setBands(int, int);
        void setBands(const arma::ivec&);
        void setQ(const arma::rowvec&);
        void setParameters(const arma::rowvec&);
        void setParameters(double, double, double);
        void setCutoff(double);
        void setGauge(std::string);
        void setMode(std::string);
        void setReciprocalVectors(int);

    private:
        // Methods for BSE matrix initialization
        void STVH0(double, double*);
        double potential(double);
        std::complex<double> fourierTransform(arma::rowvec k, const arma::mat&, bool useApproximation = true);
        double analyticFourierTransform(arma::rowvec);
        double fourierTransformFromCoefs(const arma::vec&, const arma::vec&, const arma::rowvec&, int);
        std::complex<double> motifFourierTransform(int, int, const arma::rowvec&, const arma::mat&);
        arma::cx_mat extendMotifFT(const arma::cx_mat&);
        std::complex<double> blochCoherenceFactor(const arma::cx_vec&, const arma::cx_vec&, 
                                                  const arma::rowvec&, const arma::rowvec&,
                                                  const arma::rowvec&);

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
        std::complex<double> exactInteractionTermMFT(const arma::cx_vec&, 
                                                const arma::cx_vec&,
                                                const arma::cx_vec&, 
                                                const arma::cx_vec&,
                                                const arma::cx_mat&);
        std::complex<double> interactionTermFT(const arma::cx_vec&, 
                                                const arma::cx_vec&,
                                                const arma::cx_vec&, 
                                                const arma::cx_vec&,
                                                const arma::rowvec&, 
                                                const arma::rowvec&,
                                                const arma::rowvec&, 
                                                const arma::rowvec&,
                                                int nrcells = 15);

        // Initializers
        void initializeExcitonAttributes(int, const arma::ivec&, const arma::rowvec&, const arma::rowvec&);
        void initializeExcitonAttributes(const ExcitonConfiguration&);
        void initializeBasis();
        void initializeResultsH0();
        void initializePotentialMatrix();
        void initializeMotifFT(int, const arma::mat&);
        
        // Utilities
        void generateBandDictionary();
        void createMesh();
        
        // Gauge fixing
        arma::cx_vec latticeToAtomicGauge(const arma::cx_vec&, const arma::rowvec&);
        arma::cx_vec atomicToLatticeGauge(const arma::cx_vec&, const arma::rowvec&);
        arma::cx_mat fixGlobalPhase(arma::cx_mat&);

        // Routines to compute Fermi Golden Rule
        arma::cx_vec ehPairCoefs(double, const arma::vec&, bool zone = true);

    public:
        arma::imat createBasis(const arma::ivec&, const arma::ivec&);
        arma::imat specifyBasisSubset(const arma::ivec& bands);
        void useSpinfulBasis();

        // Symmetries
        arma::mat C3ExcitonBasisRep();
        
        // BSE initialization and energies
        void initializeHamiltonian();
        void BShamiltonian(const arma::imat& basis = {});
        Result diagonalize();

        // Fermi golden rule        
        arma::cx_vec wavePacket(double, double);
        double pairDensityOfStates(const arma::ivec&, const arma::ivec&, double, double);
        double fermiGoldenRule(const arma::cx_vec&, double);
};


