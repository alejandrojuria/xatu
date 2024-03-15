#pragma once
#include <complex>
#include <omp.h>
#include <stdlib.h>
#include <memory>
#include "xatu/SystemTB.hpp"
#include "xatu/ExcitonConfiguration.hpp"
#include "xatu/Result.hpp"
#include "xatu/forward_declaration.hpp"
#include "xatu/utils.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

class ExcitonTB : public Exciton {

    // ----------------------------------- Attributes -----------------------------------
    // Read-only parameters
    private:
        
        std::shared_ptr<SystemTB> system_;

        // Keldysh potential constants
        double eps_m_, eps_s_, r0_;

        // Flags
        std::string gauge_ = "lattice";
        std::string mode_  = "realspace";

        // Internals for BSE
        arma::cx_mat ftMotifQ;
        arma::cx_cube ftMotifStack;
        std::complex<double> ftX;
        arma::mat potentialMat;
        int nReciprocalVectors_ = 1;
        

    public:
        const std::shared_ptr<SystemTB> system = system_;
        // Returns dielectric constant of embedding medium
        const double& eps_m = eps_m_;
        // Returns dielectric constante of substrate
        const double& eps_s = eps_s_;
        // Returns effective screening length r0
        const double& r0 = r0_;
        // Returns gauge for Bloch states
        const std::string gauge = gauge_;
        // Return type of interaction matrix elements
        const std::string& mode = mode_;
        // Return number of reciprocal lattice vectors to use in summations (mode="reciprocalspace")
        const int& nReciprocalVectors = nReciprocalVectors_;

    // ----------------------------------- Methods -----------------------------------
    // Constructor & Destructor
    private:
        // Private constructor to leverage to another method parameter initialization
        ExcitonTB(int, const arma::ivec&, const arma::rowvec&, const arma::rowvec&);

    public:
        // Specify number of bands participating (int)
        ExcitonTB(const SystemConfiguration&, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Specify which bands participate (vector with band numbers)
        ExcitonTB(const SystemConfiguration&, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});
        
        // Use two files: the mandatory one for system config., and one for exciton config.
        ExcitonTB(const SystemConfiguration&, const ExcitonConfiguration&);

        // Initialize exciton passing directly a System object instead of a file using removed bands
        ExcitonTB(SystemTB&, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Initialize exciton passing directly a System object instead of a file using bands vector
        ExcitonTB(SystemTB&, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        ~ExcitonTB(){};

        // Setters
        void setParameters(const arma::rowvec&);
        void setParameters(double, double, double);
        void setGauge(std::string);
        void setMode(std::string);
        void setReciprocalVectors(int);

    private:
        // Methods for BSE matrix initialization
        void STVH0(double, double*);
        double potential(double);
        std::complex<double> fourierTransform(arma::rowvec k, const arma::mat&, bool useApproximation = true);
        double analyticFourierTransform(arma::rowvec);
        std::complex<double> motifFourierTransform(int, int, const arma::rowvec&, const arma::mat&);
        arma::cx_mat motifFTMatrix(const arma::rowvec&, const arma::mat&);
        arma::cx_mat extendMotifFT(const arma::cx_mat&);
        std::complex<double> blochCoherenceFactor(const arma::cx_vec&, const arma::cx_vec&, 
                                                  const arma::rowvec&, const arma::rowvec&,
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
        void initializeResultsH0(bool triangular = false);
        void initializeMotifFT(int, const arma::mat&);
        arma::imat specifyBasisSubset(const arma::ivec& bands);

        
        // Gauge fixing
        arma::cx_mat fixGlobalPhase(arma::cx_mat&);


    public:
        // BSE initialization and energies
        void initializeHamiltonian(bool triangular = false);
        void BShamiltonian(const arma::imat& basis = {});
        Result diagonalize(std::string method = "diag", int nstates = 8);

        // Fermi golden rule       
        double pairDensityOfStates(double, double) const;
        void writePairDOS(FILE*, double delta, int n = 100);
        arma::cx_vec ehPairCoefs(double, const arma::vec&, std::string side = "right");
        double fermiGoldenRule(const Exciton&, const arma::cx_vec&, const arma::cx_vec&, double);
        double edgeFermiGoldenRule(const Exciton&, const arma::cx_vec&, double, std::string side = "right", bool increasing = false);

        // Auxiliary routines for Fermi golden rule
        arma::rowvec findElectronHolePair(const Exciton&, double, std::string, bool);
};

}