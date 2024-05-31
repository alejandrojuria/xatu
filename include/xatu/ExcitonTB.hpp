#pragma once
#include "xatu/ConfigurationExciton.hpp"
#include "xatu/SystemTB.hpp"
#include "xatu/Exciton.hpp"
#include "xatu/ResultTB.hpp"
#include "xatu/utils.hpp"
#include "xatu/davidson.hpp"

#ifndef constants
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

// Definitions

// Define pointer to potential function type (f: R -> R) within those in ExcitonTB
typedef double (ExcitonTB::*potptr)(double);
class ResultTB;

class ExcitonTB : public Exciton<SystemTB> {

    // ----------------------------------- Attributes -----------------------------------
    // Read-only parameters
    private:
        
        // Keldysh potential constants
        double eps_m_, eps_s_, r0_;
        double regularization_;

        // Flags
        std::string gauge_ = "lattice";
        std::string mode_  = "realspace";
        std::string potential_ = "keldysh";
        std::string exchangePotential_ = "keldysh";

        // Internals for BSE
        arma::cx_mat ftMotifQ;
        arma::cx_cube ftMotifStack;
        std::complex<double> ftX;
        arma::mat potentialMat;
        int nReciprocalVectors_ = 1;
        

    public:
        // Returns dielectric constant of embedding medium
        const double& eps_m = eps_m_;
        // Returns dielectric constante of substrate
        const double& eps_s = eps_s_;
        // Returns effective screening length r0
        const double& r0 = r0_;
        // Returns regularization distance
        const double& regularization = regularization_;
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
        ExcitonTB(std::shared_ptr<SystemTB>, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Initialize exciton passing directly a System object instead of a file using bands vector
        ExcitonTB(std::shared_ptr<SystemTB>, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1}, const arma::rowvec& Q = {0., 0., 0.});

        ~ExcitonTB();

        // Setters
        void setParameters(const arma::rowvec&);
        void setParameters(double, double, double);
        void setGauge(std::string);
        void setMode(std::string);
        void setReciprocalVectors(int);
        void setRegularization(double);

    private:
        // Potentials
        double keldysh(double);
        void STVH0(double, double*);
        double coulomb(double);
        potptr selectPotential(std::string);

        // Fourier transforms
        double keldyshFT(arma::rowvec);
        std::complex<double> motifFourierTransform(int, int, const arma::rowvec&, const arma::mat&, potptr);
        arma::cx_mat motifFTMatrix(const arma::rowvec&, const arma::mat&, potptr);
        arma::cx_mat extendMotifFT(const arma::cx_mat&);

        // Interaction matrix elements
        std::complex<double> realSpaceInteractionTerm(const arma::cx_vec&, const arma::cx_vec&,
                                                      const arma::cx_vec&, const arma::cx_vec&,
                                                      const arma::cx_mat&);
        std::complex<double> reciprocalInteractionTerm(const arma::cx_vec&, const arma::cx_vec&,
                                                       const arma::cx_vec&, const arma::cx_vec&,
                                                       const arma::rowvec&, const arma::rowvec&,
                                                       const arma::rowvec&, const arma::rowvec&, int nrcells = 15);
        std::complex<double> blochCoherenceFactor(const arma::cx_vec&, const arma::cx_vec&, 
                                                  const arma::rowvec&, const arma::rowvec&,
                                                  const arma::rowvec&);

        // Initializers
        void initializeExcitonAttributes(int, const arma::ivec&, const arma::rowvec&, const arma::rowvec&);
        void initializeExcitonAttributes(const ExcitonConfiguration&);
        void initializeResultsH0();
        void initializeMotifFT(int, const arma::mat&, potptr);
        arma::imat specifyBasisSubset(const arma::ivec& bands);

        // Gauge fixing
        arma::cx_mat fixGlobalPhase(arma::cx_mat&);

        // Diagonalization
        ResultTB* diagonalizeRaw(std::string method = "diag", int nstates = 8) override;

    public:
        // BSE initialization and energies
        void initializeHamiltonian();
        void BShamiltonian();
        void BShamiltonian(const arma::imat& basis);
        std::unique_ptr<ResultTB> diagonalize(std::string method = "diag", int nstates = 8);
        
        // Fermi golden rule       
        double pairDensityOfStates(double, double) const;
        void writePairDOS(FILE*, double delta, int n = 100);
        arma::cx_vec ehPairCoefs(double, const arma::vec&, std::string side = "right");
        double fermiGoldenRule(const ExcitonTB&, const arma::cx_vec&, const arma::cx_vec&, double);
        double edgeFermiGoldenRule(const ExcitonTB&, const arma::cx_vec&, double, std::string side = "right", bool increasing = false);

        // Auxiliary routines for Fermi golden rule
        arma::rowvec findElectronHolePair(const ExcitonTB&, double, std::string, bool);

        // Print information
        void printInformation();
};

}