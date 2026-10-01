#pragma once
#include <complex>
#include <omp.h>
#include <stdlib.h>
#include "xatu/SystemTB.hpp"
#include "xatu/Exciton.hpp"
#include "xatu/ExcitonConfiguration.hpp"
#include "xatu/ScreeningConfiguration.hpp"
#include "xatu/ResultTB.hpp"
#include <array>
#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

/** 
 *Method to prompt user to continue
 *@param message
 *@return void 
*/
void continueprompt(std::string);

namespace xatu {

// Definitions

// Define pointer to potential function type (f: R -> R) within those in ExcitonTB
typedef double (ExcitonTB::*potptr)(arma::rowvec);
class ResultTB;

class ExcitonTB : public Exciton<SystemTB> {

    // ----------------------------------- Attributes -----------------------------------
    // Read-only parameters
    private:
        
        // Keldysh potential constants
        double eps_m_, eps_s_, r0_, ry_, rz_;
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
        double Gc_exciton_ = 10.0;
        int nReciprocalVectors_ = 1;
        
        // Internals for dielectric function
        std::string function_;
        arma::ivec Gs_ = {0, 0};
        arma::rowvec q_ = {0.2, 0., 0.};
        uint nGs;
        double Gcutoff_ = 3.0;
        int Nqpoints = 1;   // Number of q points to compute the dielectric function at        
        double d_ = 0.0; // Effective thickness of the 2D system
        uint ncell_aux_ = 10; // Number of unit cells used to compute the dielectric function
        uint nk_aux_ = ncell_aux_*ncell_aux_; // Number of k points used to compute the dielectric function
        uint g_s_ = 1; // Spin degeneracy
        arma::mat kpoints_aux_; // Auxiliar coarser BZ mesh to compute the dielectric function
        arma::mat qpoints_list_;
        arma::mat trunreciprocalLattice_;

        // Internals for dielectric fucntion
        arma::mat eigvalkStack_, eigvalkqStack_;
        arma::cx_cube eigveckStack_, eigveckqStack_;

        arma::cx_cube Chimatrix_;
        arma::cx_cube epsilonmatrix_;
        arma::cx_cube Invepsilonmatrix_;

        double percentage_ = 0; // Percentage of norm of k0, where k0 is the k point in the BZ mesh closer to the origin, for regularization
        double W00_at_0_ = 0.0;
        arma::rowvec q0_ = {0.2, 0.0, 0.0};
        bool isotropic_ = false; // Is system isotropic

    public:
        // Returns dielectric constant of embedding medium
        const double& eps_m = eps_m_;
        // Returns dielectric constante of substrate
        const double& eps_s = eps_s_;
        // Returns effective screening length r0
        const double& r0 = r0_;
        // Returns effective screening length ry
        const double& ry = ry_;
        // Returns effective screening length rz
        const double& rz = rz_;
        // Returns regularization distance
        const double& regularization = regularization_;
        // Returns gauge for Bloch states
        const std::string gauge = gauge_;
        // Return type of interaction matrix elements
        const std::string& mode = mode_;
        // Return number of reciprocal lattice vectors to use in summations (mode="reciprocalspace")
        const int& nReciprocalVectors = nReciprocalVectors_;
        // Return Gc cutoff for reciprocal lattice vectors to use in summations (mode="reciprocalspace")
        const double& Gc_exciton = Gc_exciton_;
        // Returns auxiliar coarser BZ mesh to compute the dielectric function
        arma::mat& kpoints_aux = kpoints_aux_;
        // Returns percentage of norm of k0, where k0 is the k point in the BZ mesh closer to the origin, for regularization
        double& percentage = percentage_;
        // Returns whether system is isotropic
        bool& isotropic = isotropic_;

        // Returns spin degeneracy
        uint& g_s = g_s_;
        // Returns thickness
        double& d = d_;
        
        // Return momentum to compute the dielectric matrix at
        const arma::rowvec& q = q_;
        const arma::ivec& Gs = Gs_;

        // Return number of cells in each direction used to compute the dielectric function
        const uint& ncell_aux = ncell_aux_;
        // Return number of k points used to compute the dielectric function
        const uint& nk_aux = nk_aux_;
        
        arma::mat trunLattice_;
    // ----------------------------------- Methods -----------------------------------
    // Constructor & Destructor
    private:
        // Private constructor to leverage to another method parameter initialization
        ExcitonTB(int, const arma::ivec&, const arma::rowvec&, const arma::rowvec&);

    public:
        // Specify number of bands participating (int)
        ExcitonTB(const SystemConfiguration&, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1, 1, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Specify which bands participate (vector with band numbers)
        ExcitonTB(const SystemConfiguration&, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1, 1, 1}, const arma::rowvec& Q = {0., 0., 0.});
        
        // Use two files: the mandatory one for system config., and one for exciton config.
        ExcitonTB(const SystemConfiguration&, const ExcitonConfiguration&);

        // Use three files: the mandatory one for system config., one for the screening config. and one for exciton config.
        ExcitonTB(const SystemConfiguration&, const ExcitonConfiguration&, const ScreeningConfiguration&);

        // Initialize exciton passing directly a System object instead of a file using removed bands
        ExcitonTB(std::shared_ptr<SystemTB>, int ncell = 20, int nbands = 1, int nrmbands = 0, 
                 const arma::rowvec& parameters = {1, 5, 1, 1, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Initialize exciton passing directly a System object instead of a file using bands vector
        ExcitonTB(std::shared_ptr<SystemTB>, int ncell = 20, const arma::ivec& bands = {0, 1}, 
                 const arma::rowvec& parameters = {1, 5, 1, 1, 1}, const arma::rowvec& Q = {0., 0., 0.});

        // Initialize exciton passing directly a System object instead of a screening
        ExcitonTB(const SystemConfiguration&, int ncell, int nbands, int nrmbands,
                  const arma::rowvec& parameters, const arma::rowvec& Q, const uint ncell_aux, const uint nvbands, const uint ncbands, const double Gcutoff, const double Gc_exciton, const double d, const bool spin = true, const bool isotropic = false);

        // ~ExcitonTB();

        // Setters
        void setParameters(const arma::rowvec&);
        void setParameters(double, double, double, double, double);
        void setGauge(std::string);
        void setMode(std::string);
        void setRegularization(double);
        void setGcutoff(double);
        void setValenceBands(int);
        void setConductionBands(int);
        void setVectors(int,int);
        void setVectors(arma::ivec);
        void setPotential(std::string);
        void setExchangePotential(std::string);
        void setTrunLattice(int,double);
        void setq_points_list(arma::mat);
        void setPercentage(double);
        void setThickness(double);

        // Getters
        int getNGs() const;
        double getThickness() const;
        arma::cx_cube getChiMatrix() const;
        arma::cx_cube getInverseDielectricMatrix() const;
        std::vector<std::vector<std::vector<std::complex<double>>>> getPolarizabilityMatrix_as_vector() const;
        std::vector<std::vector<std::vector<std::complex<double>>>> getInverseDielectricMatrix_as_vector() const;

    private:
        // Potentials
        double keldysh(arma::rowvec);
        void STVH0(double, double*);
        double coulomb(arma::rowvec);
        potptr selectPotential(std::string);

        // Fourier transforms
        double coulomb_2D_FT(const arma::rowvec&) const;
        double keldyshFT(arma::rowvec) const;
        std::complex<double> rpaFT(int g, int g2, arma::rowvec) const;
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
                                                       const arma::rowvec&, const arma::rowvec&, 
                                                       std::string, int nrcells = 15);
        std::complex<double> blochCoherenceFactor(const arma::cx_vec&, const arma::cx_vec&, 
                                                  const arma::rowvec&, const arma::rowvec&,
                                                  const arma::rowvec&);
        std::complex<double> blochCoherenceFactor(const arma::cx_vec &, const arma::cx_vec &,
                                                  const arma::rowvec &, const arma::rowvec &,
                                                  const arma::rowvec &, const double) const;

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

        // Static dielectric function
        int fetchReciprocalLatticeVector(arma::rowvec);
        std::complex<double> computesinglePolarizabilityMatrixElement(arma::rowvec &, arma::rowvec &, arma::rowvec &);
        std::complex<double> compute_2D_PolarizabilityMatrixElement(const arma::rowvec&, const arma::rowvec&, const arma::rowvec&);        
        std::complex<double> compute_quasi2D_PolarizabilityMatrixElement(const arma::rowvec&, const arma::rowvec&, const arma::rowvec&, const double);
        std::complex<double> compute_2D_PolarizabilityMatrixElement(const arma::rowvec&, const arma::rowvec&, const arma::rowvec&, const arma::rowvec&, const double);        
        std::complex<double> compute_2D_DielectricMatrixElement(const arma::rowvec&, const arma::rowvec&, const arma::rowvec&);
        std::complex<double> compute_quasi2D_DielectricMatrixElement(const arma::rowvec &G, const arma::rowvec &G2, const arma::rowvec &q, const double);

        public:
        // Static dielectric function, BSE initialization and energies
        void initializeHamiltonian();
        void initializeScreeningAttributes(const ScreeningConfiguration&);
        void initializeScreeningAttributes(const ScreeningConfiguration&, const std::string);
        std::complex<double> computesingleDielectricFunctionMatrixElement();
        void computesingleInverseDielectricMatrix(std::string);
        void PolarizabilityMesh();
        void compute_2D_DielectricMatrix();
        void compute_2D_DielectricMatrix(std::string);
        void augment_2D_DielectricMatrix(double);
        void compute_ScreenedPotential_regularization(bool);
        void compute_2D_PolarizabilityMatrix(std::string);
        void invertDielectricMatrix();
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
        void printReciprocalLattice();
     
        // Write BZ mesh in a file
        void writeBZtofile() const;
        // Write polarizability matrix in a file
        void writePolarizabilityMatrix(std::string) const;
        // Read inverse of dielectric matrix in a file
        void readInverseDielectricMatrix(std::string);
        // Write inverse of dielectric matrix in a file
        void writeInverseDielectricMatrix(std::string); // if the user did not compute the inverse dielectric matrix, then the method tries to do it
        // Write dielectric matrix in a file
        void writeDielectricMatrix(std::string) const;

        // Verifies if potential chosen is 'rpa' and if a screening file was not provided
        void verifypotential();
};

}