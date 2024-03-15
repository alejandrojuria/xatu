#pragma once
#include <complex>
#include <omp.h>
#include <stdlib.h>
#include <memory>
#include "xatu/System.hpp"
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

class Exciton {

    // ----------------------------------- Attributes -----------------------------------
    // Read-only parameters
    protected:

        std::shared_ptr<System> system_;

        // General Exciton attributes
        int ncell_, totalCells_, nbands_, nrmbands_, excitonbasisdim_;
        double scissor_;
        arma::ivec bands_, valenceBands_, conductionBands_;
        arma::uvec bandList_;
        arma::imat basisStates_;
        arma::rowvec Q_;
        double cutoff_;
        arma::cx_mat HBS_;

        // Flags
        bool exchange = false;

        // Internals for BSE
        arma::mat eigvalKStack_, eigvalKQStack_;
        arma::cx_cube eigvecKStack_, eigvecKQStack_;
        arma::mat HK_;

    public:
        const std::shared_ptr<System> system = system_;
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
        // Returns cutoff for potential
        const double& cutoff = cutoff_;
        // Returns scissor cut value
        const double& scissor = scissor_;

        const arma::mat& eigvalKStack = eigvalKStack_;
        const arma::mat& eigvalKQStack = eigvalKQStack_;
        const arma::cx_cube& eigvecKStack = eigvecKStack_;
        const arma::cx_cube& eigvecKQStack = eigvecKQStack_;
        const arma::imat& basisStates = basisStates_;

        // BEWARE: This dictionary had to be exposed to be able to access it,
        // do not call outside of class methods.
        std::map<int, int> bandToIndex;

    // ----------------------------------- Methods -----------------------------------
    public:
        // Constructor & Destructor
        virtual ~Exciton(){};

        // Setters
        void setUnitCells(int);
        void setBands(int, int);
        void setBands(const arma::ivec&);
        void setQ(const arma::rowvec&);
        void setCutoff(double);
        void setScissor(double);
        void setExchange(bool);

    protected:
        void initializeBasis();
        
        // Utilities
        void generateBandDictionary();
        
        // Gauge fixing
        arma::cx_mat fixGlobalPhase(arma::cx_mat&);

    public:
        arma::imat createBasis(const arma::ivec&, const arma::ivec&);
        virtual void printInformation();
        
        // BSE initialization and energies
        virtual void initializeHamiltonian() = 0;
        virtual void BShamiltonian() = 0;
        virtual Result diagonalize(std::string method = "diag", int nstates = 8) = 0;
};

}