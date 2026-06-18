#pragma once
#include <armadillo>
#include "xatu/ConfigurationBase.hpp"

namespace xatu {

/**
 * ScreeningConfiguration is a specialization of ConfigurationBase to parse
 * dielectric function configuration files. 
 */
class ScreeningConfiguration : public ConfigurationBase{

    struct configuration {
        // Simulation label
        std::string label;
        // Number of unit cells in the calculation
        int ncell;
        // Number of valence bands included calculation.
        int nvbands = 0;
        // Number of conduction bands included calculation.
        int ncbands = 0;
        // Momentum to compute the dielectric function.
        arma::rowvec q = {0.2, 0., 0.};
        // Reduction factor of the BZ mesh. Defaults to 1.
        int submeshFactor = 1;
        // Pair of reciprocal vectors (G,G') to compute the reciprocal space dielectric function at
        arma::ivec Gs = {0, 0};
        // Pair of motif vectors (i,j) to compute the real space dielectric function at
        arma::ivec ts = {0, 0};
        // Number of Bravais lattice vectors to use in real space calculation
        int nLatticeVectors = 0;
        // Gc cutoff for the reciprocal lattice vectors
        double Gcutoff = 10.0;
        // Regularization distance
        double regularization = 0.0;
        // Compute dielectric function or polarizability or none
        std::string function = "dielectric";
        // Bool that is true/false when the dielectric function is/is not used
        bool isscreeningset = "false";
        // Bool that is true/false if the model is spinful/spinless
        bool spin = "true";
        // Number of k points in each direction for the auxiliary coarser BZ mesh
        int ncell_aux = 17;
        // Thickness of the 2D system
        double d = 0.0;
        // Is system isotropic
        bool isotropic = false;
    };

    public:
        configuration screeningInfo;
    
    public:
        ScreeningConfiguration();
        ScreeningConfiguration(std::string);
    
    private:
        virtual void parseContent();
        void checkContentCoherence();
};

}