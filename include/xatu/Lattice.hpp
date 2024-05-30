#pragma once
#include <numeric>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include "xatu/ConfigurationSystem.hpp"

#ifndef constants
#define PI 3.141592653589793
#endif

namespace xatu {

/** 
 *  The Lattice class serves as a base for the family of classes that manipulates information about reciprocal and direct space. 
 */
class Lattice {
    
    protected:
        // Lattice dimensionality (0,1,2,3)
        int ndim_;
        // Basis of Bravais vectors (R1,R2,R3), stored by columns: (3,ndim)
        arma::mat Rbasis_;
        // Unit cell volume, in the same units as Rbasis
        double unitCellVolume_;
        // Basis of reciprocal vectors (G1,G2,G3), stored by columns: (3,ndim)
        arma::mat Gbasis_;
        
    public:  // Const references to attributes (read-only)
        const int& ndim = ndim_;
        const arma::mat& Rbasis = Rbasis_;
        const double& unitCellVolume = unitCellVolume_;
        const arma::mat& Gbasis = Gbasis_;

    protected:
        Lattice() = default;
    public:
        Lattice(const ConfigurationSystem&);
        Lattice(const Lattice&) = default; 
        virtual ~Lattice() = default;

    public: 
        // Method to generate a kronecker-like list of integer combinations, to be used with Bravais vectors
        arma::mat generateCombinations(const int32_t n1, const int32_t n2, const int32_t n3, const bool centered = false);
        arma::mat generateCombinations(const int32_t n1, const int32_t n2, const bool centered = false);
        arma::mat generateCombinations(const int32_t n1, const bool centered = false);
        // Compute a Monkhorst-Pack grid in the interval [0 G1)x...x[0 Gn_dim), and return the k-points by columns in arma::mat (3,nk)
        arma::mat gridMonkhorstPack(const int32_t shrink1, const int32_t shrink2, const int32_t shrink3, const bool containsGamma = true);
        arma::mat gridMonkhorstPack(const int32_t shrink1, const int32_t shrink2, const bool containsGamma = true);
        arma::mat gridMonkhorstPack(const int32_t shrink1, const bool containsGamma = true);
        // Method to create the matrix of the first nR (at least) 3-component Bravais vectors, stored by columns and ordered by ascending norm.
        // The number of returned vectors is at least nR because full stars are given. Rbasis is (3,ndim) and contains R1,R2,R3 by columns.
        // It basically substitutes Rlist for the integrals because they generally use more R-vectors than those contained in the .outp
        arma::mat generateRlist(const uint32_t nR, const std::string& IntegralType);
        // Returns a map where each entry is the index of the direct lattice vector in the input generated_Rlist (generalization of Rlist for
        // an arbitrary number of direct lattice vectors) opposite to the lattice vector whose index is the corresponding map's key. 
        std::map<uint32_t,uint32_t> generateRlistOpposite(const arma::mat& generated_Rlist);
        
    protected:
        // Method to compute the unit cell volume, area or length (depending on the lattice dimensionality)
        void computeUnitCellVolume();
        // Compute the reciprocal lattice vectors {G_1,..,G_ndim} and return them by columns in arma::mat (3,ndim)
        void calculateGbasis();

};

}
