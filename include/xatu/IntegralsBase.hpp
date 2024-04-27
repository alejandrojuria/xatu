#pragma once
#include <string>
#include <vector>
#include <array>
#include <iostream>
#include <math.h>
#include <algorithm> 
#include <map> 
#include <unordered_map> 
#include <omp.h>
#include "xatu/GTFConfiguration.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ANG2AU 1.889726        
#endif

namespace xatu {

/// @brief The IntegralsBase family of classes is designed to compute the integrals in the
/// DFT and auxiliary basis sets for a given metric. This one contains the parameters FAC12, FAC3 and g^{l,m}_{i,j,k}
/// which are needed in all integrals, in addition to the method for E^{i,i'}_{t} (i,i'<=4) which is needed in all but
// the 2-center Coulomb integrals.
class IntegralsBase {

    //// Attributes copied from GTFConfiguration
    protected:
        int ndim_, natoms_, nspecies_, ncells_;
        arma::mat motif_, RlistAU_;
        std::map<int,int> RlistOpposites_;
        std::vector<int> R_unpaired_;

        // Const references to attributes (read-only)
    public:
        // Lattice dimensionality (0,1,2,3).
        const int& ndim = ndim_;
        // Number of atoms in the motif (unit cell).
        const int& natoms = natoms_;
        // Number of atomic species in the motif (unit cell).
        const int& nspecies = nspecies_;
        // Number of lattice vectors in which the integrals will be computed.
        const int& ncells = ncells_;
        // Matrix containing the positions of the atoms of the motif.
        // by rows. Each row has the format {x,y,z; species}, where x,y,z are in Angstrom and
        // the last element is an index representing the atomic species.
        const arma::mat& motif = motif_;
        // (3 x ncells) matrix of Bravais vectors in atomic units, stored by columns.
        const arma::mat& RlistAU = RlistAU_;
        // Map of ncells entries, the n-th (n=0,..,ncells-1) of which is the index of -R_{n} (minus the n-th Bravais vector in Rlist).
        const std::map<int,int>& RlistOpposites = RlistOpposites_;
        // Vector with the indices of the unpaired Bravais vectors, i.e. those with their opposites not present in Rlist.
        const std::vector<int>& R_unpaired = R_unpaired_;
 
    //// Attributes from this class
    public:
        // Location where the integrals files generated in this family of classes are stored.
        std::string IntFiles_Dir; 
        // Matrix dimension (or number of orbitals, with fixed l,m) in the SCF and AUX basis sets.
        uint32_t dimMat_SCF, dimMat_AUX;
        // Maximum l orbital quantum number among all shells and among both basis sets.
        int maxL;
        // The first dimension spans the orbitals (size = dimMat). The entries in the second dimension are: 1) atomic 
        // species number, 2) shell number for the species, 3) l quantum number, 4) m quantum number, 5) number of 
        // contracted Gaussians in the orbital.
        std::vector<std::vector<int>> orbitals_info_int_SCF, orbitals_info_int_AUX;
        // The first dimension spans the orbitals (size = dimMat). The entries in the second dimension are:
        // 1-3) atomic coordinates in atomic units, 4-end) {alpha_1,d_1,..,alpha_n,d_n}, where alpha are exponents and d 
        // are contraction coefficients (n being the number of contracted Gaussians).
        std::vector<std::vector<double>> orbitals_info_real_SCF, orbitals_info_real_AUX;    
        // Orbital normalization prefactor. Each entry is an orbital (size = dimMat): FAC12[orb] = FAC1[l(shell)][m(orb)] * FAC2[shell]
        std::vector<double> FAC12_SCF, FAC12_AUX;
        // Orbital normalization prefactor. The first dimension spans the orbitals (size = dimMat). 
        // The second dimension spans the contracted Gaussians.
        std::vector<std::vector<double>> FAC3_SCF, FAC3_AUX;
        // Coefficients g^{l,m}_{i,j,k} in the expansion G_{l,m}=\sum_{i,j,k=0}^{l}g^{l,m}_{i,j,k} x^i y^j z^k exp(-ar^2).
        // The key values are in a bijection with the pair (l,m) as key(l,m) = l^2 + l + m. The entries in the returned vector 
        // are: 1) number of nonzero g^{l,m}_{i,j,k} for fixed (l,m), 2-end) {i,j,k,g} for each nonzero coefficient.
        std::unordered_map<int,std::vector<int>> g_coefs;
        // Inverse of the bijective index reduction s(i,j) = j + i(i+1)/2 employed in triangular matrices. The key values are 
        // s = [0, 0.5*dimMat_AUX(dimMat_AUX+1) ) because it is to be used to re-index the 2-centers matrices in the auxiliary basis.
        // The entries in the returned array are (i,j), which for lower triangular matrices is (row,column). 
        std::unordered_map<uint64_t,std::array<uint32_t,2>> triangInd_to_rowcol;

    //// Methods
    public:
        IntegralsBase(const GTFConfiguration&, const std::string& integrals_directory = "Integrals_Files/"); 
        IntegralsBase(const IntegralsBase&); // Copy constructor
        ~IntegralsBase(){};
    
    protected:
        void initializeBasesAttributes(const GTFConfiguration&);
        void buildOrbitalsInfo(const GTFConfiguration&);
        int factorial(int);
        int doubleFactorial(int);
        // Method to build the matrix of the normalization prefactor FAC1(m,l)->FAC1[l][m], used in FAC12fun.
        std::vector<std::vector<double>> FAC1fun(const int maxL);
        // Method to build the vector of the normalization prefactor FAC2(shell,l)->FAC2[shell], used in FAC12fun.
        std::vector<double> FAC2fun(const GTFConfiguration&, const bool basis_id);
        // Method to build the vector FAC12[orb] attribute.
        void FAC12fun(const GTFConfiguration&, const int maxL, const bool basis_id);
        // Method to build the vector of vectors FAC3[orb][gaussian] attribute. 
        void FAC3fun(const bool basis_id);
        // Method to build the unordered_map method g_coefs, which contains the g^{l,m}_{i,j,k} expansion coefficients.
        void gfun(const int maxL);
        // Method to build the unordered_map method triangInd_to_rowcol, which contains the inverse of the (bijective) function: s(i,j) = j + i(i+1)/2.
        void triangIndfun(const uint32_t dimMat_AUX);
        // Method to compute the E^{i,i'}_{t} coefficients in the Hermite Gaussian expansion 
        // G_{i}G_{i'}=\sum_{t=0}^{i+i'}E^{i,i'}_{t}\Lambda_{t}, for i,i'<=4. Returns the vector for 0 <= t <= (i+i'), in increasing t order.
        // The values of the pair (i,i') with i>=i' are in a bijection with the "index" argument as index(i,i')= i' + i(i+1)/2.
        // For i<i', it must be called as Efun(index(i',i),p,PB,PA).  
        arma::colvec Efun(const int index, const double p, const double PA, const double PB);
        
};

}