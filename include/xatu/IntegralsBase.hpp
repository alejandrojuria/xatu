#pragma once
#include <array>
#include <unordered_map> 
#include "xatu/ConfigurationCRYSTAL.hpp"
#include "xatu/ConfigurationGTF.hpp"
#include "xatu/Lattice.hpp"

#ifndef constants
#define ANG2AU 1.8897261     
#endif

namespace xatu {

/** 
 *  The IntegralsBase class serves as a base for the family of classes that is designed to compute the integrals in the DFT and 
 *  AUXILIARY basis sets for a given metric. This class contains (as attributes) the parameters FAC12, FAC3 and g^{l,m}_{i,j,k}
 *  which are needed in all integrals, in addition to the method for E^{i,i'}_{t} (i,i'<=4) which is needed in all but
 *  the 2-center Coulomb integrals. Exclusive to the Gaussian mode
 */
class IntegralsBase : public Lattice {

    protected:
        // Location where the integrals and list of Bravais vectors files generated in this family of classes are stored
        std::string IntFiles_Dir_; 
        // Matrix dimension (or number of orbitals, with fixed l,m) in the SCF and AUX basis sets
        uint32_t dimMat_SCF_, dimMat_AUX_;
        // Maximum l orbital quantum number among all shells and among both basis sets
        int maxL_;
        // The first dimension spans the orbitals (size = dimMat). The entries in the second dimension are: 1) atomic 
        // species number, 2) shell number for the species, 3) l quantum number, 4) m quantum number, 5) number of 
        // contracted Gaussians in the orbital
        std::vector<std::vector<int>> orbitals_info_int_SCF_, orbitals_info_int_AUX_;
        // The first dimension spans the orbitals (size = dimMat). The entries in the second dimension are:
        // 1-3) atomic coordinates in atomic units, 4-end) {alpha_1,d_1,..,alpha_n,d_n}, where alpha are exponents and d 
        // are contraction coefficients (n being the number of contracted Gaussians)
        std::vector<std::vector<double>> orbitals_info_real_SCF_, orbitals_info_real_AUX_;    
        // Orbital normalization prefactor. Each entry is an orbital (size = dimMat): FAC12[orb] = FAC1[l(shell)][m(orb)] * FAC2[shell]
        std::vector<double> FAC12_SCF_, FAC12_AUX_;
        // Orbital normalization prefactor. The first dimension spans the orbitals (size = dimMat). 
        // The second dimension spans the contracted Gaussians
        std::vector<std::vector<double>> FAC3_SCF_, FAC3_AUX_;
        // Coefficients g^{l,m}_{i,j,k} in the expansion G_{l,m}=\sum_{i,j,k=0}^{l}g^{l,m}_{i,j,k} x^i y^j z^k exp(-ar^2).
        // The key values are in a bijection with the pair (l,m) as key(l,m) = l^2 + l + m. The entries in the returned vector 
        // are: 1) number of nonzero g^{l,m}_{i,j,k} for fixed (l,m), 2-end) {i,j,k,g} for each nonzero coefficient
        std::unordered_map<int,std::vector<int>> g_coefs_;
        // Inverse of the bijective index reduction s(i,j) = j + i(i+1)/2 employed in triangular matrices. The key values are 
        // s = [0, 0.5*dimMat_AUX(dimMat_AUX+1) ) because it is to be used to re-index the 2-centers matrices in the auxiliary basis.
        // The entries in the returned array are (i,j), which for lower triangular matrices is (row,column) 
        std::unordered_map<uint64_t,std::array<uint32_t,2>> triangInd_to_rowcol_;

    public:  // Const references to attributes (read-only)
        const std::string& IntFiles_Dir = IntFiles_Dir_;
        const uint32_t& dimMat_SCF = dimMat_SCF_;
        const uint32_t& dimMat_AUX = dimMat_AUX_;
        const int& maxL = maxL_;
        const std::vector<std::vector<int>>& orbitals_info_int_SCF = orbitals_info_int_SCF_;
        const std::vector<std::vector<int>>& orbitals_info_int_AUX = orbitals_info_int_AUX_;
        const std::vector<std::vector<double>>& orbitals_info_real_SCF = orbitals_info_real_SCF_;
        const std::vector<std::vector<double>>& orbitals_info_real_AUX = orbitals_info_real_AUX_;
        const std::vector<double>& FAC12_SCF = FAC12_SCF_;
        const std::vector<double>& FAC12_AUX = FAC12_AUX_;
        const std::vector<std::vector<double>>& FAC3_SCF = FAC3_SCF_;
        const std::vector<std::vector<double>>& FAC3_AUX = FAC3_AUX_;
        const std::unordered_map<int,std::vector<int>>& g_coefs = g_coefs_;
        const std::unordered_map<uint64_t,std::array<uint32_t,2>>& triangInd_to_rowcol = triangInd_to_rowcol_;

    protected:
        IntegralsBase();
        IntegralsBase(const IntegralsBase&) = default; 
    public:
        IntegralsBase(const ConfigurationCRYSTAL&, const ConfigurationGTF&, const std::string& integrals_directory = "Integrals_Files/"); 
        virtual ~IntegralsBase() = default;
    
    private:
        // Method to build the arma::mat orbitals_info attributes
        void buildOrbitalsInfo(const ConfigurationGTF&, const int natoms, const arma::mat& motif);
        // Method to build the matrix of the normalization prefactor FAC1(m,l)->FAC1[l][m], used in FAC12fun
        std::vector<std::vector<double>> FAC1fun(const int maxL);
        // Method to build the vector of the normalization prefactor FAC2(shell,l)->FAC2[shell], used in FAC12fun
        std::vector<double> FAC2fun(const ConfigurationGTF&, const bool basis_id);
        // Method to build the vector FAC12[orb] attribute
        void FAC12fun(const ConfigurationGTF&, const int maxL, const bool basis_id);
        // Method to build the vector of vectors FAC3[orb][gaussian] attribute 
        void FAC3fun(const bool basis_id);
        // Method to build the unordered_map attribute g_coefs, which contains the g^{l,m}_{i,j,k} expansion coefficients
        void gfun(const int maxL);
        // Method to build the unordered_map attribute triangInd_to_rowcol, which contains the inverse of the (bijective) function: s(i,j) = j + i(i+1)/2
        void triangIndfun();

    protected:
        int factorial(int);
        int doubleFactorial(int);
        // Method to compute the E^{i,i'}_{t} coefficients in the Hermite Gaussian expansion 
        // G_{i}G_{i'}=\sum_{t=0}^{i+i'}E^{i,i'}_{t}\Lambda_{t}, for i,i'<=4. Returns the vector for 0 <= t <= (i+i'), in increasing t order.
        // The values of the pair (i,i') with i>=i' are in a bijection with the "index" argument as index(i,i')= i' + i(i+1)/2.
        // For i<i', it must be called as Efun(index(i',i),p,PB,PA)  
        arma::colvec Efun(const int index, const double p, const double PA, const double PB);
        
};

}