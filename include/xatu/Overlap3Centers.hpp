#pragma once
#include "xatu/Overlap2Centers.hpp"

namespace xatu {

/// @brief The Overlap3Centers class is designed to compute the three-center overlap integrals in the mixed
/// SCF and auxiliary basis sets. It is called if and only if the OVERLAP METRIC is chosen.
class Overlap3Centers : public virtual Overlap2Centers {

    //// Methods
    public:
        Overlap3Centers(const GTFConfiguration&, const int tol, const uint32_t nR2, const std::string& o3Mat_name = "o3Mat", const std::string& o2Mat_name = "o2Mat");
        Overlap3Centers(const IntegralsBase&, const int tol, const uint32_t nR2, const std::string& o3Mat_name = "o3Mat", const std::string& o2Mat_name = "o2Mat_AUX"); 
        ~Overlap3Centers(){};

    protected:
    // Analogous to Efunt0 in the parent class Overlap2Centers, but instead including E^{i,i'}_{0} for (5<=i<=8,i'<=4). 
    // The values of the pair (i,i') with i>=i' are in a bijection with the "index" argument as index(i,i')= 5*(i-5) + i'.
    // For i<i', it must be called as Efun(index(i',i),p,PB,PA).
    double EfunTriplet0(const int index, const double p, const double PA, const double PB);
    // Coefficients D^{t}_{l}(p) in the expansion \Lambda_{t}(x,p,Px) = D^{t}_{l}(p) *(x-Px)^l *exp(-p(x-Px)^2).
    // The argument t spans 0 <= t <= 8, and the entries of the returned vector are the nonzero D^{t}_{l} for the given t.  
    // There are ceil((t+1)/2) of such entries, each corresponding to l = t, t-2, t-4, ... 1 (or 0)
    std::vector<double> Dfun(const int t, const double p);
    // Method to compute the rectangular overlap matrices <P,0|mu,R;mu',R'> in the mixed SCF and auxiliary basis sets for the first nR2 Bravais
    // vectors R and R' (nR2^2 pairs of vectors). These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
    // Each entry above a certain tolerance (10^-tol) is stored in an entry of a vector (of arrays) along with the corresponding indices:
    // P,mu,mu',R,R',value; in that order. The vector is saved in the o3Mat_name.o3c file.
    void overlap3Cfun(const uint32_t nR2, const int tol, const std::string& o3Mat_name);

};

}