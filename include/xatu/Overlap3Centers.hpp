#pragma once
#include "xatu/Overlap2Centers.hpp"

namespace xatu {

/// @brief The Overlap3Centers class is designed to compute the three-center overlap integrals in the mixed
/// SCF and auxiliary basis sets. It is called if and only if the OVERLAP METRIC is chosen.
class Overlap3Centers : public virtual Overlap2Centers {

    //// Methods
    public:
        Overlap3Centers(const GTFConfiguration&, const std::string& o3Mat_name = "o3Mat", const std::string& o2Mat_name = "o2Mat", const std::string& o2Mat_read_SCF = "o2Mat_SCF", const std::string& o2Mat_read_Mix = "o2MixMat");
        Overlap3Centers(const IntegralsBase&, const std::string& o3Mat_name = "o3Mat", const std::string& o2Mat_name = "o2Mat_AUX", const std::string& o2Mat_read_SCF = "o2Mat_SCF", const std::string& o2Mat_read_Mix = "o2MixMat"); 
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
    // Method to compute the rectangular overlap matrices <P,0|mu,R;mu',R'> in the mixed SCF and auxiliary basis sets for the first nR Bravais
    // vectors R and R' (nR^2 in total), where nR <= ncells (attribute of IntegralsBase and third argument of GTFConfiguration 's
    // constructor). The resulting cube (third dimension spans auxiliary basis {P}, while the Bravais lattice vectors R and R' vary by rows and columns
    // (respectively) once the block for {mu,mu'} has been completed) is saved in the o3Mat_name.o2c file.
    void overlap3Cfun(const int nR, const std::string& o3Mat_name, const std::string& o2Mat_read_SCF, const std::string& o2Mat_read_Mix);

};

}