#pragma once
#include "xatu/IntegralsBase.hpp"

namespace xatu {

/// @brief The Overlap2Centers class is designed to compute the two-center overlap integrals in the
/// auxiliary basis set. It is called if and only if the OVERLAP METRIC is chosen.
class Overlap2Centers : public virtual IntegralsBase {

    //// Methods
    public:
        Overlap2Centers(const GTFConfiguration&, const std::string& overlap2Matrices_filename = "overlap2Matrices", const bool comp = true);
        Overlap2Centers(const IntegralsBase&, const std::string& overlap2Matrices_filename = "overlap2Matrices"); 
        ~Overlap2Centers(){};

    protected:
    // Analogous to Efun in the parent class IntegralsBase, but only returns the t=0 component. 
    double Efunt0(const int index, const double p, const double PA, const double PB);
    // Computes the overlap matrices in the auxiliary basis <P,0|P',R> for the first nR Bravais vectors R, 
    // where nR <= ncells (attribute of IntegralsBase and third argument of GTFConfiguration 's constructor).
    // The resulting cube (third dimension spans the Bravais vectors) is saved in the overlap2Matrices_filename.o2c file.
    void overlap2Cfun(const int nR, const std::string& overlap2Matrices_filename);

};

}