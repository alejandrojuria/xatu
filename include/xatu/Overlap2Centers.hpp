#pragma once
#include "xatu/IntegralsBase.hpp"

namespace xatu {

/// @brief The Overlap2Centers class is designed to compute the two-center overlap integrals in the
/// SCF or auxiliary basis set. It will be called for the SCF basis (via its children class Overlap2MixCenters) irrespective of the
/// metric, and it will also be called for the auxiliary basis with the OVERLAP METRIC.
class Overlap2Centers : public virtual IntegralsBase {

    //// Methods
    public:
        Overlap2Centers(const GTFConfiguration&, const std::string& o2Mat_name = "o2Mat", const bool comp = true);
        Overlap2Centers(const IntegralsBase&, const std::string& o2Mat_name = "o2Mat_AUX", const bool basis_id = false); 
        ~Overlap2Centers(){};

    protected:
    // Analogous to Efun in the parent class IntegralsBase, but only returns the t=0 component. 
    double Efunt0(const int index, const double p, const double PA, const double PB);
    // Method to compute the overlap matrices in the auxiliary (if basis_id == false) or SCF (if basis_id == true) basis (<P,0|P',R> or <mu,0|mu',R>) 
    // for the first nR Bravais vectors R, where nR <= ncells (attribute of IntegralsBase and third argument of GTFConfiguration 's constructor).
    // The resulting cube (third dimension spans the Bravais vectors) is saved in the o2Mat_name.o2c file. 
    // The basis_id parameter determines the basis for which the integrals will be computed: true => SCF, false => AUX.
    void overlap2Cfun(const int nR, const std::string& o2Mat_name, const bool basis_id = false);

};

}