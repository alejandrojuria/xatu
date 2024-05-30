#pragma once
#include "xatu/IntegralsBase.hpp"

namespace xatu {

/**
 * The IntegralsOverlap2C class is designed to compute and store the two-center overlap integrals in the SCF or AUXILIARY basis set. 
 * It will be called for the AUXILIARY basis via its children class IntegralsOverlap3C basis with the OVERLAP METRIC. 
 * Exclusive to the Gaussian mode
 */
class IntegralsOverlap2C : public virtual IntegralsBase {

    protected:
        IntegralsOverlap2C() = default;
    public:
        IntegralsOverlap2C(const IntegralsBase&, const uint32_t nR, const std::string& intName = "", const bool basis_id = false); 

    private:
        // Method to compute the overlap matrices in the auxiliary (if basis_id == false) or SCF (if basis_id == true) basis 
        // (<P,0|P',R> or <mu,0|mu',R>) for the first nR Bravais vectors R. These first nR (at least, until the star of vectors is 
        // completed) are generated with IntegralsBase::generateRlist. The resulting cube (third dimension spans the Bravais vectors)
        // is saved in the o2Mat_intName.o2c file, and the list of Bravais vectors in atomic units is saved in the RlistAU_intName.o2c 
        // file. The basis_id parameter determines the basis for which the integrals will be computed: true => SCF, false => AUX.
        void overlap2Cfun(const uint32_t nR, const std::string& intName, const bool basis_id = false);

    protected:
        // Analogous to Efun in the parent class IntegralsBase, but only returns the t=0 component. 
        double Efunt0(const int index, const double p, const double PA, const double PB);

};

}