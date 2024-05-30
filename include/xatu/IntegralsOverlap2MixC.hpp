#pragma once
#include "xatu/IntegralsOverlap2C.hpp"

namespace xatu {

/**
 * The IntegralsOverlap2MixC class is designed to compute the two-center overlap integrals in the
 * mixed SCF and auxiliary basis set: <P,0|mu,R>. It will be called irrespective of the metric.
 * Exclusive to the Gaussian mode, but in principle not needed (included due to possible future use for truncations)
 */
 class IntegralsOverlap2MixC final : public IntegralsOverlap2C {

    public:
        IntegralsOverlap2MixC(const IntegralsBase&, const uint32_t nR, const std::string& intName = ""); 

    private:
        // Method to compute the overlap matrices in the mixed SCF and auxiliary basis sets (<P,0|mu,R>) for the first nR Bravais 
        // vectors R. These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
        // The resulting cube (third dimension spans the Bravais vectors) is saved in the o2MixMat_intName.o2mc file, and the list of 
        // Bravais vectors in atomic units is saved in the RlistAU_intName.o2mc file
        void overlap2MixCfun(const uint32_t nR, const std::string& intName);

};

}