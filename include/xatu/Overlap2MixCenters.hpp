#pragma once
#include "xatu/Overlap2Centers.hpp"

namespace xatu {

/// @brief The Overlap2Centers class is designed to compute the two-center overlap integrals in the
/// mixed SCF and auxiliary basis set: <P,0|mu,R>. It will be called irrespective of the metric.
class Overlap2MixCenters : public virtual Overlap2Centers {

    //// Methods
    public:
        Overlap2MixCenters(const IntegralsBase&, const uint32_t nR, const std::string& o2MixMat_name = "o2MixMat", const std::string& o2Mat_name = "o2Mat_SCF"); 
        ~Overlap2MixCenters(){};

    protected:
    // Method to compute the overlap matrices in the mixed SCF and auxiliary basis sets (<P,0|mu,R>) for the first nR Bravais vectors R. 
    // These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
    // The resulting cube (third dimension spans the Bravais vectors) is saved in the o2MixMat_name.o2mc file.
    void overlap2MixCfun(const uint32_t nR, const std::string& o2MixMat_name);

};

}