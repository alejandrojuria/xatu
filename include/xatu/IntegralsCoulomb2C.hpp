#pragma once
#include "xatu/IntegralsBase.hpp"
#include "xatu/asa239.hpp"

namespace xatu {

/**
 * The IntegralsCoulomb2C class is designed to compute the two-center Coulomb integrals in the AUXILIARY basis set. 
 * It is called irrespective of the metric. Exclusive to the Gaussian mode
 */
class IntegralsCoulomb2C : public IntegralsBase {

    protected:
        IntegralsCoulomb2C() = default;
    public:
        IntegralsCoulomb2C(const IntegralsBase&, const uint32_t nR, const std::string& intName = ""); 

    private:
        // Method to compute the Coulomb matrices in the auxiliary basis (<P,0|V_c|P',R>) for the first nR Bravais vectors R. 
        // These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
        // The resulting cube (third dimension spans the Bravais vectors) is saved in the C2Mat_intName.C2c file, and the list of 
        // Bravais vectors in atomic units is saved in the RlistAU_intName.C2c file 
        void Coulomb2Cfun(const uint32_t nR, const std::string& intName);

    protected:
        // Analogous to Efun in the parent class IntegralsBase, but restricted to i'=0 and setting PA=0. These are the 
        // expansion coefficients of a single (cartesian) GTF in Hermite Gaussians with the same exponent and center
        arma::colvec Efun_single(const int i, const double p);
        // Method to compute and return the Boys function F_{n}(arg) = \int_{0}^{1}t^{2n}exp(-arg*t^2)dt. It is computed
        // with the lower incomplete Gamma function as: F_{n}(arg) = Gamma(n+0.5)*IncGamma(n+0.5,arg)/(2*arg^(n+0.5)), see (9.8.20)-Helgaker
        double Boysfun(const int n, const double arg);
        // Method to compute and return the auxiliary Hermite Coulomb integral R^{n}_{0,0,0} = (-2p)^n *F_{n}(arg), see (9.9.14)-Helgaker
        double Rn000(const int n, const double p, const double arg);
        // Method to compute and return the Hermite Coulomb integral R^{0}_{t,u,v}(r,(X,Y,Z)), see (9.9.9)-Helgaker
        double HermiteCoulomb(const int t, const int u, const int v, const double p, const double X, const double Y, const double Z);

};

}