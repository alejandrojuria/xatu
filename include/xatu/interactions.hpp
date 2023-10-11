#include <armadillo>
#include <complex>

#include "xatu.hpp"


namespace xatu {

    // ------ Potentials ------

    /* Keldysh potential */
    double keldysh(double, double, double, double, double, double);
    void STVH0(double, double*);  // Auxiliary routine for Struve function

    /* Coulomb potential */
    double coulomb(double);

    /* Keldysh potential FT */
    double keldyshFT(const arma::rowvec&, double, double, double, double, int);

    /* Coulomb potential FT */
    double coulombFT(double);

    // ------ Interaction matrix elements ------

    /* Real-space interaction */
    std::complex<double> motifFourierTransform(const arma::rowvec&, const arma::rowvec&, const arma::rowvec&, const arma::mat&);
    arma::cx_mat motifFTMatrix(const arma::rowvec&, const arma::mat&, const arma::mat&);
    arma::cx_mat extendMotifFT(const arma::cx_mat&, int, const arma::mat&, const arma::urowvec&);
    arma::cx_vec sumStateOverOrbitals(const arma::rowvec&, const arma::urowvec&);

    std::complex<double> realSpaceInteractionElement(const arma::cx_vec&,
                                                     const arma::cx_vec&,
                                                     const arma::cx_vec&,
                                                     const arma::cx_vec&,
                                                     const arma::cx_mat&,
                                                     const arma::mat&,
                                                     const arma::urowvec&);

    /* Reciprocal space interaction */
    std::complex<double> blochCoherenceFactor(const arma::cx_vec&, const arma::cx_vec&,
                                              const arma::rowvec&, const arma::rowvec&,
                                              const arma::rowvec&, const arma::mat&, const arma::urowvec&);

    std::complex<double> reciprocalSpaceInteractionElement(const arma::cx_vec&,
                                                           const arma::cx_vec&,
                                                           const arma::cx_vec&,
                                                           const arma::cx_vec&,
                                                           const arma::rowvec&,
                                                           const arma::rowvec&,
                                                           const arma::rowvec&,
                                                           const arma::rowvec&,
                                                           int,
                                                           const arma::mat&);

    
}