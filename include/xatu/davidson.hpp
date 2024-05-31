#include <armadillo>

namespace xatu {
    void davidson_method(arma::vec&, arma::cx_mat&, const arma::cx_mat&, const int64_t neigval = 4, const double tol = 1E-8);
}
