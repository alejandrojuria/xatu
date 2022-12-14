#include <armadillo>

namespace xatu {
    void davidson_method(arma::vec&, arma::cx_mat&, const arma::cx_mat&, int neigval = 4, double tol = 1E-8);
}
