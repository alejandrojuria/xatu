// Header file with Davidson method to extract the first eigenvalues

#include <armadillo>

void davidson_method(
    arma::vec&, 
    arma::cx_mat&, 
    const arma::cx_mat&, 
    int k = 4, 
    int neigval = 4,
    double tol = 1E-7);