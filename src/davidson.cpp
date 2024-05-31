#include "xatu/davidson.hpp"

namespace xatu {

void davidson_method(
    arma::vec& eigval, 
    arma::cx_mat& eigvec, 
    const arma::cx_mat& mat, 
    const int64_t neigval, 
    const double tol){

    int64_t max_iterations = mat.n_rows/2;
    int64_t k = neigval * 2;
        
    if(!mat.is_hermitian()){
        throw std::invalid_argument("davidson_method: provided matrix must be hermitian");
    }
    
    arma::cx_mat guess_eigvec = arma::eye<arma::cx_mat>(mat.n_rows, k);
    arma::cx_mat proyected_matrix;
    arma::vec aux_eigval = arma::ones(neigval);
    arma::cx_mat Q, R;

    for(int64_t i = 0; i < max_iterations; i++){

        // QR decomposition
        Q.clear();
        arma::qr_econ(Q, R, guess_eigvec);

        // Proyect matrix on
        proyected_matrix = Q.t()*mat*Q;
        arma::eig_sym(eigval, eigvec, proyected_matrix);

        // Add new eigvec to guess
        eigvec = Q.cols(0, (i + 1)*k - 1) * eigvec.cols(0, k - 1);
        for (int64_t j = 0; j < k; j++){
            arma::cx_mat w = (mat - eigval(j)*arma::eye<arma::cx_mat>(mat.n_rows, mat.n_cols))*eigvec.col(j);
            arma::cx_vec normalized_w = w/(eigval(j) - mat(j, j));
            Q.insert_cols(Q.n_cols, normalized_w);
        }

        // Check convergence
        if(arma::norm(eigval.subvec(0, neigval - 1) - aux_eigval) < tol){
            break;
        }

        // Prepare data for next iteration
        aux_eigval = eigval.subvec(0, neigval - 1);
        guess_eigvec.clear();
        guess_eigvec = Q;

        if (i == max_iterations - 1){
            std::cout << "Reached maximum number of iterations" << std::endl;
        }
    }

    // Store final eigenvectors
    eigval = eigval.subvec(0, neigval - 1);
    arma::cout << eigvec.n_rows << arma::endl;
    arma::cout << eigvec.n_cols << arma::endl;
}

}