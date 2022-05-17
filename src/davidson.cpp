#include "davidson.hpp"

void davidson_method(
    arma::vec& eigval, 
    arma::cx_vec& eigvec, 
    const arma::cx_mat& mat,
    int k,
    int neigval,
    double tol
){

    if(!mat.is_hermitian()){
        throw std::invalid_argument("davidson_method: provided matrix must be hermitian");
    }
    
    arma::cx_mat guess_eigvec = arma::eye<cx_mat>(mat.n_rows, k);
    bool notConverged = true;
    arma::cx_mat proyected_matrix, aux_eigvec;
    arma::vec aux_eigval;

    while (notConverged){
        

        // Proyect matrix on 
        proyected_matrix = arma::dot(guess_eigvec.t(), arma::dot(mat, guess_eigvec));
        arma::eig_sym(eigval, aux_eigvec, proyected_matrix);

        if(arma::norm(eigval - aux_eigval) < tol){
            notConverged = false;
        }

        eigval = aux_eigval;


        
        

    }
    
    


}
