#include <armadillo>
#include <stdlib.h>

using namespace arma;

int main(){

    mat A = arma::zeros(3,3);
    mat* B;

    B = &A;

    *B = arma::eye(3,3);

    cout << A << endl;
    cout << B << endl;

    return 0;
};