#include <armadillo>
#include <xatu.hpp>

using namespace xatu;

int main(){

    int nbands = 1;
    int nrmbands = 0;
    int ncell = 40;
    arma::rowvec parameters = {1., 1., 10.};
    int nstates = 8;

    auto config = SystemConfiguration("./examples/material_models/hBN.model");
    auto exciton = Exciton(config, ncell, nbands, nrmbands, parameters);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    return 0;
};