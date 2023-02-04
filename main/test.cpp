#include <complex>
#include <armadillo>
#include <GExciton.hpp>

extern "C" {
    void skubo_w_(int* nR, int* norb, int* norb_ex, int* nv, int* nRvec, double* bravaisLattice, double* motif, 
                  double* hhop, double* shop, int* nk, double* rkx, double* rky, double* rkz, 
                  std::complex<double>* fk_ex, double* e_ex);
}

int main(int argc, char* argv[]){
    
    int nbands = 1;
    int nrmbands = 0;
    int ncell = 50;
    arma::rowvec Q = {0., 0., 0.};
    double eps_s = 3.9;
    double eps   = 40;
    double r0    = 1.589*eps/(1 + eps_s);
    
    arma::rowvec parameters = {1., 1., 10.};
    std::string modelfile = argv[1];    

    // -------------------------- Main body ---------------------------

    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                  Parameters                               |" << std::endl;
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "N. cells: " << ncell*ncell << std::endl;
    std::cout << "#bands: " << nbands << std::endl;
    std::cout << "#removed bands: " << nrmbands << std::endl;
    std::cout << "System configuration file: " << modelfile << "\n" << std::endl;

    SystemConfiguration config = SystemConfiguration(modelfile);
    GExciton bulkExciton = GExciton(config, ncell, nbands, nrmbands, parameters);
    
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                Initialization                             |" << std::endl;
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;

    bulkExciton.brillouinZoneMesh(ncell);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize();
    printEnergies(results);

    int nR = bulkExciton.unitCellList.n_rows;
    int norb = bulkExciton.basisdim;
    int norb_ex = bulkExciton.excitonbasisdim;
    int nv = bulkExciton.valenceBands.n_elem;
    int nc = bulkExciton.conductionBands.n_elem;
    //arma::Mat<int> nRvec = {{0, 0, 0}, {-1, 0, 0}, {-1, 1, 0}, {0, -1, 0}, 
    //                        {0, 1, 0}, {1, -1, 0}, {1, 0, 0}};
    arma::Mat<int> nRvec = {{0,0,0}, {-1,0,0}, {0, -1, 0}, {0, 1, 0}, {1,0,0}};
    // Extend bravais lattice to 3x3 matrix
    arma::mat R = arma::zeros(3, 3);
    for (int i = 0; i < bulkExciton.bravaisLattice.n_rows; i++){
        R.row(i) = bulkExciton.bravaisLattice.row(i);
    }

    arma::mat B = bulkExciton.motif.cols(0, 2);
    arma::cube hhop = arma::real(bulkExciton.hamiltonianMatrices);
    arma::cube shop(arma::size(hhop));
    if (bulkExciton.overlapMatrices.empty()){
        for (int i = 0; i < hhop.n_slices; i++){
            shop.slice(i) = arma::eye(size(hhop.slice(i)));
        }
    }
    else{
        shop = arma::real(bulkExciton.overlapMatrices);
    }
    int nk = bulkExciton.nk;
    arma::vec rkx = bulkExciton.kpoints.col(0);
    arma::vec rky = bulkExciton.kpoints.col(1);
    arma::vec rkz = bulkExciton.kpoints.col(2);

    arma::cx_mat eigvec = results.eigvec;
    arma::vec eigval = results.eigval;


    skubo_w_(&nR, &norb, &norb_ex, &nv, nRvec.memptr(), R.memptr(), B.memptr(), hhop.memptr(), shop.memptr(), &nk, rkx.memptr(),
             rky.memptr(), rkz.memptr(), eigvec.memptr(), eigval.memptr());
}

