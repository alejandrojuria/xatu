#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "GExciton.hpp"
#include "System.hpp"
#include "libwavefunction.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    if (argc != 2){
		cout << "Error: One argument must be given" << endl;
		exit(1);
	}
    std::string modelfile = argv[1];

    auto start = high_resolution_clock::now();
    
    int nbands = 4;
    int nrmbands = 2;
    int Ncell = 31;
    double filling = 5./8;
    int nQ = 1;
    System ribbon = System(modelfile);
    bool storeAllVectors = true;

    arma::vec QArray = arma::linspace(-PI/ribbon.a + 0.01, PI/ribbon.a - 0.01, nQ);
    arma::mat Qpoints = arma::zeros(nQ, 3);
    Qpoints.col(1) = QArray;
    Qpoints = {{0., 0., 0.}};

    // ----------------- Model parameters & Output --------------------

    std::string filename = "bi_entanglement_spectrum";
    FILE* textfile_en = fopen(filename.c_str(), "w");
    bool writeEntanglement = true;

    // ----------------------- Main body -----------------------

    cout << "##############################################################################" << endl;
    cout << "--------------------------------- Parameters ---------------------------------" << endl;
    cout << "N. cells: " << Ncell << endl;
    cout << "#bands: " << nbands << endl;
    cout << "#removed bands: " << nrmbands << "\n" << endl;
    
    for (int n = 0; n < nQ; n++){
        arma::rowvec Q = Qpoints.row(n);
        cout << "Q: " << Q << endl;
        arma::cx_mat eigvecX;
        vec eigvalX;

        GExciton bulkExciton = GExciton(modelfile, Ncell, Q, nbands, nrmbands, filling, true, storeAllVectors);
        int nk = bulkExciton.nk;
        cout << nk << endl;
        cout << "Basis size: " << bulkExciton.basisStates.n_rows << endl;

        bulkExciton.BShamiltonian();
        arma::eig_sym(eigvalX, eigvecX, bulkExciton.HBS);
        cout << "Ground state energy: " << eigvalX(0) << endl;
        cout << "First excited state energy: " << eigvalX(1) << endl;
        cout << "Second excited state energy: " << eigvalX(2) << endl;
        cout << "Thirds excited state energy: " << eigvalX(3) << endl;
        cx_vec BSEcoefs = eigvecX.col(0);
        cout << "TB basis dim: " << bulkExciton.basisdim << endl;
        cout << "Fermi level: " << bulkExciton.fermiLevel << endl;

        for(int kIndex = 0; kIndex < nk; kIndex++){
            arma::cx_mat density = arma::zeros<arma::cx_mat>(bulkExciton.basisdim/2, bulkExciton.basisdim/2);
            for (int i = 0; i < bulkExciton.basisdim/2; i++){
                for (int j = 0; j < bulkExciton.basisdim/2; j++){
                    density(i, j) = densityMatrixK(kIndex, bulkExciton, BSEcoefs, i, j);
                }
            }
            arma::cx_mat eigvecEntang;
            arma::vec eigvalEntangl;
            arma::eig_sym(eigvalEntangl, eigvecEntang, density);

            if(writeEntanglement){
                    for(int i = 0; i < eigvalEntangl.n_elem; i++){
                        fprintf(textfile_en, "%10.7lf\t", eigvalEntangl(i));
                    }
                    fprintf(textfile_en, "\n");
            };
        }
        
        cout << "##############################################################" << endl;
    }

    fclose(textfile_en);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};