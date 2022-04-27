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
#include "Crystal.hpp"
#include "Result.hpp"
#include "utils.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    if (argc != 2){
		throw std::invalid_argument("Error: One input file is expected");
	}

    auto start = high_resolution_clock::now();
    
    int nbands = 1;
    int nrmbands = 0;
    int ncell = 60;
    double filling = 1./2;
    arma::rowvec Q = {0., 0., 0.};
    arma::rowvec parameters = {1., 1., 10};
    std::string modelfile = argv[1];
    
    bool useCenteredMeshK = false;

    // ----------------- Model parameters & Output --------------------

    std::string filename = "eigval_reduced.out";
    FILE* textfile_en = fopen(filename.c_str(), "w");
    bool writeEigvals = false;

    std::string filename_wf = "k_wf_spectrum_antonio.out";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
    bool writeWF = true;

    bool calculateSpin = false;
    std::string filename_spin = "spin.out";
    FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
    bool writeSpin = false;

    std::string filename_cells = "cells.out";
    FILE* textfile_cells = fopen(filename_cells.c_str(), "w");

    // ----------------------- Main body ----------------ExcitonState-------

    fprintf(textfile_en, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");

    cout << "##############################################################################" << endl;
    cout << "--------------------------------- Parameters ---------------------------------" << endl;
    cout << "N. cells: " << ncell*ncell << endl;
    cout << "#bands: " << nbands << endl;
    cout << "#removed bands: " << nrmbands << "\n" << endl;

    GExciton bulkExciton = GExciton(modelfile, ncell, nbands, nrmbands, parameters);

    cout << "Valence bands:" << bulkExciton.valenceBands << endl;
    cout << "Conduction bands:" << bulkExciton.conductionBands << endl;
    //arma::rowvec K = (bulkExciton.reciprocalLattice.row(0) 
    //                - bulkExciton.reciprocalLattice.row(1))/2;

    double norm = arma::norm(bulkExciton.reciprocalLattice.row(0));
    arma::rowvec K = norm/sqrt(3)*(bulkExciton.reciprocalLattice.row(0)/2 - bulkExciton.reciprocalLattice.row(1)/2)/(
                    arma::norm(bulkExciton.reciprocalLattice.row(0)/2 - bulkExciton.reciprocalLattice.row(1)/2));
    
    bulkExciton.brillouinZoneMesh(ncell);
    //bulkExciton.reducedBrillouinZoneMesh(ncell, ncell*3);
    //bulkExciton.shiftBZ(K);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize();
    int stateindex = 6;

    cout << "#----------------------- Energy -----------------------" << endl;
    cout << "Fundamental state energy is: " << results.eigval(0) << endl;
    cout << "Binding energy: " << results.bindingEnergy(0) << endl;
    cout << "Kinectic energy: "  << results.kineticEnergy(0) << endl;
    cout << "Potential energy: " << results.potentialEnergy(0) << "\n" << endl;
    cout << "#----------------------- Excited states -----------------------" << endl;
    cout << "Second state energy is: " << results.eigval(1) << endl;
    cout << "Second state binding energy is: " << results.eigval(1) - 7.25 << endl;

    cout << "Third state energy is: " << results.eigval(2) << endl;
    cout << "Third state binding energy is: " << results.eigval(2) - 7.25 << endl;

    cout << "Fourth state energy is: " << results.eigval(3) << endl;
    cout << "Fourth state binding energy is: " << results.eigval(3) - 7.25 << endl;

    cout << "Fifth state energy is: " << results.eigval(4) << endl;
    cout << "Fifth state binding energy is: " << results.eigval(4) - 7.25 << endl;

    cout << "Sixth state energy is: " << results.eigval(5) << endl;
    cout << "Sixth state binding energy is: " << results.eigval(5) - 7.25 << endl;

    cout << "Seventh state energy is: " << results.eigval(6) << endl;
    cout << "Seventh state binding energy is: " << results.eigval(6) - 7.25 << endl;

    cout << "Delta E: " << results.eigval(0) - results.eigval(1) << endl;
    
    if(writeWF){
        writeVectorsToFile(K, textfile_wf);
        std::cout << "Writing " << stateindex << "th k wavefunction..." << std::endl;
        results.writeExtendedReciprocalAmplitude(0, textfile_wf);
        results.writeExtendedReciprocalAmplitude(1, textfile_wf);
        arma::cx_mat combinations = results.diagonalizeC3(arma::vec{2, 3});
        arma::cx_vec coefs = combinations.col(0)(0) * results.eigvec.col(2) + 
                             combinations.col(0)(1) * results.eigvec.col(3);
        results.writeExtendedReciprocalAmplitude(coefs, textfile_wf);
        coefs = combinations.col(1)(0) * results.eigvec.col(2) + 
                combinations.col(1)(1) * results.eigvec.col(3);
        results.writeExtendedReciprocalAmplitude(coefs, textfile_wf);
        combinations = results.diagonalizeC3(arma::vec{4, 5});
        coefs = combinations.col(0)(0) * results.eigvec.col(4) + 
                combinations.col(0)(1) * results.eigvec.col(5);
        results.writeExtendedReciprocalAmplitude(coefs, textfile_wf);
        coefs = combinations.col(1)(0) * results.eigvec.col(4) + 
                combinations.col(1)(1) * results.eigvec.col(5);
        results.writeExtendedReciprocalAmplitude(coefs, textfile_wf);
        results.writeExtendedReciprocalAmplitude(6, textfile_wf);
        results.writeExtendedReciprocalAmplitude(7, textfile_wf);
        arma::mat C3 = bulkExciton.C3ExcitonBasisRep();
        arma::cx_mat commutator = bulkExciton.HBS*C3 - C3*bulkExciton.HBS;
        arma::cout << "Norm [HBS, C3]: " << arma::norm(commutator) << arma::endl;
        arma::cout << "Norm [HK, C3]: " << arma::norm(bulkExciton.HK*C3 - C3*bulkExciton.HK) << arma::endl;
    };

    fclose(textfile_en);
    fclose(textfile_wf);
    fclose(textfile_spin);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};