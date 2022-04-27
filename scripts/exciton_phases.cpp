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
    int ncell = 20;
    double filling = 1./2;
    arma::rowvec Q = {0., 0., 0.};
    arma::rowvec parameters = {1., 1., 10};
    std::string modelfile = argv[1];    

    // ----------------- Model parameters & Output --------------------
    bool writeEigvals = false;
    std::string filename = "eigval.out";
    FILE* textfile_en = fopen(filename.c_str(), "w");
        
    bool writeWF = true;
    std::string filename_wf = "phases.out";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");

    // -------------------------- Main body ---------------------------

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "N. cells: " << ncell*ncell << endl;
    cout << "#bands: " << nbands << endl;
    cout << "#removed bands: " << nrmbands << "\n" << endl;

    GExciton bulkExciton = GExciton(modelfile, ncell, nbands, nrmbands, parameters);
    //bulkExciton.setGauge("atomic");

    cout << "Valence bands:\n" << bulkExciton.valenceBands << endl;
    cout << "Conduction bands:\n" << bulkExciton.conductionBands << endl;
    cout << "Gauge used: " << bulkExciton.gauge << "\n" << endl;
    
    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    bulkExciton.brillouinZoneMesh(ncell);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian({}, false);
    auto results = bulkExciton.diagonalize();
    int stateindex = 0;

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    cout << "#--------------------------------- Ground state -----------------------------" << endl;
    cout << "Fundamental state energy is: " << results.eigval(0) << endl;
    cout << "Binding energy: " << results.bindingEnergy(0) << endl;
    cout << "Kinectic energy: "  << results.kineticEnergy(0) << endl;
    cout << "Potential energy: " << results.potentialEnergy(0) << "\n" << endl;
    cout << "#-------------------------------- Excited states ----------------------------" << endl;
    cout << "Second state energy is: " << results.eigval(1) << endl;
    cout << "Third state energy is: " << results.eigval(2) << endl;
    cout << "Fourth state energy is: " << results.eigval(3) << endl;
    cout << "Fifth state energy is: " << results.eigval(4) << endl;
    cout << "Sixth state energy is: " << results.eigval(5) << endl;
    cout << "Seventh state energy is: " << results.eigval(6) << endl;
    cout << "Eigth state energy is: " << results.eigval(7) << endl;
    cout << "Delta E: " << results.eigval(0) - results.eigval(1) << endl;

    if(writeEigvals){
        std::cout << "Writing eigvals to file..." << std::endl;
        fprintf(textfile_en, "%d\n", ncell);
        results.writeEigenvalues(textfile_en, 7);
    }
    
    if(writeWF){
        std::cout << "Writing " << stateindex << "th k wavefunction..." << std::endl;
        arma::cx_mat states = results.symmetrizeStates(results.eigvec.col(0), results.eigvec.col(1));
        results.writeExtendedPhase(states.col(0), textfile_wf);

    };

    fclose(textfile_en);
    fclose(textfile_wf);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};