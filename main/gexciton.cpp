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

    bool writeStates = true;
    std::string filename_st = "states.out";
    FILE* textfile_st = fopen(filename_st.c_str(), "w");
        
    bool writeWF = false;
    std::string filename_kwf = "kwf.out";
    FILE* textfile_kwf = fopen(filename_kwf.c_str(), "w");

    bool writeRSWF = false;
    std::string filename_rswf = "rswf.out";
    FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");

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

    cout << "-----------------------------------------------------------------------------" << endl;
    cout << "|                                    Output                                 |" << endl;
    cout << "-----------------------------------------------------------------------------" << endl;

    if(writeEigvals){
        std::cout << "Writing eigvals to file: " << filename << std::endl;
        fprintf(textfile_en, "%d\n", ncell);
        results.writeEigenvalues(textfile_en, 7);
    }

    if(writeStates){
        std::cout << "Writing states to file: " << filename_st << std::endl;
        results.writeStates(textfile_st);
    }
    
    if(writeWF){
        std::cout << "Writing k w.f. to file: " << filename_kwf << std::endl;
        int nstates = 8;
        for(int stateindex = 0; stateindex < nstates; stateindex++){
            results.writeExtendedReciprocalAmplitude(stateindex, textfile_kwf);        
        }
    }

    if(writeRSWF){
        std::cout << "Writing real space w.f. to file: " << filename_rswf << std::endl;
        int holeIndex = 0;
        arma::rowvec holeCell = {0., 0., 0.};
        results.writeRealspaceAmplitude(0, holeIndex, holeCell, textfile_rswf);
    }

    fclose(textfile_en);
    fclose(textfile_st);
    fclose(textfile_kwf);
    fclose(textfile_rswf);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};