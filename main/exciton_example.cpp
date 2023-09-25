#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include <xatu.hpp>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    // Parse console stdin

    if (argc != 2){
		throw std::invalid_argument("Error: One input file is expected");
	}
    else if (argc < 2){
        throw std::invalid_argument("Error: At least one input file is required (system config.)");
    };

    auto start = high_resolution_clock::now();
    
    int nbands = 1;
    int nrmbands = 0;
    int ncell = 40;
    arma::rowvec Q = {0., 0., 0.};
    
    arma::rowvec parameters = {1., 1., 10.};
    std::string modelfile = argv[1];    

    int nstates = 8;

    // ----------------- Model parameters & Output --------------------
    bool writeEigvals = false;
    std::string filename = "eigval.out";
    FILE* textfile_en = fopen(filename.c_str(), "w");

    bool writeStates = false;
    std::string filename_st = "states.out";
    FILE* textfile_st = fopen(filename_st.c_str(), "w");
        
    bool writeWF = true;
    std::string filename_kwf = "kwf.out";
    FILE* textfile_kwf = fopen(filename_kwf.c_str(), "w");

    bool writeRSWF = true;
    std::string filename_rswf = "rswf.out";
    FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");

    // -------------------------- Main body ---------------------------

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "N. cells: " << ncell*ncell << endl;
    cout << "#bands: " << nbands << endl;
    cout << "#removed bands: " << nrmbands << endl;
    cout << "System configuration file: " << modelfile << "\n" << endl;

    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);
    xatu::Exciton bulkExciton = xatu::Exciton(config, ncell, nbands, nrmbands, parameters);
    arma::cout << "Orbitals: " << bulkExciton.orbitals << arma::endl;
    bulkExciton.setMode("realspace");

    cout << "Valence bands:\n" << bulkExciton.valenceBands << endl;
    cout << "Conduction bands:\n" << bulkExciton.conductionBands << endl;
    cout << "Gauge used: " << bulkExciton.gauge << "\n" << endl;
    
    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    bulkExciton.brillouinZoneMesh(ncell);
    bulkExciton.initializeHamiltonian();
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize("diag", nstates);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    printEnergies(results, nstates, 6);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Output                                 |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    if(writeEigvals){
        std::cout << "Writing eigvals to file: " << filename << std::endl;
        fprintf(textfile_en, "%d\n", ncell);
        results.writeEigenvalues(textfile_en, 8);
    }

    if(writeStates){
        std::cout << "Writing states to file: " << filename_st << std::endl;
        results.writeStates(textfile_st, nstates);
    }
    
    if(writeWF){
        std::cout << "Writing k w.f. to file: " << filename_kwf << std::endl;
        int nstates = 8;
        for(int stateindex = 0; stateindex < nstates; stateindex++){
            results.writeExtendedReciprocalAmplitude(stateindex, textfile_kwf);   
        }
    }

    if(writeRSWF){
        arma::uvec statesToWrite = arma::regspace<arma::uvec>(0, 7);
        std::cout << "Writing real space w.f. to file: " << filename_rswf << std::endl;
        arma::rowvec holeCell = {0., 0., 0.};
        int holeIndex = 0;
        
        for(unsigned int i = 0; i < statesToWrite.n_elem; i++){
            std::cout << "Writing state " << i + 1 << " out of " << statesToWrite.n_elem << std::endl;
            results.writeRealspaceAmplitude(statesToWrite(i), holeIndex, holeCell, textfile_rswf, 10);
        }
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