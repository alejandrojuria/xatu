#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>
#include <iomanip>

#include <tclap/CmdLine.h>

#include "GExciton.hpp"
#include "System.hpp"
#include "Crystal.hpp"
#include "Result.hpp"
#include "utils.hpp"
#include "SystemConfiguration.hpp"
#include "CrystalDFTConfiguration.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    auto start = high_resolution_clock::now();

    // Parse CLI arguments
    TCLAP::CmdLine cmd("Command line interface options of the Xatu binary. For a more detailed description, refer to the user guide or the API documentation.", ' ', "1.0");

    TCLAP::ValueArg<int>    statesArg("n", "states", "Specify number of exciton states to show.", false, 8, "No. states", cmd);
    TCLAP::ValueArg<int>    precisionArg("p", "precision", "Desired energy precision. Used to compute degeneracies.", false, 6, "No. decimals", cmd);
    TCLAP::SwitchArg        spinArg("s", "spin", "Compute exciton spin and write it to file.", cmd, false);
    TCLAP::ValueArg<int>    dftArg("d", "dft", "Indicates that the system file is a .outp CRYSTAL file.", false, -1, "No. Fock matrices", cmd);
    TCLAP::ValueArg<double> absorptionArg("a", "absorption", "Computes the absorption spectrum.", false, 0.01, "Broadening", cmd);

    TCLAP::AnyOf         outputOptions;
    TCLAP::SwitchArg     energyArg("e", "energy", "Write energies.", false);
    TCLAP::SwitchArg     eigenstatesArg("c", "eigenstates", "Write eigenstates.", false);
    TCLAP::SwitchArg     reciprocalArg("k", "kwf", "Write reciprocal wavefunction.", false);
    TCLAP::MultiArg<int> realspaceArg("r", "rswf", "Write real-space wavefunction.", false, "Atom index, [no. unit cells]");
    outputOptions.add(energyArg).add(eigenstatesArg).add(reciprocalArg).add(realspaceArg);
    cmd.add(outputOptions);
    TCLAP::SwitchArg outputArg("o", "output", "Write to file information about the excitons.", cmd, false);

    std::vector<std::string> methods = {"diag", "davidson"};
    TCLAP::ValuesConstraint<std::string> allowedMethods(methods);
    TCLAP::ValueArg<std::string> methodArg("m", "method", "Method to solve the Bethe-Salpeter equation.", false, "diag", &allowedMethods, cmd);
    TCLAP::ValueArg<std::string> bandsArg("b", "bands", "Computes the bands of the system on the specified kpoints.", false, "kpoints.txt", "Filename");
    
    TCLAP::UnlabeledValueArg<std::string> systemArg("systemfile", "System file", true, "system.txt", "filename", cmd);
    TCLAP::UnlabeledValueArg<std::string> excitonArg("excitonfile", "Exciton file", false, "exciton.txt", "filename", cmd);

    cmd.parse(argc, argv);

    // Extract information from parsed CLI options
    int nstates        = statesArg.getValue();
    int ncells         = dftArg.getValue();
    bool triangular    = dftArg.isSet();
    int decimals       = precisionArg.getValue();
    float broadening   = absorptionArg.getValue();
    std::string method = methodArg.getValue();
    std::vector<int> rsInfo = realspaceArg.getValue();
    int holeIndex = 0, ncellsRSWF = 8;
    if (rsInfo.size() == 1){
        holeIndex = rsInfo[0];
    }
    else if(rsInfo.size() == 2){
        holeIndex  = rsInfo[0];
        ncellsRSWF = rsInfo[1];
    }
    else if(rsInfo.size() > 2){
        throw std::invalid_argument("-r takes at most two values, holeIndex and ncells");
    }

    std::string systemfile  = systemArg.getValue();
    std::string excitonfile = excitonArg.getValue();

    std::string kpointsfile = bandsArg.getValue();

    // Init. configurations
    std::unique_ptr<SystemConfiguration> systemConfig;
    std::unique_ptr<ExcitonConfiguration> excitonConfig;

    if (dftArg.isSet()){
        systemConfig.reset(new CrystalDFTConfiguration(systemfile, ncells));
    }
    else{
        systemConfig.reset(new SystemConfiguration(systemfile));
    }

    // If bands flag is present, compute bands and exit.
    // Otherwise, init. exciton configuration.
    if (bandsArg.isSet()){
        System system = System(*systemConfig);
        system.solveBands(kpointsfile, triangular);

        return 0;
    }
    else{
        if (!excitonArg.isSet()){
            throw std::invalid_argument("Must provide exciton file.");
        }
        excitonConfig.reset(new ExcitonConfiguration(excitonfile));
    }

    std::string output = systemConfig->systemInfo.name;

    // --------------------------- Output ---------------------------
    bool writeEigvals = energyArg.isSet();
    std::string filename_en = output + ".eigval";
    FILE* textfile_en = fopen(filename_en.c_str(), "w");

    bool writeStates = eigenstatesArg.isSet();
    std::string filename_st = output + ".states";
    FILE* textfile_st = fopen(filename_st.c_str(), "w");
        
    bool writeWF = reciprocalArg.isSet();
    std::string filename_kwf = output + ".kwf";
    FILE* textfile_kwf = fopen(filename_kwf.c_str(), "w");

    bool writeRSWF = realspaceArg.isSet();
    std::string filename_rswf = output + ".rswf";
    FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");

    bool writeAbs = absorptionArg.isSet();
    std::string filename_abs = output + ".abs";
    FILE* textfile_abs = fopen(filename_abs.c_str(), "w");

    // -------------------------- Main body ---------------------------

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;
    
    GExciton bulkExciton = GExciton(*systemConfig, *excitonConfig);
    bulkExciton.setMode(excitonConfig->excitonInfo.mode);

    cout << "\n";
    cout << std::left << std::setw(30) << "System configuration file: " << std::setw(10) << systemfile << endl;
    cout << std::left << std::setw(30) << "Exciton configuration file: " << std::setw(10) << excitonfile << "\n" << endl;
    bulkExciton.printInformation();
    
    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    bulkExciton.brillouinZoneMesh(excitonConfig->excitonInfo.ncell);
    bulkExciton.initializeHamiltonian(triangular);
    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize(method, nstates);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    printEnergies(results, nstates, decimals);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Output                                 |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    if(writeEigvals){
        std::cout << "Writing eigvals to file: " << filename_en << std::endl;
        fprintf(textfile_en, "%d\n", excitonConfig->excitonInfo.ncell);
        results.writeEigenvalues(textfile_en, nstates);
    }

    if(writeStates){
        std::cout << "Writing states to file: " << filename_st << std::endl;
        results.writeStates(textfile_st, nstates);
    }
    
    if(writeWF){
        std::cout << "Writing k w.f. to file: " << filename_kwf << std::endl;
        for(int stateindex = 0; stateindex < nstates; stateindex++){
            results.writeExtendedReciprocalAmplitude(stateindex, textfile_kwf);        
        }
    }

    if(writeRSWF){
        arma::uvec statesToWrite = arma::regspace<arma::uvec>(0, nstates - 1);
        std::cout << "Writing real space w.f. to file: " << filename_rswf << std::endl;
        arma::rowvec holeCell = {0., 0., 0.};
        
        for(unsigned int i = 0; i < statesToWrite.n_elem; i++){
            std::cout << "Writing state " << i + 1 << " out of " << statesToWrite.n_elem << std::endl;
            results.writeRealspaceAmplitude(statesToWrite(i), holeIndex, holeCell, textfile_rswf, ncellsRSWF);
        }
    }

    if(writeAbs){
        std::cout << "Writing absorption spectrum fo file: " << filename_abs << std::endl;
        results.writeAbsorptionSpectrum(textfile_abs, broadening);
    }

    fclose(textfile_en);
    fclose(textfile_st);
    fclose(textfile_kwf);
    fclose(textfile_rswf);
    fclose(textfile_abs);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};