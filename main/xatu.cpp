#include <math.h>
#include <chrono>
#include <iomanip>
#include <tclap/CmdLine.h>
#include "xatu.hpp"


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
    TCLAP::ValueArg<int>    w90Arg("w", "w90", "Indicates that the system file is a tb.dat Wannier90 file.", false, -1, "No. electrons  ", cmd);
    TCLAP::SwitchArg        absorptionArg("a", "absorption", "Computes the absorption spectrum.", cmd, false);
    TCLAP::ValueArg<std::string> formatArg("f", "format", "Format of the input system file.", false, "model", "model or hdf5", cmd);

    TCLAP::AnyOf         outputOptions;
    TCLAP::SwitchArg     energyArg("e", "energy", "Write energies.", false);
    TCLAP::SwitchArg     eigenstatesArg("c", "eigenstates", "Write eigenstates.", false);
    TCLAP::SwitchArg     reciprocalArg("k", "kwf", "Write reciprocal wavefunction.", false);
    TCLAP::MultiArg<int> realspaceArg("r", "rswf", "Write real-space wavefunction.", false, "Atom index, [no. unit cells]");
    outputOptions.add(energyArg).add(eigenstatesArg).add(reciprocalArg).add(realspaceArg);
    cmd.add(outputOptions);
    TCLAP::SwitchArg outputArg("o", "output", "Write to file information about the excitons.", cmd, false);

    std::vector<std::string> methods = {"diag", "davidson", "sparse"};
    TCLAP::ValuesConstraint<std::string> allowedMethods(methods);
    TCLAP::ValueArg<std::string> methodArg("m", "method", "Method to solve the Bethe-Salpeter equation.", false, "diag", &allowedMethods, cmd);
    TCLAP::ValueArg<std::string> bandsArg("b", "bands", "Computes the bands of the system on the specified kpoints.", false, "kpoints.txt", "Filename", cmd);
    TCLAP::ValueArg<std::string> screeningArg("z", "screening", "Provides input parameters for computing the screening.", false, "model.screening", "Filename", cmd);
    TCLAP::UnlabeledValueArg<std::string> systemArg("systemfile", "System file", true, "system.txt", "filename", cmd);
    TCLAP::UnlabeledValueArg<std::string> excitonArg("excitonfile", "Exciton file", false, "exciton.txt", "filename", cmd);
    cmd.parse(argc, argv);

    // Extract information from parsed CLI options
    int nstates        = statesArg.getValue();
    int ncells         = dftArg.getValue();
    int electronNum    = w90Arg.getValue();
    int decimals       = precisionArg.getValue();
    std::string method = methodArg.getValue();
    std::vector<int> rsInfo = realspaceArg.getValue();
    std::string format = formatArg.getValue();
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
    std::string screeningfile = screeningArg.getValue();
    // Init. configurations
    std::unique_ptr<xatu::SystemConfiguration> systemConfig;
    std::unique_ptr<xatu::ExcitonConfiguration> excitonConfig;
    std::unique_ptr<xatu::ScreeningConfiguration> screeningConfig;

    if (dftArg.isSet()){
        systemConfig.reset(new xatu::CRYSTALConfiguration(systemfile, ncells));
    } else if (w90Arg.isSet()){
        cout << "Parsing w90 file..." << std::endl;
        systemConfig.reset(new xatu::Wannier90Configuration(systemfile, electronNum));
        cout << "Parsed Wannier90 _tb.out" << std::endl;
    }
    else{
        if (format == "hdf5"){
            systemConfig.reset(new xatu::HDF5Configuration(systemfile));
        }
        else if (format == "model"){
            systemConfig.reset(new xatu::SystemConfiguration(systemfile));
        }
        else{
            throw std::invalid_argument("Format not recognized. Use 'model' or 'hdf5'.");
        }
    }

    // If bands flag is present, compute bands and exit.
    // Otherwise, init. exciton configuration.
    if (bandsArg.isSet()){
        xatu::SystemTB system = xatu::SystemTB(*systemConfig);
        system.setAU(dftArg.isSet());

        system.solveBands(kpointsfile);

        return 0;
    } else if (screeningArg.isSet() && !excitonArg.isSet()){
        throw std::invalid_argument("For screening, both screening and exciton files must be provided.");
    } else if (!screeningArg.isSet() && excitonArg.isSet()) {
        if (!excitonArg.isSet()){
            throw std::invalid_argument("Must provide exciton file.");
        }

        excitonConfig.reset(new xatu::ExcitonConfiguration(excitonfile));

    } else {
        screeningConfig.reset(new xatu::ScreeningConfiguration(screeningfile));
        excitonConfig.reset(new xatu::ExcitonConfiguration(excitonfile));
    }

    // -------------------------- Main body ---------------------------

    xatu::printHeader();

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                  Parameters                               |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    xatu::ExcitonTB bulkExciton = xatu::ExcitonTB(*systemConfig, *excitonConfig);
    bulkExciton.setMode(excitonConfig->excitonInfo.mode);
    bulkExciton.system->setAU(dftArg.isSet());

    if (screeningArg.isSet()){
        bulkExciton.initializeScreeningAttributes(*screeningConfig,excitonConfig->excitonInfo.mode);
    }

    bulkExciton.verifypotential();

    cout << std::left << std::setw(30) << "System configuration file: " << std::setw(10) << systemfile << endl;
    cout << std::left << std::setw(30) << "Exciton configuration file: " << std::setw(10) << excitonfile << endl;
    if (screeningArg.isSet()){
        cout << std::left << std::setw(30) << "Screening configuration file: " << std::setw(10) << screeningfile << "\n" << endl;
    } else {
        cout << std::left << std::setw(30) << "Screening configuration file: " << std::setw(10) << "Not provided.\n" << endl;
    }

    bulkExciton.printInformation();
    
    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                Initialization                             |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    if(excitonConfig->excitonInfo.submeshFactor != 1){
        bulkExciton.system->reducedBrillouinZoneMesh(excitonConfig->excitonInfo.ncell, excitonConfig->excitonInfo.submeshFactor);   
    }
    else{
        bulkExciton.brillouinZoneMesh(excitonConfig->excitonInfo.ncell);
    }

    if(!excitonConfig->excitonInfo.shift.is_empty()){
        bulkExciton.system->shiftBZ(excitonConfig->excitonInfo.shift);
    }


    bulkExciton.initializeHamiltonian();

    // If screening flag is present, computes the static dielectric function first
    if(screeningArg.isSet()){

        if (screeningConfig->screeningInfo.function == "dielectric") {
            
            std::complex<double> epsilon = bulkExciton.computesingleDielectricFunctionMatrixElement();            
            
            return 0;
        } else if (screeningConfig->screeningInfo.function == "polarizability") {

            if (bulkExciton.mode == "reciprocalspace")
            {
                bulkExciton.writeBZtofile();
            }

            bulkExciton.PolarizabilityMesh();

            return 0;

        } else if (screeningConfig->screeningInfo.function == "inversedielectric") {
            std::cout << "Proceeding with computation of the inverse of a dielectric matrix...\n" << std::endl;

            bulkExciton.computesingleInverseDielectricMatrix(excitonConfig->excitonInfo.label);

            return 0;

        } else if (screeningConfig->screeningInfo.function == "exciton"){
            
            std::cout << "Proceeding with computation of the permittivity matrix..." << std::endl;

            if (bulkExciton.mode == "reciprocalspace")
            {
                bulkExciton.writeBZtofile();
            }

            bulkExciton.compute_2D_DielectricMatrix();

            bulkExciton.invertDielectricMatrix();
            
            std::string filename_dielectric = excitonConfig->excitonInfo.label;

            bulkExciton.writeInverseDielectricMatrix(filename_dielectric + "_invepsilon.dat");

        } else {
            std::cout << "\nThe value '" + screeningConfig->screeningInfo.function + "' assigned to 'function' not recognized. Terminating program.\n" << std::endl;
            return 0;
        }
    }

    bulkExciton.BShamiltonian();
    auto results = bulkExciton.diagonalize(method, nstates);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Results                                |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    xatu::printEnergies(results, nstates, decimals);

    cout << "+---------------------------------------------------------------------------+" << endl;
    cout << "|                                    Output                                 |" << endl;
    cout << "+---------------------------------------------------------------------------+" << endl;

    std::string output = excitonConfig->excitonInfo.label;

    // --------------------------- Output ---------------------------
    bool writeEigvals = energyArg.isSet();
    if(writeEigvals){
        std::string filename_en = output + ".eigval";
        FILE* textfile_en = fopen(filename_en.c_str(), "w");

        std::cout << "Writing eigvals to file: " << filename_en << std::endl;
        fprintf(textfile_en, "%d\n", excitonConfig->excitonInfo.ncell);
        results->writeEigenvalues(textfile_en, nstates);

        fclose(textfile_en);
    }
    
    bool writeStates = eigenstatesArg.isSet();
    if(writeStates){
        std::string filename_st = output + ".states";
        FILE* textfile_st = fopen(filename_st.c_str(), "w");

        std::cout << "Writing states to file: " << filename_st << std::endl;
        results->writeStates(textfile_st, nstates);

        fclose(textfile_st);
    }
    
    bool writeWF = reciprocalArg.isSet();
    if(writeWF){
        std::string filename_kwf = output + ".kwf";
        FILE* textfile_kwf = fopen(filename_kwf.c_str(), "w");

        std::cout << "Writing k w.f. to file: " << filename_kwf << std::endl;
        for(int stateindex = 0; stateindex < nstates; stateindex++){
            if (excitonConfig->excitonInfo.submeshFactor != 1){
                results->writeReciprocalAmplitude(stateindex, textfile_kwf);
            }
            else{
                results->writeExtendedReciprocalAmplitude(stateindex, textfile_kwf);
            }
        }

        fclose(textfile_kwf);
    }
    
    bool writeRSWF = realspaceArg.isSet();
    if(writeRSWF){
        std::string filename_rswf = output + ".rswf";
        FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");

        arma::uvec statesToWrite = arma::regspace<arma::uvec>(0, nstates - 1);
        std::cout << "Writing real space w.f. to file: " << filename_rswf << std::endl;
        arma::rowvec holeCell = {0., 0., 0.};
        
        for(unsigned int i = 0; i < statesToWrite.n_elem; i++){
            std::cout << "Writing state " << i + 1 << " out of " << statesToWrite.n_elem << std::endl;
            results->writeRealspaceAmplitude(statesToWrite(i), holeIndex, holeCell, textfile_rswf, ncellsRSWF);
        }

        fclose(textfile_rswf);
    }

    bool writeAbs = absorptionArg.isSet();
    if(writeAbs){
        std::cout << "Writing absorption spectrum fo file... " << std::endl;
        results->writeAbsorptionSpectrum();
    }

    bool writeSpin = spinArg.isSet();
    if(writeSpin){
        std::string filename_spin = output + ".spin";
        FILE* textfile_spin = fopen(filename_spin.c_str(), "w");

        std::cout << "Writing excitons spin fo file: " << filename_spin << std::endl;
        results->writeSpin(nstates, textfile_spin);
    }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};