#include <armadillo>
#include <xatu.hpp>
#include <iostream>
#include <set>

using namespace std::chrono;

// run command: ./read_screening ../examples/material_models/DFT/hBN_base_HSE06.outp ../examples/excitonconfig/hBN_reciprocal.txt ../examples/screeningconfig/hBN_DFT_screening.txt <name_of_inv_dielectric_matrix_file>.dat

int main(int argc, char* argv[]){

    // Parse console stdin
    int n_args = 5;

    if (argc != n_args){
		throw std::invalid_argument("Error: Four input files are expected");
	}
    else if (argc < n_args){
        throw std::invalid_argument("Error: At least two input file are required (system config, exciton config, screening config and screening data file).");
    };

    std::string modelfile = argv[1];
    std::string excitonfile = argv[2];
    std::string screeningfile = argv[3];
    std::string inv_epsilon_file = argv[4];

    int nstates = 8;
    int decimals = 6;

    std::unique_ptr<xatu::SystemConfiguration> systemConfig;
    std::unique_ptr<xatu::ExcitonConfiguration> excitonConfig;
    std::unique_ptr<xatu::ScreeningConfiguration> screeningConfig;

    if (modelfile.find(".outp") != std::string::npos){
        systemConfig.reset(new xatu::CRYSTALConfiguration(modelfile, 100));
    }
    else if (modelfile.find(".model") != std::string::npos){
        systemConfig.reset(new xatu::SystemConfiguration(modelfile));
    } else {
        throw std::invalid_argument("Error: Unsupported system configuration file format. Use .outp or .model files.");
    }

    screeningConfig.reset(new xatu::ScreeningConfiguration(screeningfile));
    excitonConfig.reset(new xatu::ExcitonConfiguration(excitonfile));


    xatu::ExcitonTB exciton = xatu::ExcitonTB(*systemConfig, *excitonConfig, *screeningConfig);

    exciton.setMode(excitonConfig->excitonInfo.mode);
    
    if (modelfile.find(".outp") != std::string::npos){
        exciton.system->setAU(true); // if input model is CRYSTAL
    }

    exciton.brillouinZoneMesh(exciton.ncell);
    exciton.initializeHamiltonian();

    exciton.printInformation();

    exciton.readInverseDielectricMatrix(inv_epsilon_file);
       
    exciton.BShamiltonian();

    auto results = exciton.diagonalize("diag", nstates);

    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                    Results                                |" << std::endl;
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;

    xatu::printEnergies(results, nstates, decimals);
    
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                    Output                                 |" << std::endl;
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;

    std::string output = excitonConfig->excitonInfo.label;

    // --------------------------- Output ---------------------------
    bool writeEigvals = true;
    if(writeEigvals){
        std::string filename_en = output + ".eigval";
        FILE* textfile_en = fopen(filename_en.c_str(), "w");

        std::cout << "Writing eigvals to file: " << filename_en << std::endl;
        fprintf(textfile_en, "%d\n", excitonConfig->excitonInfo.ncell);
        results->writeEigenvalues(textfile_en, nstates);

        fclose(textfile_en);
    }
    
    bool writeStates = true;
    if(writeStates){
        std::string filename_st = output + ".states";
        FILE* textfile_st = fopen(filename_st.c_str(), "w");

        std::cout << "Writing states to file: " << filename_st << std::endl;
        results->writeStates(textfile_st, nstates);

        fclose(textfile_st);
    }
    
    bool writeWF = true;
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
    
    bool writeAbs = true;
    if(writeAbs){
        std::cout << "Writing absorption spectrum fo file... " << std::endl;
        results->writeAbsorptionSpectrum();
    }

    return 0;
}