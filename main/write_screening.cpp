#include <armadillo>
#include <xatu.hpp>
#include <iostream>
#include <set>

using namespace std::chrono;

// run command: ./write_screening ../examples/material_models/DFT/hBN_base_HSE06.outp ../examples/excitonconfig/hBN_test.txt ../examples/screeningconfig/hBN_DFT_screening.txt <name_of_q_points_file>.dat

int main(int argc, char* argv[]){

    // Parse console stdin

    if (argc != 5){
		throw std::invalid_argument("Error: Four input files are expected");
	}
    else if (argc < 5){
        throw std::invalid_argument("Error: At least two input file are required (system config, exciton config, screening config and q points file).");
    };

    std::string modelfile = argv[1];
    std::string excitonfile = argv[2];
    std::string screeningfile = argv[3];
    std::string q_points_file = argv[4];
    
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
    
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    std::cout << "|                                  Parameters                               |" << std::endl;
    std::cout << "+---------------------------------------------------------------------------+" << std::endl;
    
    std::cout << "System configuration file: " << modelfile << "\n" << std::endl;
    std::cout << "q points file file:        " << q_points_file << "\n" << std::endl;


    exciton.brillouinZoneMesh(exciton.ncell);
    exciton.initializeHamiltonian();    

    exciton.compute_2D_DielectricMatrix(q_points_file);

    size_t lastindex = q_points_file.find_last_of("."); 
    std::string rawname = q_points_file.substr(0, lastindex); 

    exciton.invertDielectricMatrix();    

    exciton.writeInverseDielectricMatrix(rawname + "_invepsilon.dat");

    return 0;
}