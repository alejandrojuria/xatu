#include <complex>
#include <math.h>
#include <stdlib.h>
#include "xatu/System.hpp"

namespace xatu {

// -------------------- Constructors and destructor --------------------

/**
 * Copy constructor.
 * @details Initializes one default System objects, and afterwards copies
 * its attributes to be the same as those of the given System object to copy.
 * @param system SystemTB object to copy.
 */
System::System(const System& system) : Lattice(system){

	systemName			  = system.systemName;
	orbitalsPerSpecies_   = system.orbitalsPerSpecies;
	hamiltonianMatrices_  = system.hamiltonianMatrices;
	overlapMatrices_      = system.overlapMatrices;
	filling_			  = system.filling;
	fermiLevel_			  = filling_ - 1;
	basisdim_ 			  = system.basisdim;
}

/**
 * Configuration constructor.
 * @details Constructor which takes in a SystemConfiguration object, i.e.
 * to init a System from a configuration file.
 * @param configuration SystemConfiguration object obtained from config. file.
 */
System::System(const SystemConfiguration& configuration) : Lattice(){

	initializeLatticeAttributes(configuration);
	initializeSystemAttributes(configuration);
};

// --------------------------------- Methods ---------------------------------
/**
 * Routine to extract the information contained in the SystemConfiguration object from
 * parsing the input text file. To be used in the TB mode.
 * @details To be called from the configuration constructor.
 * @param configuration SystemConfiguration object.
 */
void System::initializeSystemAttributes(const SystemConfiguration& configuration){
	
	orbitalsPerSpecies_   = configuration.systemInfo.orbitalsPerSpecies;
	hamiltonianMatrices_  = configuration.systemInfo.hamiltonian;
	overlapMatrices_      = configuration.systemInfo.overlap;
	filling_			  = configuration.systemInfo.filling;
	fermiLevel_			  = filling_ - 1;

    int basisdim = 0;
    for(int i = 0; i < natoms; i++){
        int species = this->motif.row(i)(3);
        basisdim += orbitalsPerSpecies(species);
    }
	basisdim_   = basisdim;
}


/**
 * Filling setter.
 * @details Sets both the filling and the Fermi level, which is defined as filling - 1.
 * Filling must be a positive integer.
 * @param filling Number of electrons of the system.
*/
void System::setFilling(int filling){
	if (filling > 0){
		filling_ = filling;
		fermiLevel_ = filling_ - 1;
	}
	else{
		std::cout << "Filling must be a positive integer" << std::endl;
	}
}


/**
 * System name setter.
 * @details Sets the name of the system.
 * @param name Name of the system.
*/
void System::setSystemName(std::string name){
	systemName = name;
}


/**
 * Method to write to a file the energy bands evaluated on a set of kpoints specified on a file.
 * @details It creates a file with the name "[systemName].bands" where the bands are stores.
 * @param kpointsfile File with the kpoints where we want to obtain the bands.
 * @param triangular To specify if the Hamiltonian is triangular.
*/
void System::solveBands(std::string kpointsfile) const {
	std::ifstream inputfile;
	std::string line;
	double kx, ky, kz;
	arma::vec eigval;
	arma::cx_mat eigvec;
	std::string outputfilename = kpointsfile + ".bands";
	FILE* bandfile = fopen(outputfilename.c_str(), "w");
	try{
		inputfile.open(kpointsfile.c_str());
		while(std::getline(inputfile, line)){
			std::istringstream iss(line);
			iss >> kx >> ky >> kz;
			arma::rowvec kpoint{kx, ky, kz};
			solveBands(kpoint, eigval, eigvec);
			for (int i = 0; i < eigval.n_elem; i++){
				fprintf(bandfile, "%12.6f\t", eigval(i));
			}
			fprintf(bandfile, "\n");
		}
	}
	catch(const std::exception& e){
		std::cerr << e.what() << std::endl;
	}
	fclose(bandfile);
	arma::cout << "Done" << arma::endl;
}

/**
 * Method to add a Zeeman term to the Hamiltonian. 
 * @details This method assumes that the Hamiltonian incorporates spin in 
 * the following way: |i1,up>,|i1,down>,...,|in,up>,|in,down>, where i runs over orbitals.
 * @param amplitude Strength of the Zeeman term.
 * @returns void
*/
void System::addZeeman(double amplitude){

	arma::cx_vec zeeman_values = {amplitude, -amplitude};
	arma::cx_mat zeeman_matrix = arma::diagmat<arma::cx_mat>(arma::kron(arma::ones<arma::cx_vec>(basisdim/2), zeeman_values));

	// Identify hamiltonian slice for R=0
	int idx;
	for (unsigned int i = 0; i < unitCellList.n_rows; i++){
		if (arma::norm(unitCellList.row(i)) < 1E-5){
			arma::cout << unitCellList.row(i) << arma::endl;
			idx = i;
			break;
		}
	}

	hamiltonianMatrices_.slice(idx) += zeeman_matrix;
}

}