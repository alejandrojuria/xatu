#include "xatu/System.hpp"

namespace xatu {

/**
 * Default constructor.
 * @details The default constructor throws an error as the class must always be initialized from either a ConfigurationSystem 
 * or a ConfigurationCRYSTAL object, or with a copy constructor.
 */
System::System(){
	
    throw std::invalid_argument("System must be called with either a ConfigurationSystem, a ConfigurationCRYSTAL or another System object");

}

/**
 * Constructor from ConfigurationSystem, to be used in the TB mode only.
 * @param SystemConfig ConfigurationSystem object obtained from any configuration file.
 */
System::System(const ConfigurationSystem& SystemConfig) : Lattice{SystemConfig}{

	this->highestValenceBand_	    = SystemConfig.filling - 1;
	this->Rlist_                    = SystemConfig.Rlist;
	this->ncells_                   = SystemConfig.ncells;
	this->norbitals_   			    = SystemConfig.hamiltonianMatrices.n_cols;
	this->ptr_hamiltonianMatrices   = &SystemConfig.hamiltonianMatrices;
	this->ptr_overlapMatrices       = &SystemConfig.overlapMatrices;
	
}

/**
 * Constructor from ConfigurationCRYSTAL, to be used in CRYSTAL calculations (compatible with both modes). 
 * @param CRYSTALconfig ConfigurationCRYSTAL object obtained from a CRYSTAL .outp file.
 */
System::System(const ConfigurationCRYSTAL& CRYSTALconfig) : Lattice{CRYSTALconfig}{

	this->highestValenceBand_   = CRYSTALconfig.filling - 1;
	this->Rlist_                = CRYSTALconfig.Rlist;
	this->ncells_               = CRYSTALconfig.ncells;
	this->norbitals_            = CRYSTALconfig.norbitals;
	this->ptr_overlapMatrices   = &CRYSTALconfig.overlapMatrices;
	
	if(!CRYSTALconfig.MAGNETIC_FLAG || CRYSTALconfig.SOC_FLAG){
		this->ptr_hamiltonianMatrices = &CRYSTALconfig.hamiltonianMatrices;
		this->ptr_alphaMatrices       = nullptr;
		this->ptr_betaMatrices        = nullptr;
	} 
	else { // CRYSTALconfig.MAGNETIC_FLAG && !CRYSTALconfig.SOC_FLAG (magnetism without SOC)
		if(CRYSTALconfig.overlapMatrices.n_cols == CRYSTALconfig.alphaMatrices.n_cols){ // <=> GTF mode
			this->ptr_hamiltonianMatrices = nullptr;
			this->ptr_alphaMatrices       = &CRYSTALconfig.alphaMatrices;
			this->ptr_betaMatrices        = &CRYSTALconfig.betaMatrices;
		} 
		else { // <=> TB mode
			this->ptr_hamiltonianMatrices = &CRYSTALconfig.hamiltonianMatrices;
			this->ptr_alphaMatrices       = nullptr;
			this->ptr_betaMatrices        = nullptr;
		}
	}

}

/**
 * System name setter.
 * @param systemName Name of the system.
 * @return void.
*/
void System::setSystemName(const std::string& systemName){
	
	this->systemName = systemName;

}

/**
 * Method to write to a file the energy bands evaluated on a set of kpoints specified on a file. This method is just a general envelope, 
 * and the core diagonalization method must be defined in the relevant derived class, depending on the mode (TB or Gaussian).
 * @details Each k-point must occupy a row in the file, with no blank lines. The format for each k-point is: kx ky kz ,
 * whose values are in Angstrom. It creates a file with the name "kpointsfile.bands" where the bands are stored.
 * @param kpointsfile Name of the file with the kpoints where we want to obtain the bands.
 * @return void.
*/
void System::solveBands(const std::string& kpointsfile) const {
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
			for (uint i = 0; i < eigval.n_elem; i++){
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


}