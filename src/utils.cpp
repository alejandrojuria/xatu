#include <armadillo>
#include <fstream>
#include <complex>


#include "utils.hpp"

void writeVectorsToFile(const arma::mat& vectors, FILE* textfile){
    for(unsigned int i = 0; i < vectors.n_rows; i++){
        arma::rowvec vector = vectors.row(i);
        fprintf(textfile, "%10.7lf\t%10.7lf\t%10.7lf\n", vector(0), vector(1), vector(2));
    }
}

arma::vec readVectorFromFile(std::string filename){
    std::ifstream file(filename.c_str());
    std::string line;
    std::vector<double> vector;
    double value;
    while(std::getline(file, line)){
        std::istringstream iss(line);
        while(iss >> value){
            vector.push_back(value);
        };
    };

    arma::vec coefs(vector);
    return coefs;
};

/* Definition of non-interacting retarded Green function */
std::complex<double> rGreenF(double energy, double delta, double eigEn){

	std::complex<double> i(0,1);
	return 1./((energy + i*delta) - eigEn);
};

/* Routine to calcule the density of states at a given energy,
associated to a given set of eigenvalues (e.g. bulk or edge).
NB: It is NOT normalized. */
double densityOfStates(double energy, double delta, const arma::mat& energies){

		double dos = 0;
		for(int i = 0; i < (int)energies.n_rows; i++){
			for(int j = 0; j < (int)energies.n_cols; j++){
				double eigEn = energies(i,j);
				dos += -PI*imag(rGreenF(energy, delta, eigEn));
			};
		};
        // Divide by number of k's and length a (currently a is missing)
		dos /= energies.n_cols; 

		return dos;
}

/* Routine to calculate and write the density of states associated 
to a given set of eigenenergies. Density of states is normalized to 1
(integral of DOS over energies equals 1)
Input: mat energies (eigenvalues), double delta (convergence parameter),
FILE* dosfile output file
Output: void (write results to output file) */
void writeDensityOfStates(const arma::mat& energies, double delta, FILE* dosfile){

	double minE = energies.min();
	double maxE = energies.max();
	int nE = 2000; // Number of points in energy mesh

	arma::vec energyMesh = arma::linspace(minE - 0.5, maxE + 0.5, nE);

	// Main loop
	double totalDos = 0;
	double dos = 0;
	// First loop over energies to normalize
	for(int n = 0; n < (int)energyMesh.n_elem; n++){
		double energy = energyMesh(n);
		double deltaE = energyMesh(1) - energyMesh(0);

		dos = densityOfStates(energy, delta, energies);
		totalDos += dos*deltaE;
	};
	// Main loop to write normalized DOS
	for(int n = 0; n < (int)energyMesh.n_elem; n++){
		double dos = 0;
		double energy = energyMesh(n);

		dos = densityOfStates(energy, delta, energies);
		dos = dos/totalDos; // Normallise
		fprintf(dosfile, "%lf\t%lf\n", energy, dos);
	};
	return;
};