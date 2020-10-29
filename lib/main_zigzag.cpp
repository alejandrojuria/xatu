#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "Zigzag.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){

    auto start = high_resolution_clock::now();

	// ------------------------ Model output into files ------------------------

	std::string filename_wf = "bands_TB_wavefunction";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
	bool writeTBwf = true;

	std::string filename = "bands_N4_k2001";
	FILE* textfile = fopen(filename.c_str(), "w");
	bool writeBands = true;

	std::string filename_pos = "bands_position_edge";
	FILE* textfile_pos = fopen(filename_pos.c_str(), "w");
	bool writeEdgeLocalization = false;

	std::string filename_over = "bands_overlap";
	FILE* textfile_over = fopen(filename_over.c_str(), "w");
	bool writeOverlap = false;

	std::string filename_der = "bands_derivative";
	FILE* textfile_der = fopen(filename_der.c_str(), "w");
	bool writeBandDerivative = false;

	std::string filename_spin = "bands_spin_4bands";
	FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
	bool writeSpin = true;

	std::string filename_dos = "bands_dos_bulk";
	FILE* textfile_dos = fopen(filename_dos.c_str(), "w");
	bool writeDoS = false;

	// ------------------------ Initialization ------------------------

    int N = 15;
    Zigzag system = Zigzag(N);

	int Ncell = 300;
    int nk = 2*Ncell; // Avoid k=PI/a for edge stats
	//double zeeCenter = 0.693268 - PI/a;
	vec kpoints = arma::linspace(0, 2*PI/system.a, nk);
	int Qindex = 0;
	double Q = 0.0;

	// Indices for file writing
	int tbWFindex = 53;
	int localizationIndex = nk/2;

    cx_mat eigenvec, eigenvecQ;
	vec eigenval, eigenvalQ;
	int vband = 2*(N+1)*5-2;
	int vband2 = 2*(N+1)*5-1;
	int cband = 2*(N+1)*5+0;
	int cband2 = 2*(N+1)*5+1;

	bool calculateMaxValence = false;
	// Variables for valence maximum calculation
	double maxEnergy, maxK;

	// Variables for overlap calculation 
	cx_vec coefv1, coefv2;
	cx_vec coefc1, coefc2;
	double overlapv, overlapc, directTerm, exchangeTerm;

	// Variables for derivative calculation
	double derivative, prevEnergy;

	// Variable for DoS calculation
	int excludedBands = 2*(N+1)*8 - 4;
	int dimDoS = 2*(N+1)*8 - excludedBands;
	mat energies(dimDoS, nk);
	double FermiEnergy;

	// ---------------------------- Main loop ----------------------------
	for (int i = 0; i < nk; i++) {

		cx_mat h = system.hamiltonian(kpoints(i));
		arma::eig_sym(eigenval, eigenvec, h);

		cout << "k = : " << kpoints(i) << endl;
		cout << "Q = " << 0 << endl;
		double spinV = system.expectedSpinZValue(eigenvec.col(vband));
		double spinC = system.expectedSpinZValue(eigenvec.col(cband + 1));
		double spinV_1 = system.expectedSpinZValue(eigenvec.col(vband + 1));
		double spinC_1 = system.expectedSpinZValue(eigenvec.col(cband));
		cout << "Spin valence: " << spinV << endl;
		cout << "Spin valence deg.: " << spinV_1 << endl;
		cout << "Spin conduction: " << spinC << endl;
		cout << "Spin conduction deg.: " << spinC_1 << endl;

		if(i == (int)nk/2){
			FermiEnergy = eigenval(vband + 1);
		}

		if(writeDoS){
			//vec eigenval1 = eigenval(arma::span(0, vband - 1));
			//vec eigenval2 = eigenval(arma::span(cband + 1, 2*(N+1)*8-1));
			//eigenval = arma::join_cols(eigenval1, eigenval2);
			eigenval = eigenval(arma::span(vband, cband));

			energies.col(i) = eigenval;
		};

		if(Q != 0){
			cx_mat hQ = system.hamiltonian(kpoints(i) + Q);
			arma::eig_sym(eigenvalQ, eigenvecQ, hQ);

			cout << "Q = :" << Q << endl;
			double spinV = system.expectedSpinZValue(eigenvecQ.col(vband - 1));
			double spinC = system.expectedSpinZValue(eigenvecQ.col(cband + 1));
			double spinV_1 = system.expectedSpinZValue(eigenvecQ.col(vband));
			double spinC_1 = system.expectedSpinZValue(eigenvecQ.col(cband));
			cout << "Spin valence: " << spinV << endl;
			cout << "Spin valence deg.: " << spinV_1 << endl;
			cout << "Spin conduction: " << spinC << endl;
			cout << "Spin conduction deg.: " << spinC_1 << endl;
		};

		if(calculateMaxValence == true){
			if(i == 0){
				maxEnergy = eigenval(vband);
			}
			else if(eigenval(vband) > maxEnergy){
				maxEnergy = eigenval(vband);
				maxK = kpoints(i);
			};
		};

		// ----------------- Write TB real w.f. into file -----------------
		// Specific (n,k) state
		if(writeTBwf && (i == tbWFindex)){
			int band = 2*(N+1)*5 + 1;
			for(int j = 0; j < (int)system.motif.n_rows; j++){

				std::complex<double> coef = 0;
				rowvec position = system.motif.row(j);
				for(int m = 0; m < 8; m++){
					coef += abs(eigenvec.col(band)(8*j+m))*abs(eigenvec.col(band)(8*j+m));
				};
				fprintf(textfile_wf, "%lf\t%lf\t%lf\n", position(0), 
						position(1), coef);
			};
		};

		// ----------------- Write TB overlap into file -----------------
		// k-dependent (band)
		if(writeOverlap){
			if(i == 0){
				coefv1 = eigenvec.col(vband);
				coefv2 = eigenvec.col(vband);
				coefc1 = eigenvec.col(cband);
				coefc2 = eigenvec.col(cband);
				if(Q != 0){
					coefc1 = eigenvecQ.col(cband);
					coefc2 = eigenvecQ.col(cband);
				}
			}
			coefv2 = eigenvec.col(vband);
			if(Q != 0){
				coefc2 = eigenvecQ.col(cband);
			}
			else{
				coefc2 = eigenvec.col(cband);
			}

			overlapv = pow(abs(arma::cdot(coefv1, coefv2)),2);
			overlapc = pow(abs(arma::cdot(coefc1, coefc2)),2);

			directTerm = abs(arma::cdot(coefv1, coefv2)*arma::cdot(coefc1, coefc2));
			exchangeTerm = abs(arma::cdot(coefc1, coefv1)*arma::cdot(coefv2, coefc2));

			fprintf(textfile_over, "%lf\t%lf\t%lf\t%lf\t%lf\n", 
					kpoints(i), overlapv, overlapc, directTerm, exchangeTerm);
		};
		
		// -------------------- Calculate & write state localization --------------------
		// k-dependent (band)
		if(writeEdgeLocalization){
			double coefv = 0;
			double coefc = 0;

			// ----- Select right or left edges -----
			//int j = motif.n_rows - 1;
			//int j2 = motif.n_rows - 2;
			int j = 0;
			int j2 = 1;
			
			for(int m = 0; m < 8; m++){
				coefv += abs(eigenvec.col(vband)(8*j+m))*abs(eigenvec.col(vband)(8*j+m)) + 
				abs(eigenvec.col(vband)(8*j2+m))*abs(eigenvec.col(vband)(8*j2+m));
				coefc += abs(eigenvec.col(cband)(8*j+m))*abs(eigenvec.col(cband)(8*j+m)) + 
				abs(eigenvec.col(cband)(8*j2+m))*abs(eigenvec.col(cband)(8*j2+m));
			};
			fprintf(textfile_pos, "%lf\t%lf\n", coefv, coefc);
		};

		// -------------------- Calculate & write expected spin value (components) --------------------
		// k-dependent (band)
		if(writeSpin){
			double szV = system.expectedSpinZValue(eigenvec.col(vband));
			double szV_2 = system.expectedSpinZValue(eigenvec.col(vband2));
			double szC = system.expectedSpinZValue(eigenvec.col(cband));
			double szC_2 = system.expectedSpinZValue(eigenvec.col(cband2));
			double syV = system.expectedSpinYValue(eigenvec.col(vband));
			double syV_2 = system.expectedSpinYValue(eigenvec.col(vband2));
			double syC = system.expectedSpinYValue(eigenvec.col(cband));
			double syC_2 = system.expectedSpinYValue(eigenvec.col(cband2));			
			double sxV = system.expectedSpinXValue(eigenvec.col(vband));
			double sxV_2 = system.expectedSpinXValue(eigenvec.col(vband2));
			double sxC = system.expectedSpinXValue(eigenvec.col(cband));
			double sxC_2 = system.expectedSpinXValue(eigenvec.col(cband2));
			fprintf(textfile_spin, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					kpoints(i), szV, szV_2, szC, szC_2, syV, syV_2, syC, syC_2, sxV, sxV_2, sxC, sxC_2);
		};


		// ----------------------- Compute derivative for a certain band -----------------------
		// k-dependent (band)
		if(writeBandDerivative){
			if(i==0){
				prevEnergy = eigenval(vband);
				continue;
			};
			derivative = (eigenval(vband) - prevEnergy)/(kpoints(i) - kpoints(i-1));

			fprintf(textfile_der, "%lf\t%lf\n", kpoints(i-1), derivative);
			prevEnergy = eigenval(vband);
		};

		if(writeBands){
			system.writeEigenvaluesToFile(textfile, eigenval, kpoints[i]);
		};
	};

	// Output results after main loop
	if(calculateMaxValence){
		cout << "Valence band maximum (edge)" << endl;
		cout << "k points (A-1) : " << maxK << endl;
		cout << "Energy (eV): " << maxEnergy << endl;
	};

	if(writeDoS){
		energies = energies - FermiEnergy;
		double delta = 0.01;
		system.writeDensityOfStates(energies, delta, textfile_dos);
	};

	// Close all text files
	fclose(textfile);
	fclose(textfile_wf);
	fclose(textfile_pos);
	fclose(textfile_over);
	fclose(textfile_der);
	fclose(textfile_spin);
	fclose(textfile_dos);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << duration.count() << " ms" << std::endl;

	return 0;
};