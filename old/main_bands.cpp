#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "zigzag.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

// Shared variable declaration
double a, c;
double Es, Ep, Vsss, Vsps, Vpps, Vppp;
double lambda, zeeman, onsiteEdge;
vec a1, a2, tau;
vec n1, n2, n3;
vec Gamma, K, M;
mat M0, M1, M2p, M2m, Mzeeman;
cx_mat Mso;
cx_mat H0, Ha, Hsoc, Hzeeman;
cx_mat spinOperator;

cx_mat HBS;
mat HK;
mat states;


// Header for band crossing
bool detectBandCrossing(vec, vec, vec, int, int, vec);

/* Routine to fix band crossing for two given bands. Compares slope of one band against the other to 
detect abrupt changes. Thus it measures \delta(actual_slope - previous_slope)*/
bool detectBandCrossing(vec actualEigenval, vec previousEigenval, vec previousEigenval2, int firstBand, int secondBand, vec kpoints){

	double initialSlope = (previousEigenval(firstBand) - previousEigenval2(firstBand))/(kpoints(1) - kpoints(0));
	double actualSlopeFirst = (actualEigenval(firstBand) - previousEigenval(firstBand))/(kpoints(1) - kpoints(0));
	double actualSlopeSecond = (actualEigenval(secondBand) - previousEigenval(firstBand))/(kpoints(1) - kpoints(0));
	cout << "first slope: " << initialSlope << endl; 
	cout << "second slope: " << actualSlopeFirst << endl; 
	cout << "second slope p: " << actualSlopeSecond << endl; 


	double deltaFirst = abs((actualSlopeFirst - initialSlope)/initialSlope);
	double deltaSecond = abs((actualSlopeSecond - initialSlope)/initialSlope);
	
	if (deltaSecond < deltaFirst){
		return true;
	}
	else{
		return false;
	};
};

int main(){

	vec previousEigenvalues;
	vec previousEigenvalues2;
	vec auxEigval;
	cx_mat auxEigvec;

    auto start = high_resolution_clock::now();

	// ------------------------ Model output into files ------------------------

	std::string filename_wf = "bands_TB_wavefunction";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
	bool writeTBwf = false;

	std::string filename = "bands_N15_k2001";
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

	std::string filename_spin = "bands_spin";
	FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
	bool writeSpin = true;

	std::string filename_dos = "bands_dos_bulk";
	FILE* textfile_dos = fopen(filename_dos.c_str(), "w");
	bool writeDoS = false;

	// ------------------------ Initialization ------------------------

    initializeConstants();
    //lambda = 0;
	// Number of cells in the finite direction
	int N = 15;

	initializeBlockMatrices();
	prepareHamiltonian(N);

	int Ncell = 1000;
    int nk = 2*Ncell + 1; // Avoid k=PI/a for edge stats
	//double zeeCenter = 0.693268 - PI/a;
	vec kpoints = arma::linspace(0.0, 2*PI/a + 0.0, nk);

	// Fine mesh to resolve gap with infinitesimal zeeman
	bool useFineMesh = false;
	if(useFineMesh){
		int nfine = 2000;
		vec kpoints_fine = arma::linspace(PI/a-PI/(500*a), PI/a+PI/(500*a), nfine);
		kpoints = arma::sort(arma::join_cols(kpoints, kpoints_fine));
		nk = kpoints.n_elem;
	};
	int Qindex = 0;
	double Q = 0.0;

	// Indices for file writing
	int tbWFindex = nk/2;
	int localizationIndex = nk/2;

    cx_mat eigenvec, eigenvecQ;
	vec eigenval, eigenvalQ;
	mat motif = createMotiv(N);
	int vband = 2*(N+1)*5-4;
	int vband2 = 2*(N+1)*5-3;
	int cband = 2*(N+1)*5+2;
	int cband2 = 2*(N+1)*5+3;

	bool doSwitchBandCross = false;
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

	int count = 0;

	// ---------------------------- Main loop ----------------------------
	for (int i = 0; i < nk; i++) {

		cx_mat h = hamiltonian(kpoints(i), H0, Ha, Hsoc);
		arma::eig_sym(eigenval, eigenvec, h);

		cout << "k = " << kpoints(i) << endl;
		cout << "Q = " << 0 << endl;

		if(i != 0 && i != 1 && i == -1){
			bool detectV1 = detectBandCrossing(eigenval, previousEigenvalues, previousEigenvalues2, vband, vband - 2, kpoints);
			bool detectV2 = detectBandCrossing(eigenval, previousEigenvalues, previousEigenvalues2, vband + 1, vband - 1, kpoints);
			bool detectC1 = detectBandCrossing(eigenval, previousEigenvalues, previousEigenvalues2, cband, cband + 2, kpoints); 
			bool detectC2 = detectBandCrossing(eigenval, previousEigenvalues, previousEigenvalues2, cband + 1, cband + 3, kpoints); 

			if(detectV1){
				cout << "Crossing!" << endl;
				cout << "Counts: " << count << endl;
				auxEigval = eigenval;
				eigenval(vband) = eigenval(vband - 2);
				eigenval(vband - 2) = auxEigval(vband);

				auxEigvec = eigenvec;
				eigenvec.col(vband) = eigenvec.col(vband - 2);
				eigenvec.col(vband - 2) = auxEigvec.col(vband);
				count++;
			}
			if(detectV2){
				auxEigval = eigenval;
				eigenval(vband + 1) = eigenval(vband - 1);
				eigenval(vband - 1) = auxEigval(vband + 1);

				auxEigvec = eigenvec;
				eigenvec.col(vband + 1) = eigenvec.col(vband - 1);
				eigenvec.col(vband - 1) = auxEigvec.col(vband + 1);
			}
			if(detectC1){
				auxEigval = eigenval;
				eigenval(cband) = eigenval(cband + 2);
				eigenval(cband + 2) = auxEigval(cband);

				auxEigvec = eigenvec;
				eigenvec.col(cband) = eigenvec.col(cband + 2);
				eigenvec.col(cband + 2) = auxEigvec.col(cband);
			}
			if(detectC2){
				auxEigval = eigenval;
				eigenval(cband + 1) = eigenval(cband + 3);
				eigenval(cband + 3) = auxEigval(cband + 1);

				auxEigvec = eigenvec;
				eigenvec.col(cband + 1) = eigenvec.col(cband  + 3);
				eigenvec.col(cband + 3) = auxEigvec.col(cband + 1);
			}
		};

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
			cx_mat hQ = hamiltonian(kpoints(i) + Q, H0, Ha, Hsoc);
			arma::eig_sym(eigenvalQ, eigenvecQ, hQ);

			cout << "Q = :" << Q << endl;
			double spinV = expectedSpinZValue(eigenvecQ.col(vband - 1), N);
			double spinC = expectedSpinZValue(eigenvecQ.col(cband + 1), N);
			double spinV_1 = expectedSpinZValue(eigenvecQ.col(vband), N);
			double spinC_1 = expectedSpinZValue(eigenvecQ.col(cband), N);
			cout << "Spin valence: " << spinV << endl;
			cout << "Spin valence deg.: " << spinV_1 << endl;
			cout << "Spin conduction: " << spinC << endl;
			cout << "Spin conduction deg.: " << spinC_1 << endl;
		};

		if(doSwitchBandCross && (kpoints(i) > PI/a) && (kpoints(i) < 2*PI/a)){
			// Valence band
			double auxEigenval = eigenval(vband);
			eigenval(vband) = eigenval(vband + 1);
			eigenval(vband + 1) = auxEigenval;
			cx_vec auxEigenvec = eigenvec.col(vband);
			eigenvec.col(vband) = eigenvec.col(vband + 1);
			eigenvec.col(vband + 1) = auxEigenvec;

			// Conduction band
			auxEigenval = eigenval(cband);
			eigenval(cband) = eigenval(cband - 1);
			eigenval(cband - 1) = auxEigenval;
			auxEigenvec = eigenvec.col(cband);
			eigenvec.col(cband) = eigenvec.col(cband - 1);
			eigenvec.col(cband - 1) = auxEigenvec;

			if(Q != 0){
				auxEigenvec = eigenvecQ.col(vband);
				eigenvecQ.col(vband) = eigenvecQ.col(vband + 1);
				eigenvecQ.col(vband + 1) = auxEigenvec;

				auxEigenvec = eigenvecQ.col(cband);
				eigenvecQ.col(cband) = eigenvecQ.col(cband - 1);
				eigenvecQ.col(cband - 1) = auxEigenvec;
			};
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
			int band = cband + 1;
			for(int j = 0; j < (int)motif.n_rows; j++){

				std::complex<double> coef = 0;
				rowvec position = motif.row(j);
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
			double szV = expectedSpinZValue(eigenvec.col(vband), N);
			double szV_2 = expectedSpinZValue(eigenvec.col(vband2), N);
			double szV_3 = expectedSpinZValue(eigenvec.col(vband - 1), N);
			double szV_4 = expectedSpinZValue(eigenvec.col(vband - 2), N);
			double sxV = expectedSpinXValue(eigenvec.col(vband), N);
			double sxV_2 = expectedSpinXValue(eigenvec.col(vband2), N);
			double sxV_3 = expectedSpinXValue(eigenvec.col(vband - 1), N);
			double sxV_4 = expectedSpinXValue(eigenvec.col(vband - 2), N);
			double syV = expectedSpinYValue(eigenvec.col(vband), N);
			double syV_2 = expectedSpinYValue(eigenvec.col(vband2), N);
			double syV_3 = expectedSpinYValue(eigenvec.col(vband - 1), N);
			double syV_4 = expectedSpinYValue(eigenvec.col(vband - 2), N);
			fprintf(textfile_spin, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					kpoints(i), szV, szV_2, szV_3, szV_4, sxV, sxV_2, sxV_3, sxV_4, syV, syV_2, syV_3, syV_4);
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
			writeEigenvaluesToFile(textfile, eigenval, kpoints[i]);
		};

		if(i > 0){
			previousEigenvalues2 = previousEigenvalues;
		}
		previousEigenvalues = eigenval;
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
		writeDensityOfStates(energies, delta, textfile_dos);
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