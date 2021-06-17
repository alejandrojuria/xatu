#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "Exciton.hpp"
#include "Zigzag.hpp"
#include "libwavefunction_old.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){

    auto start = high_resolution_clock::now();
    
    int nBulkBands = 2;
    int nEdgeBands = 0;
    int N = 15;
    //NcellArray = arma::join_cols(NcellArray, arma::regspace(80, 20, 200));
    Zigzag ribbon = Zigzag(N);
    vec QArray = {0.0};

    cout << "c: " << ribbon.c << endl;
    cout << "a: " << ribbon.a << endl;

    // ----------------- Model parameters & Output --------------------

    std::string filename_en = "exciton_spectrum_exact";
    FILE* textfile_en = fopen(filename_en.c_str(), "w");
    bool writeEigvals = false;

    std::string filename_rswf = "exciton_bulk_wf_N15_nosoc";
    FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");
    bool writeRSWF = true;

    std::string filename_kwf = "exciton_bulk_wf_N15";
    FILE* textfile_kwf = fopen(filename_kwf.c_str(), "w");
    bool writeKWF = true;

    // ----------------------- Main body -----------------------

    fprintf(textfile_en, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");
    int Ncell = 15;
    int nk = 2*Ncell + 1;
    double Q = QArray(0);
    cx_mat eigvecX;
    vec eigvalX;

    cout << "##############################################################################" << endl;
    cout << "--------------------------------- Parameters ---------------------------------" << endl;
    cout << "N. cells: " << Ncell << endl;
    cout << "N: " << N << endl;
    cout << "nk = " << nk << endl;
    cout << "#bulk bands: " << nBulkBands << endl;
    cout << "#edge bands: " << nEdgeBands << "\n" << endl;

    Exciton bulkExciton = Exciton(N, Ncell, Q, nBulkBands, nEdgeBands);
    cout << "Valence bands:" << bulkExciton.valenceBands << endl;
    cout << "Conduction bands:" << bulkExciton.conductionBands << endl;
    //cout << bulkExciton.kpoints << endl;
    //cout << bulkExciton.potentialMat << endl;
    bulkExciton.BShamiltonian();
    arma::eig_sym(eigvalX, eigvecX, bulkExciton.HBS);

    cout << "#----------------------- Energy -----------------------" << endl;
    cout << "Fundamental state energy is: " << eigvalX(0) << endl;
    cout << "Second state energy is: " << eigvalX(1) << endl;
    cout << "Third state energy is: " << eigvalX(2) << endl;
    vec energies = bulkExciton.computeEnergies(eigvecX.col(0));
    cout << "Kinectic energy: "  << energies(0) << endl;
    cout << "Potential energy: " << energies(1) << "\n" << endl;

    if(writeEigvals){
        //writeEigenvaluesToFile(textfile, eigvalX, Q);
        fprintf(textfile_en, "%10.7lf\t", Q);
        for(int i = 0; i < eigvalX.n_elem; i++){
            fprintf(textfile_en, "%10.7lf\t", eigvalX(i));
        }
        fprintf(textfile_en, "\n");
    };

    //cx_mat states = bulkExciton.fixDegeneracy(eigvecX.col(0), eigvecX.col(1), 20);
    //cx_vec state = states.col(0);

    if(writeKWF){
        fprintf(textfile_kwf, "%d\n", N);
        cx_vec state = eigvecX.col(0);
        int nbands = bulkExciton.nBulkBands + bulkExciton.nEdgeBands;
        int nbands2 = nbands*nbands;
        for (int i = 0; i < nk - 1; i++){
            double coef = 0;
            for(int nband = 0; nband < nbands2; nband++){
                coef += abs(state(nbands2*i + nband))*abs(state(nbands2*i + nband));
                //coef += abs(eigvecX.col(state + 1)(nbands2*i + nband))*abs(eigvecX.col(state + 1)(nbands2*i + nband));
            };
            coef /= (bulkExciton.kpoints(1) - bulkExciton.kpoints(0)); // L2 norm instead of l2
            fprintf(textfile_kwf, "%11.8lf\t%11.8lf\n", bulkExciton.kpoints(i), coef);
        };
        fprintf(textfile_kwf, "#\n");
    };

    // Print real-space wavefunction
    if(writeRSWF){
        for(int i = 0; i < 1; i++){

            cx_vec excitonCoefficients = eigvecX.col(0);

            int holeIndex = bulkExciton.motif.n_rows/2 - 1;
            double holeCell = 0;
            rowvec holePosition = bulkExciton.motif.row(holeIndex) + holeCell*arma::rowvec{0., ribbon.a, 0.};

            cx_mat RScoefs = RScoefficientCalc(bulkExciton, 
                                                excitonCoefficients, 
                                                holeIndex);

            // Write hole position
            fprintf(textfile_rswf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            holePosition(0), holePosition(1), 0.0);

            double coefSum = 0;
            vec coefs = arma::zeros(bulkExciton.kpoints.n_elem*bulkExciton.motif.n_rows);
            int it = 0;


            for(int cell = -Ncell/2; cell <= Ncell/2; cell++){
    
                for (int n = 0; n < bulkExciton.motif.n_rows; n++){

                    coefs(it) = 
                    atomCoefficientSquared(n, cell, holeCell, RScoefs, bulkExciton);
                                
                    coefSum += coefs(it);
                    it++;
                };
            };

            it = 0;
            cout << "Writing w.f. to file" << endl;
            std::complex<double> imag(0,1);
            for(int cell = -Ncell/2; cell <= Ncell/2; cell++){
                    
                for (int n = 0; n < bulkExciton.motif.n_rows; n++){

                    rowvec position = bulkExciton.motif.row(n) + cell*arma::rowvec{0., ribbon.a, 0.};

                    fprintf(textfile_rswf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            position(0), position(1), coefs(it)/coefSum);
                    it++;
                };
            };
            fprintf(textfile_rswf, "#\n");
        }
    };

    fclose(textfile_en);
    fclose(textfile_rswf);
    fclose(textfile_kwf);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

	return 0;
};