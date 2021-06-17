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

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){

    auto start = high_resolution_clock::now();
    
    //int N = 15;
    int nBulkBands = 2;
    int nEdgeBands = 0;
    int N = 5;
    //NcellArray = arma::join_cols(NcellArray, arma::regspace(80, 20, 200));
    Zigzag ribbon = Zigzag(N);
    vec QArray = {0.0};

    // ----------------- Model parameters & Output --------------------

    std::string filename = "exciton_spectrum_exact";
    FILE* textfile = fopen(filename.c_str(), "w");
    bool writeEigvals = true;

    std::string filename_wf = "exciton_bulk_wf_N15";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
    bool printWF = false;

    bool calculateSpin = false;
    std::string filename_spin = "exciton_bulk_spin_6b_z1em7";
    FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
    bool writeSpin = false;

    // ----------------------- Main body -----------------------

    fprintf(textfile, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");
    for(int n = 0; n < (int)QArray.n_elem; n++){
        int Ncell = 1;
        int nk = 2*Ncell + 1;
        double Q = QArray(n);
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
            fprintf(textfile, "%10.7lf\t", Q);
            for(int i = 0; i < eigvalX.n_elem; i++){
                fprintf(textfile, "%10.7lf\t", eigvalX(i));
            }
            fprintf(textfile, "\n");
        };

        //cx_mat states = bulkExciton.fixDegeneracy(eigvecX.col(0), eigvecX.col(1), 20);
        //cx_vec state = states.col(0);

        if(printWF){
            fprintf(textfile_wf, "%d\n", N);
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
                fprintf(textfile_wf, "%11.8lf\t%11.8lf\n", bulkExciton.kpoints(i), coef);
            };
            fprintf(textfile_wf, "#\n");
        };

        cx_vec spin;
        if(calculateSpin){
            spin = bulkExciton.spinX(eigvecX.col(0));
            cout << "#------------------------ Spin ------------------------" << endl;
            cout << "Hole Sz is: " << spin(0) << endl;
            cout << "Electron Sz is: " << spin(1) << endl;
            cout << "Total spin is: " << spin(2) << "\n" << endl;
                
            if(writeSpin){
                fprintf(textfile_spin, "%d\t%13.10lf\n", nk, real(spin(2)));
            };
        };
    };

    fclose(textfile);
    fclose(textfile_wf);
    fclose(textfile_spin);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};