#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "GExciton.hpp"
#include "System.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    if (argc > 3){
		cout << "Error: No more than three arguments can be given" << endl;
		exit(1);
	}

    auto start = high_resolution_clock::now();
    
    int nbands = 2;
    int nrmbands = 0;
    int Ncell = 21;
    double filling = 5./8;
    arma::vec NcellArray = {21};
    arma::rowvec Q = {0., 0., 0.};
    std::string modelfile = argv[1];
    try{
        std::string overlapfile = argv[2];
    }
    catch(const std::exception& e){}

    
    bool useCenteredMeshK = false;

    // ----------------- Model parameters & Output --------------------

    std::string filename = "e_conv_approx_K_30f10";
    FILE* textfile_en = fopen(filename.c_str(), "w");
    bool writeEigvals = false;

    std::string filename_wf = "bi_ribbon_kwf_approx_N40";
    FILE* textfile_wf = fopen(filename_wf.c_str(), "w");
    bool printWF = true;

    bool calculateSpin = false;
    std::string filename_spin = "gexciton_bulk_spin_6b_z1em7";
    FILE* textfile_spin = fopen(filename_spin.c_str(), "w");
    bool writeSpin = false;

    // ----------------------- Main body -----------------------

    fprintf(textfile_en, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");
    for(int n = 0; n < NcellArray.n_elem; n++){
        int Ncell = NcellArray(n);
        cx_mat eigvecX;
        vec eigvalX;

        cout << "##############################################################################" << endl;
        cout << "--------------------------------- Parameters ---------------------------------" << endl;
        cout << "N. cells: " << Ncell << endl;
        cout << "#bands: " << nbands << endl;
        cout << "#removed bands: " << nrmbands << "\n" << endl;

        GExciton bulkExciton = GExciton(modelfile, Ncell, Q, nbands, nrmbands, filling, true);
        int nk = bulkExciton.nk;
        Ncell = bulkExciton.Ncell;
        /*double norm = arma::norm(bulkExciton.bravais_lattice.row(0)) * Ncell/2.5;
        cout << norm << endl;
        arma::mat combinations = bulkExciton.generate_combinations_gamma(Ncell, 2);
        for (int i = 0; i < combinations.n_rows; i++){
            arma::rowvec lattice_vector = combinations.row(i)(0)*bulkExciton.bravais_lattice.row(0) + 
            combinations.row(i)(1)*bulkExciton.bravais_lattice.row(1);
            if (arma::norm(lattice_vector) < norm + 1E-5){
                fprintf(textfile_wf, "%10.7lf\t%10.7lf\t%d\n", 
                    lattice_vector(0), lattice_vector(1), 1);
            }
            else{
                fprintf(textfile_wf, "%10.7lf\t%10.7lf\t%d\n", 
                    lattice_vector(0), lattice_vector(1), 0);
            }
        }*/
        
        
        if(useCenteredMeshK){
            double norm = arma::norm(bulkExciton.reciprocal_lattice.row(0));
            arma::rowvec K = norm/sqrt(3)*(bulkExciton.reciprocal_lattice.row(0)/2. -
                                    bulkExciton.reciprocal_lattice.row(1)/2.)/arma::norm(
                                    bulkExciton.reciprocal_lattice.row(0)/2. -
                                    bulkExciton.reciprocal_lattice.row(1)/2.);
                                    
            int reductionFactor = 4;
            for(int i = 0; i < nk; i++){
                bulkExciton.kpoints.row(i) /= reductionFactor;
                bulkExciton.kpoints.row(i) += K;
            };
            bulkExciton.Ncell *= reductionFactor;
            bulkExciton.reinitializeInternals();
            cout << bulkExciton.Ncell << ", " << bulkExciton.nk << endl;
            bulkExciton.nk = pow(bulkExciton.Ncell, bulkExciton.ndim);
        }

        cout << "Valence bands:" << bulkExciton.valenceBands << endl;
        cout << "Conduction bands:" << bulkExciton.conductionBands << endl;
        bulkExciton.BShamiltonian();
        arma::eig_sym(eigvalX, eigvecX, bulkExciton.HBS);

        cout << "#----------------------- Energy -----------------------" << endl;
        cout << "Fundamental state energy is: " << eigvalX(0) << endl;
        cout << "Second state energy is: " << eigvalX(1) << endl;
        cout << "Third state energy is: " << eigvalX(2) << endl;
        cout << "Fourth state energy is: " << eigvalX(3) << endl;
        cout << "Delta E: " << eigvalX(1) - eigvalX(0) << endl;
        vec energies = bulkExciton.computeEnergies(eigvecX.col(0));
        cout << "Kinectic energy: "  << energies(0) << endl;
        cout << "Potential energy: " << energies(1) << "\n" << endl;

        if(writeEigvals){
            //writeEigenvaluesToFile(textfile, eigvalX, Q);
            fprintf(textfile_en, "%d\t", Ncell);
            for(int i = 0; i < eigvalX.n_elem; i++){
                fprintf(textfile_en, "%10.7lf\t", eigvalX(i));
            }
            fprintf(textfile_en, "\n");
        };

        //cx_mat states = bulkExciton.fixDegeneracy(eigvecX.col(0), eigvecX.col(1), 20);
        //cx_vec state = states.col(0);

        if(printWF){
            fprintf(textfile_wf, "kx\tky\tkz\tProb.\n");
            cx_vec state = eigvecX.col(0);
            //cx_vec state = (eigvecX.col(0) + eigvecX.col(1))/sqrt(2);
            cx_vec deg_state = eigvecX.col(0);
            int nbandsCombinations = bulkExciton.conductionBands.n_elem*
                                    bulkExciton.valenceBands.n_elem;
            for (int i = 0; i < bulkExciton.kpoints.n_rows; i++){
                double coef = 0;
                for(int nband = 0; nband < nbandsCombinations; nband++){
                    coef += abs(state(nbandsCombinations*i + nband))*
                            abs(state(nbandsCombinations*i + nband));
                    //coef += abs(eigvecX.col(state + 1)(nbandsCombinations*i + nband))*abs(eigvecX.col(state + 1)(nbands2*i + nband));
                };
                coef /= arma::norm(bulkExciton.kpoints.row(1) - bulkExciton.kpoints.row(0)); // L2 norm instead of l2
                fprintf(textfile_wf, "%11.8lf\t%11.8lf\t%11.8lf\t%11.8lf\n", 
                        bulkExciton.kpoints.row(i)(0), bulkExciton.kpoints.row(i)(1), 
                        bulkExciton.kpoints.row(i)(2), coef);
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
    }
    

    fclose(textfile_en);
    fclose(textfile_wf);
    fclose(textfile_spin);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

    return 0;
};