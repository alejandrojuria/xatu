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
#include "libwavefunction.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]){

    cout << argv[0] << endl;
    if (argc != 2){
		cout << "Error: One argument must be given" << endl;
		exit(1);
	}

    auto start = high_resolution_clock::now();

    // ----------------- Model parameters & Output --------------------
    
    int nbands = 1;
    int nrmbands = 0;
    int Ncell = 20;
    double filling = 5./8;
    arma::rowvec Q = {0., 0., 0.};
    std::string modelfile = argv[1];
    bool useKcenteredMesh = false;

    std::string filename_en = "e_conv_approx_K_30f10";
    FILE* textfile_en = fopen(filename_en.c_str(), "w");
    bool writeEnergy = false;

    std::string filename_rswf = "hbn_rswf_approx_N30_hMid";
    FILE* textfile_rswf = fopen(filename_rswf.c_str(), "w");
    bool writeRSWF = true;

    std::string filename_bloch = "bloch_state_K";
    FILE* textfile_bloch = fopen(filename_bloch.c_str(), "w");

    // ----------------------- Main body -----------------------

    fprintf(textfile_en, "#k points\tTotal energy\tTB energy\tV energy (eV)\n");
    cx_mat eigvecX;
    vec eigvalX;

    cout << "##############################################################################" << endl;
    cout << "--------------------------------- Parameters ---------------------------------" << endl;
    cout << "N. cells: " << Ncell << endl;
    cout << "#bands: " << nbands << endl;
    cout << "#removed bands: " << nrmbands << "\n" << endl;

    GExciton bulkExciton = GExciton(modelfile, Ncell, Q, nbands, nrmbands, filling, true);
    cout << __LINE__ << endl;
    int nk = bulkExciton.nk;

    //double norm = arma::norm(bulkExciton.reciprocal_lattice.row(0));
    /*arma::rowvec K = norm/sqrt(3)*(bulkExciton.reciprocal_lattice.row(0)/2. -
                                bulkExciton.reciprocal_lattice.row(1)/2.)/arma::norm(
                                bulkExciton.reciprocal_lattice.row(0)/2. -
                                bulkExciton.reciprocal_lattice.row(1)/2.
                                );
    */   
   arma::rowvec K;
    if(useKcenteredMesh){
        
        int reductionFactor = 2;
        for(int i = 0; i < nk; i++){
            bulkExciton.kpoints.row(i) /= reductionFactor;
            bulkExciton.kpoints.row(i) += K;
        };
        bulkExciton.Ncell *= reductionFactor;
        bulkExciton.reinitializeInternals();
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
    cout << "Fifth state energy is: " << eigvalX(4) << endl;
    cout << "Delta E: " << eigvalX(1) - eigvalX(0) << endl;
    arma::vec energies = bulkExciton.computeEnergies(eigvecX.col(0));
    cout << "Kinectic energy: "  << energies(0) << endl;
    cout << "Potential energy: " << energies(1) << "\n" << endl;

    cx_vec groundStateExciton = eigvecX.col(0);
    cx_vec degenerateGSExciton = eigvecX.col(1);
    cx_vec totalExciton = groundStateExciton;

    // Print ground state kinetic, potential and total energies
    if(writeEnergy){
        fprintf(textfile_en, "#bands\t#k points\tTotal energy\tTB energy\tV energy (eV)\n");
        fprintf(textfile_en, 
                "%d\t%10.7lf\t%10.7lf\t%10.7lf\n", 
                nk, eigvalX(0), energies(0), energies(1));
    };

    // Print real-space wavefunction
    if(writeRSWF){
        for(int i = 0; i < 1; i++){

            //cx_vec excitonCoefficients = (eigvecX.col(0) + eigvecX.col(1) + eigvecX.col(2) + eigvecX.col(3))/2;
            //cx_vec excitonCoefficientsDeg = (eigvecX.col(4) + eigvecX.col(5) + eigvecX.col(6) + eigvecX.col(7))/2;
            //cx_vec exCoefs = (excitonCoefficients + excitonCoefficientsDeg)/sqrt(2);
            cx_vec excitonCoefficients = eigvecX.col(0);
            cx_vec excitonCoefficientsDeg = eigvecX.col(1);
            cx_vec excitonCoefficientsDeg2 = eigvecX.col(3);

            int holeIndex = bulkExciton.motif.n_rows/2;
            //arma::rowvec holeCell = Ncell/2*bulkExciton.bravais_lattice.row(0) + 
            //                        Ncell/2*bulkExciton.bravais_lattice.row(1); 
            //arma::rowvec holeCell = (Ncell/2)*bulkExciton.bravais_lattice.row(0);
            arma::rowvec holeCell = arma::rowvec{0., 0., 0.};
            rowvec holePosition = bulkExciton.motif.row(holeIndex) + holeCell;
            //excitonCoefficients.subvec(0, 5) = excitonCoefficients(0)*arma::ones<vec>(6);

            cx_mat RScoefs = RScoefficientCalc(bulkExciton, 
                                                excitonCoefficients, 
                                                holeIndex);
            /*cx_mat RScoefsDeg = RScoefficientCalc(bulkExciton, 
                                                excitonCoefficientsDeg, 
                                                holeIndex);*/

            // Write hole position
            fprintf(textfile_rswf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            holePosition(0), holePosition(1), 0.0);

            double coefSum = 0;
            
            int it = 0;

            //arma::mat cellCombinations = bulkExciton.generate_combinations_gamma(bulkExciton.Ncell, bulkExciton.ndim);
            double radius = arma::norm(bulkExciton.bravais_lattice.row(0)) * bulkExciton.Ncell / 2.5;
            arma::mat cellCombinations = bulkExciton.truncate_supercell(bulkExciton.Ncell, radius);
            vec coefs = arma::zeros(cellCombinations.n_rows*bulkExciton.motif.n_rows);
            cx_vec bloch_state;
            cx_mat bloch_states;
            vec energies;
            cout << __LINE__ << endl;
            //cx_mat h = bulkExciton.hamiltonian(K);
            //arma::eig_sym(energies, bloch_states, h);
            //bloch_state = bloch_states.col(1);

            for(int cell_index = 0; cell_index < cellCombinations.n_rows; cell_index++){
                
                arma::rowvec cell = arma::zeros<arma::rowvec>(3);
                for(int j = 0; j < bulkExciton.ndim; j++){
                    cell += cellCombinations.row(cell_index)(j)*bulkExciton.bravais_lattice.row(j);
                }
                cout << __LINE__ << endl;
                cell = cellCombinations.row(cell_index);
                for (int n = 0; n < bulkExciton.motif.n_rows; n++){
                    cout << __LINE__ << endl;

                    arma::rowvec electronPosition = cell + bulkExciton.motif.row(n);

                    if(false){
                    for (int alpha = 0; alpha < bulkExciton.norbitals; alpha++){
                        for (int beta = 0; beta < bulkExciton.norbitals; beta++){
                            int electronIndex = n*bulkExciton.norbitals + alpha;
                            int holeIndexWOrb = holeIndex*bulkExciton.norbitals + beta;
                            coefs(it) += 
                            realSpaceWavefunction(bulkExciton, excitonCoefficients, electronIndex, holeIndexWOrb, cell, holeCell);
                        }
                    }
                    }

                    //fourierTransformExciton(excitonCoefficients, bulkExciton, electronPosition, holePosition);
                    coefs(it) = atomCoefficientSquared(n, cell, holeCell, RScoefs, bulkExciton);
                    //atomCoefficientSquared(n, cell, holeCell, RScoefsDeg, bulkExciton);
                            
                    coefSum += coefs(it);
                    it++;
                };
            };

            it = 0;
            cout << "Writing w.f. to file" << endl;
            std::complex<double> imag(0,1);
            for(int cell_index = 0; cell_index < cellCombinations.n_rows; cell_index++){
                
                arma::rowvec cell = arma::zeros<arma::rowvec>(3);
                for(int j = 0; j < bulkExciton.ndim; j++){
                    cell += cellCombinations.row(cell_index)(j)*bulkExciton.bravais_lattice.row(j);
                }
                cell = cellCombinations.row(cell_index);
                for (int n = 0; n < bulkExciton.motif.n_rows; n++){

                    rowvec position = bulkExciton.motif.row(n) + cell;

                    //std::complex<double> tb_coef = bloch_state(n)*exp(imag*arma::dot(K, cell));

                    fprintf(textfile_rswf, "%11.8lf\t%11.8lf\t%14.11lf\n",
                            position(0), position(1), coefs(it)/coefSum);
                    //fprintf(textfile_bloch, "%11.8lf\t%11.8lf\t%14.11lf\n",
                    //        position(0), position(1), abs(tb_coef));
                    it++;
                };
            };
            fprintf(textfile_rswf, "#\n");
        }
    };

    fclose(textfile_en);
    fclose(textfile_rswf);
    fclose(textfile_bloch);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
	std::cout << "Elapsed time: " << duration.count()/1000.0 << " s" << std::endl;

	return 0;
};