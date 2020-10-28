#include <iostream>
#include <fstream>
#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>

#include "Zigzag.hpp"
#include "Exciton.hpp"

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

using namespace arma;
using namespace std::chrono;

int main(){
	std::string filename = "dos_eh_pair";
	FILE* textfile = fopen(filename.c_str(), "w");

    int N = 15;
    int Ncell = 150;
    double Q = 0.0;
    int bulkBands = 0;
    int edgeBands = 2;
    double delta = 0.01;
    Exciton exciton = Exciton(N, Ncell, Q, bulkBands, edgeBands);

    vec energyArray = arma::linspace(-1, 1, 200);
    for(int i = 0; i < (int)energyArray.n_elem; i++){
        double energy = energyArray(i);
        double dos = exciton.pairDensityOfStates(energy, delta);

        fprintf(textfile, "%lf\t%lf\n", energy, dos);
    };
    fclose(textfile);

    return 0;
};