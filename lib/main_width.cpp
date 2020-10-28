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
	std::string filename = "ribbon_width_maxE";
	FILE* textfile = fopen(filename.c_str(), "w");

    vec Narray = arma::regspace(5, 50);

    for(int n = 0; n < (int)Narray.n_elem; n++){

        int N = Narray(n);
        Zigzag zigzag = Zigzag(N);
        vec kpoints = arma::linspace(0, 2*PI/zigzag.a, 400);
        double maxE = -100;

        for(int n = 0; n < kpoints.n_elem; n++){

            double k = kpoints(n);
            arma::cx_mat h = zigzag.hamiltonian(k);

            vec eigval;
            arma::cx_mat eigvec;
            arma::eig_sym(eigval, eigvec, h);
            double EDiff = eigval(2*(N+1)*5) - eigval(2*(N+1)*5 - 2);
            if(EDiff > maxE){
                maxE = EDiff;
            };
        };
        
        fprintf(textfile, "%d\t%lf\n", N, maxE);
    }
    fclose(textfile);

    return 0;
};