// Code to extract and calculate the real-space wavefunctions for several excitonic states from .txt file containing their expansion coefficients

#include <fstream>
#include <vector>
#include <armadillo>
#include <complex>
#include <string>

#include "../lib/libwavefunction.hpp"


int main(){

    std::string filename = "rs_wfs_0_to_5";
    FILE* file = fopen(filename.c_str(), "w");

    std::vector<std::complex<double>> vec1, vec2, vec3, vec4, vec5;

    double col1, col2, col3, col4, col5, col6, col7, col8, col9, col10;
    std::ifstream is("eigenstates_coefficients");
    while(!is.eof()){
        is >> col1 >> col2 >> col3 >> col4 >> col5 >> col6 >> col7 >> col8 >> col9 >> col10;
        std::complex<double> cpx1(col1, col2), cpx2(col3, col4), cpx3(col5, col6), cpx4(col7, col8), cpx5(col9, col10);
        vec1.push_back(cpx1);
        vec2.push_back(cpx2);
        vec3.push_back(cpx3);
        vec4.push_back(cpx4);
        vec5.push_back(cpx5);
    };
};



