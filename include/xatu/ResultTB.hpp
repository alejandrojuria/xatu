#ifndef RESULTTB_HPP
#define RESULTTB_HPP

#pragma once
#include <armadillo>
#include "xatu/SystemTB.hpp"
#include "xatu/ExcitonTB.hpp"
#include "xatu/Result.hpp"

extern "C" {
    void skubo_w_(int* nR, int* norb, int* norb_ex, int* nv, int* nc, int* filling, 
                  double* Rvec, double* bravaisLattice, double* motif, 
                  std::complex<double>* hhop, double* shop, int* nk, double* rkx, 
                  double* rky, double* rkz, std::complex<double>* fk_ex, double* e_ex, 
                  double* eigval_stack, std::complex<double>* eigvec_stack);

    void exciton_oscillator_strength_(int* nR, int* norb, int* norb_ex, int* nv, int* nc, int* filling, 
                  double* Rvec, double* bravaisLattice, double* motif, 
                  std::complex<double>* hhop, double* shop, int* nk, double* rkx, 
                  double* rky, double* rkz, std::complex<double>* fk_ex, double* e_ex, 
                  double* eigval_stack, std::complex<double>* eigvec_stack, std::complex<double>* vme,
                  std::complex<double>* vme_ex, bool* convert_to_au);
}

namespace xatu {

class ResultTB : public Result<SystemTB> {

    public:
        // Constructor
        ResultTB(ExcitonTB* exciton_, arma::vec& eigval_, arma::cx_mat& eigvec_);
        
        // Observables
        using Result<SystemTB>::spinX; // Add overload for spinX
        arma::cx_vec spinX(const arma::cx_vec&);
        arma::mat velocity(int);
        arma::cx_vec velocitySingleParticle(int, int, int, std::string);
        arma::cx_mat excitonOscillatorStrength();
        double realSpaceWavefunction(const arma::cx_vec&, int, int,
                                     const arma::rowvec&, const arma::rowvec&);

        // Output and plotting
        void writeRealspaceAmplitude(const arma::cx_vec&, int, const arma::rowvec&, FILE*, int ncells = 3);
        void writeAbsorptionSpectrum();        
        
    private:
        arma::cx_vec addExponential(arma::cx_vec&, const arma::rowvec&);
};


}

#endif