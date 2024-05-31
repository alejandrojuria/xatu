#ifndef RESULTTB_HPP
#define RESULTTB_HPP

#pragma once
#include "xatu/SystemTB.hpp"
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

class ExcitonTB;

class ResultTB : public Result<SystemTB> {

    public:
        // Imaginary unit i
        static constexpr std::complex<double> imag {0., 1.};

    public:
        //// Constructor
        ResultTB(ExcitonTB* exciton_, arma::vec& eigval_, arma::cx_mat& eigvec_);
        ~ResultTB() = default;
        
        //// Observables
        // First recover hidden methods from Result<SystemTB>
        using Result<SystemTB>::spinX; // Add overload for spinX
        using Result<SystemTB>::writeRealspaceAmplitude;

        // Define additional methods
        arma::cx_vec spinX(const arma::cx_vec&) override;
        arma::mat velocity(int);
        arma::cx_vec velocitySingleParticle(int, int, int, std::string);
        arma::cx_mat excitonOscillatorStrength();
        double realSpaceWavefunction(const arma::cx_vec&, int, int,
                                     const arma::colvec&, const arma::colvec&) override;

        // Output and plotting
        void writeRealspaceAmplitude(const arma::cx_vec&, int, const arma::colvec&, FILE*, int ncells = 3) override;
        void writeAbsorptionSpectrum() override;        
        
    private:
        arma::cx_vec addExponential(arma::cx_vec&, const arma::colvec&);
};


}

#endif