#pragma once
#include <armadillo>
#include "xatu/Exciton.hpp"
#include "xatu/forward_declaration.hpp"


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

class Result{

    private:
        // Exciton object which has been diagonalized.
        Exciton& exciton;
        // Exciton eigenenergies
        arma::vec m_eigval;
        // Exciton eigenstates
        arma::cx_mat m_eigvec;

    public:
        // Returns eigenvalues of BSE
        const arma::vec& eigval = m_eigval;
        // Returns eigenstates of BSE
        const arma::cx_mat& eigvec = m_eigvec;

    public:
        Result(Exciton& exciton_, arma::vec& eigval_, arma::cx_mat& eigvec_) : 
                exciton(exciton_),
                m_eigval(eigval_),
                m_eigvec(eigvec_){};

        // Observables
        double kineticEnergy(int);
        double potentialEnergy(int);
        double bindingEnergy(int, double gap = -1);
        double determineGap();
        arma::cx_vec spinX(int);
        arma::mat velocity(int);
        arma::cx_vec velocitySingleParticle(int, int, int, std::string);
        arma::cx_mat excitonOscillatorStrength();

        // Symmetries
        arma::cx_mat diagonalizeC3(const arma::vec&);
        arma::cx_mat symmetrizeStates(const arma::cx_vec&, const arma::cx_vec&);

        // Output and plotting
        void writeReciprocalAmplitude(const arma::cx_vec&, FILE*);
        void writeReciprocalAmplitude(int, FILE*);
        void writePhase(const arma::cx_vec&, FILE*);
        void writePhase(int, FILE*);
        void writeExtendedReciprocalAmplitude(const arma::cx_vec&, FILE*);
        void writeExtendedReciprocalAmplitude(int, FILE*);
        void writeExtendedPhase(const arma::cx_vec&, FILE*);
        void writeExtendedPhase(int, FILE*);
        void writeRealspaceAmplitude(const arma::cx_vec&, int, const arma::rowvec&, FILE*, int ncells = 3);
        void writeRealspaceAmplitude(int, int, const arma::rowvec&, FILE*, int ncells = 3);
        void writeEigenvalues(FILE*, int n = 0);
        void writeStates(FILE*, int n = 0);
        void writeAbsorptionSpectrum();
        void writeSpin(int, FILE*);
        

        double fourierTransformExciton(int, const arma::rowvec&, const arma::rowvec&);

        double realSpaceWavefunction(const arma::cx_vec&, int, int,
                                     const arma::rowvec&, const arma::rowvec&);
        
        std::complex<double> densityMatrixElement(Exciton&, int, int, int, int);

        std::complex<double> densityMatrix(Exciton&, const arma::cx_vec&, int, int);
        std::complex<double> densityMatrixK(int, Exciton&, const arma::cx_vec&, int, int);

    private:
        int findExcitonPeak(int);
        double boundingBoxBZ();
        arma::cx_vec addExponential(arma::cx_vec&, const arma::rowvec&);
};

}