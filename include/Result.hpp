#pragma once
#include <armadillo>
#include "GExciton.hpp"
#include "forward_declaration.hpp"

class Result{

    private:
        GExciton& exciton;
        arma::vec m_eigval;
        arma::cx_mat m_eigvec;

    public:
        // Returns eigenvalues of BSE
        const arma::vec& eigval = m_eigval;
        // Returns eigenstates of BSE
        const arma::cx_mat& eigvec = m_eigvec;

    public:
        Result(GExciton& exciton_, arma::vec& eigval_, arma::cx_mat& eigvec_) : 
                exciton(exciton_),
                m_eigval(eigval_),
                m_eigvec(eigvec_){};

        // Observables
        double kineticEnergy(int);
        double potentialEnergy(int);
        double bindingEnergy(int, double gap = -1);
        double determineGap();
        arma::cx_vec spinX(int);

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

        double fourierTransformExciton(int, const arma::rowvec&, const arma::rowvec&);

        double realSpaceWavefunction(const arma::cx_vec&, int, int,
                                     const arma::rowvec&, const arma::rowvec&);
        
        std::complex<double> densityMatrixElement(GExciton&, int, int, int, int);

        std::complex<double> densityMatrix(GExciton&, const arma::cx_vec&, int, int);
        std::complex<double> densityMatrixK(int, GExciton&, const arma::cx_vec&, int, int);

    private:
        int findExcitonPeak(int);
        double boundingBoxBZ();
        arma::cx_vec addExponential(arma::cx_vec&, const arma::rowvec&);
};