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

        // Output and plotting
        void writeReciprocalAmplitude(int, FILE*);
        void writeExtendedReciprocalAmplitude(int, FILE*);
        void writeRealspaceAmplitude(int, int, const arma::rowvec&, FILE*);
        void writeEigenvalues(FILE*);

        double fourierTransformExciton(int, const arma::rowvec&, const arma::rowvec&);

        double realSpaceWavefunction(GExciton&, const arma::cx_vec&, int, int, 
                                     const arma::rowvec& eCell, const arma::rowvec& hCell);
        std::complex<double> densityMatrixElement(GExciton&, int, int, int, int);

        std::complex<double> densityMatrix(GExciton&, const arma::cx_vec&, int, int);
        std::complex<double> densityMatrixK(int, GExciton&, const arma::cx_vec&, int, int);

    private:
        int findExcitonPeak(int);
        arma::cx_mat RScoefficientCalc(int, int);
        double atomCoefficientSquared(int, const arma::rowvec&, const arma::rowvec&, 
                              const arma::cx_mat&);
        double boundingBoxBZ();
};