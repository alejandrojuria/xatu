#pragma once
#include <armadillo>
#include "GExciton.hpp"

class Result{

    private:
        GExciton& exciton;
        arma::vec eigval;
        arma::cx_mat eigvec;

    public:
        Result(GExciton&, arma::vec&, arma::cx_mat&);
        ~Result();

        // Observables
        double kineticEnergy(int);
        double potentialEnergy(int);
        double bindingEnergy(int, double gap = -1);
        double determineGap();
        arma::cx_vec spinX(int);

        // Output and plotting
        void writeReciprocalAmplitude(int, FILE*);
        void writeExtendedReciprocalAmplitude(int, FILE*);
        void writeRealspaceAmplitude(int, FILE*);

        arma::cx_mat RScoefficientCalc(int, int);
        double atomCoefficientSquared(int, const arma::rowvec&, const arma::rowvec&, 
                              const arma::cx_mat&);
        double fourierTransformExciton(int, const arma::rowvec&, const arma::rowvec&);

        double realSpaceWavefunction(GExciton&, const arma::cx_vec&, int, int, const arma::rowvec& eCell, const arma::rowvec& hCell);
        std::complex<double> densityMatrixElement(GExciton&, int, int, int, int);

        std::complex<double> densityMatrix(GExciton&, const arma::cx_vec&, int, int);
        std::complex<double> densityMatrixK(int, GExciton&, const arma::cx_vec&, int, int);

    private:
        int findExcitonPeak(int);
};