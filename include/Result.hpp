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

    private:
        int findExcitonPeak(int);
};