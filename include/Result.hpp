#pragma once
#include <armadillo>
#include "GExciton.hpp"

class Result{

    private:
        GExciton& exciton;
        arma::mat eigval;
        arma::cx_cube eigvec;

    public:
        Result();
        Result(GExciton&, arma::mat&, arma::cx_cube&);
        ~Result();

        // Observables
        double computeKineticEnergy(int);
        double computePotentialEnergy(int);
        double computeBindingEnergy(int);
        double determineGap();
        arma::cx_vec spinX(int);

        // Output and plotting
        void writeReciprocalAmplitude(int, FILE*);
        void writeExtendedReciprocalAmplitude(int, FILE*);
        void writeRealspaceAmplitude(int, FILE*);
};