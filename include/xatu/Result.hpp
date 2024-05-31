#ifndef RESULT_HPP
#define RESULT_HPP

#pragma once
#include "xatu/Exciton.hpp"

namespace xatu {

template <typename T>
class Result {

    protected:
        // Exciton object which has been diagonalized.
        Exciton<T>* exciton_;
        // Pointer to System object used within Exciton
        std::shared_ptr<T> system_;
        // Exciton eigenenergies
        arma::vec m_eigval;
        // Exciton eigenstates
        arma::cx_mat m_eigvec;

    public:
        const Exciton<T>* exciton = exciton_;
        const std::shared_ptr<T>& system = system_;
        // Returns eigenvalues of BSE
        const arma::vec& eigval = m_eigval;
        // Returns eigenstates of BSE
        const arma::cx_mat& eigvec = m_eigvec;

    public:
        // Constructor
        Result(Exciton<T>* exciton_, arma::vec& eigval_, arma::cx_mat& eigvec_); 
        virtual ~Result(){};

        // Observables
        double kineticEnergy(int);
        double potentialEnergy(int);
        double bindingEnergy(int, double gap = -1);
        double determineGap(); 
        double ftExcitonEnvelope(int, const arma::rowvec&, const arma::rowvec&);
        arma::cx_vec spinX(int);
        virtual arma::cx_vec spinX(const arma::cx_vec&) = 0;
        virtual double realSpaceWavefunction(const arma::cx_vec&, int, int,
                                     const arma::colvec&, const arma::colvec&) = 0;

        // Output and plotting
        void writeReciprocalAmplitude(const arma::cx_vec&, FILE*);
        void writeReciprocalAmplitude(int, FILE*);
        void writeExtendedReciprocalAmplitude(const arma::cx_vec&, FILE*);
        void writeExtendedReciprocalAmplitude(int, FILE*);
        void writePhase(const arma::cx_vec&, FILE*);
        void writePhase(int, FILE*);
        void writeExtendedPhase(const arma::cx_vec&, FILE*);
        void writeExtendedPhase(int, FILE*);
        void writeEigenvalues(FILE*, int64_t n = 0);
        void writeStates(FILE*, int64_t n = 0);
        void writeSpin(int64_t, FILE*);
        void writeRealspaceAmplitude(int, int, const arma::colvec&, FILE*, int);
        virtual void writeRealspaceAmplitude(const arma::cx_vec&, int, const arma::colvec&, FILE*, int) = 0;
        virtual void writeAbsorptionSpectrum() = 0;

    protected:
        int findExcitonPeak(int);
        double boundingBoxBZ();
};

}

#include "xatu/Result.tpp"

#endif