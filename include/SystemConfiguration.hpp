#pragma once
#include <armadillo>
#include "ConfigurationBase.hpp"

class SystemConfiguration : public virtual ConfigurationBase {

    struct configuration {
            int ndim;
            double filling;
            arma::mat bravaisLattice;
            arma::mat motif;
            arma::cx_cube hamiltonian;
            arma::cx_cube overlap;
            arma::mat bravaisVectors;
            arma::urowvec norbitals;
            };
        
    public:
        configuration systemInfo;

    protected:
        SystemConfiguration();
    public:
        SystemConfiguration(std::string);
        
    protected:
        arma::mat parseVectors(std::vector<std::string>&);
        arma::cx_cube parseMatrices(std::vector<std::string>&);
        arma::mat parseMotif(std::vector<std::string>&);
        arma::urowvec parseOrbitals(std::vector<std::string>&);
        virtual void parseContent();
        void checkContentCoherence();

};