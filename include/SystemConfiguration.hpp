#include <armadillo>
#include "ConfigurationBase.hpp"

class SystemConfiguration : public ConfigurationBase {

    struct configuration {
            int ndim;
            arma::mat bravaisLattice;
            arma::mat motif;
            arma::cx_cube hamiltonian;
            arma::cx_cube overlap;
            arma::mat bravaisVectors;
            arma::urowvec norbitals;
            std::map<std::string, int> atomToIndex;
            };
        
    public:
        configuration systemInfo;

    public:
        SystemConfiguration();
        SystemConfiguration(std::string);
        
    private:
        arma::mat parseVectors(std::vector<std::string>&);
        arma::cx_cube parseMatrices(std::vector<std::string>&);
        arma::mat parseMotif(std::vector<std::string>&);
        arma::urowvec parseOrbitals(std::vector<std::string>&);
        virtual void parseContent();
        void checkContentCoherence();

};