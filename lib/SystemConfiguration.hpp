#include <armadillo>
#include "ConfigurationBase.hpp"

class SystemConfiguration : public ConfigurationBase {

    struct configuration {
            int ndim;
            arma::mat bravaisLattice;
            arma::mat motif;
            arma::cx_cube hamiltonian;
            arma::mat bravaisVectors;
            arma::rowvec norbitals;
            std::map<std::string, int> atomToIndex;
            };

    public:
        SystemConfiguration();
        SystemConfiguration(std::string);
        
    private:
        void mapFunctionsToArgs();
        arma::mat parseBravaisLattice(std::vector<std::string>&);
        arma::cx_cube parseHamiltonian(std::vector<std::string>&);
        arma::mat parseMotif(std::vector<std::string>&);
        arma::rowvec parseOrbitals(std::vector<std::string>&);
        virtual void parseContent();
        void checkContentCoherence();

    public:
        configuration systemInfo;


};