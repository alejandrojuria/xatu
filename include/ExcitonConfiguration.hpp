#include "ConfigurationBase.hpp"
#include <armadillo>

class ExcitonConfiguration : public ConfigurationBase{

    struct configuration {
        int ncell, nbands, nrmbands, filling;
        arma::urowvec bands = {};
        arma::rowvec Q;
        bool useApproximation;
        double cutoff;
        arma::rowvec eps = {};
        double r0, d;
    };

    public:
        configuration excitonInfo;
    
    public:
        ExcitonConfiguration();
        ExcitonConfiguration(std::string);
    
    private:
        virtual void parseContent();
        void checkContentCoherence();
        double parseFraction(std::string&);
};