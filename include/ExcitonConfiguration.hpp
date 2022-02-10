#include "ConfigurationBase.hpp"
#include <armadillo>

class ExcitonConfiguration : public ConfigurationBase{

    struct configuration {
        int Ncell, nbands, nrmbands, filling;
        arma::vec bands={};
        arma::rowvec Q;
        bool useApproximation;

    };

    public:
        configuration excitonInfo;
    
    public:
        ExcitonConfiguration();
        ExcitonConfiguration(std::string);
    
    private:
        virtual void parseContent();
        void checkContentCoherence();

};