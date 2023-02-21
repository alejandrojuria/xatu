#pragma once
#include "ConfigurationBase.hpp"
#include <armadillo>

namespace xatu {

class ExcitonConfiguration : public ConfigurationBase{

    struct configuration {
        int ncell, nbands;
        int submeshFactor = 1;
        arma::ivec bands = {};
        arma::rowvec Q = {0., 0., 0.};
        arma::rowvec shift;
        double cutoff;
        arma::vec eps = {};
        double r0, d;
        std::string mode = "realspace";
        bool exchange = false;
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

}