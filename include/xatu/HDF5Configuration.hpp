#pragma once
#include <armadillo>
#include "xatu/SystemConfiguration.hpp"

namespace xatu {

class HDF5Configuration : public SystemConfiguration {

    private:
        std::string filename;
    
    public:
        HDF5Configuration(std::string);
        ~HDF5Configuration(){};

    private:
        void parseContent();
};

}