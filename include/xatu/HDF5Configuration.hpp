#pragma once
#include "xatu/ConfigurationSystem.hpp"

namespace xatu {

class HDF5Configuration : public ConfigurationSystem {

    private:
        std::string filename;
    
    public:
        HDF5Configuration(const std::string&);
        ~HDF5Configuration(){};

    private:
        void parseContent();
};

}