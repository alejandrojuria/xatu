#pragma once
#include "xatu/ConfigurationSystem.hpp"

namespace xatu {

/**
 * The ConfigurationTB is used to parse the configuration from a custom file accepted by XATU.
 * Exclusive to the TB mode
 */
class ConfigurationTB : public ConfigurationSystem {

    protected:
        ConfigurationTB() = default;
    public:
        ConfigurationTB(const std::string&);
        
    protected:
        // Central method to parse the content from the configuration file
        void parseContent() override;
        arma::cx_cube parseMatrices(std::vector<std::string>&);
        arma::mat parseMotif(std::vector<std::string>&);
        arma::urowvec parseOrbitals(std::vector<std::string>&);
        // Method to check whether the parsed contents of the configuration file satisfy some basic requirements
        void checkContentCoherence();

};

}