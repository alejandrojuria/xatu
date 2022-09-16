#pragma once
#include <armadillo>
#include "SystemConfiguration.hpp"

class CrystalDFTConfiguration : public SystemConfiguration {

    public:
        int natoms;

    public:
        CrystalDFTConfiguration(std::string);

        void parseContent();

    private:
        void parseBravaisLattice();
        void parseAtoms();
};

