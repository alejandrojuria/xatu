#pragma once
#include "xatu/SystemConfiguration.hpp"


typedef std::vector<std::vector<std::vector<double>>> cube_vector;

namespace xatu {

class CRYSTALConfiguration : public SystemConfiguration {

    public:
        int ndim, natoms, nsh, nspecies;
        int64_t norbitals;
        int totalElectrons, coreElectrons;
        bool SOC_FLAG = false;
        bool MAGNETIC_FLAG = false;
        bool alpha_electrons = true;
        arma::mat motif;
        // Basis of Bravais vectors (R1,R2,R3), stored by columns: (ndim x 3)
        arma::mat bravaisLattice;
        // List of Bravais vectors in Angstrom, stored by rows: (ncells x 3)
        arma::mat bravaisVectors;
        
        std::vector<int> atomic_number_ordering;
        std::vector<int> shellsPerSpecies;
        std::vector<int64_t> orbitalsPerSpecies;
        // std::map<int, std::vector<std::string>> shellTypesPerSpecies;
        std::map<int, cube_vector> gaussianCoefficients;

        arma::cx_cube overlapMatrices, fockMatrices;
        arma::cx_cube alphaMatrices, betaMatrices;

    public:
        CRYSTALConfiguration(std::string, int ncells);
        ~CRYSTALConfiguration(){};


        void parseContent(int);

    private:
        void parseBravaisLattice();
        void parseAtoms();
        void parseAtomicBasis();
        arma::cx_mat parseMatrix();
        void extractDimension();
    
    protected:
        void mapContent(bool debug = false);
        
};

size_t split(const std::string&, std::vector<double>&);

}
