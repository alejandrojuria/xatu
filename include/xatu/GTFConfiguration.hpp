#pragma once
#include "xatu/CRYSTALConfiguration.hpp"

namespace xatu {

class GTFConfiguration : public CRYSTALConfiguration {

    public:
        //Vector assigning the vector of shells (itself a vector of contracted GTF) to each atomic species.
        //One for the DFT basis and one for the auxiliary basis
        cube_vector Shells_all_species_SCF, Shells_all_species_AUX;
        //Vector assigning the vector of L (ang. mom. quant. num.) to each atomic species
        std::vector<std::vector<int>> L_all_species_SCF, L_all_species_AUX;
        std::string b_filename;

    protected:
        std::ifstream b_file;

    public:
        GTFConfiguration(std::string, std::string, int ncells = 50);
        virtual ~GTFConfiguration(){};

    protected:  
        double replace_D_2_double(std::string);

    private:
        void skip_PP();
        void parseBasis(bool);
        void parseBases();
        
};

}