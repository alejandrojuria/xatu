#define CATCH_CONFIG_MAIN

#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <string>
#include <vector>

#include <xatu.hpp>
#include <catch.hpp>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


TEST_CASE("Modelfile parsing", "[model_parsing]"){
    
    std::string modelfile = "../models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    int expectedFilling = 1;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {1, 1};

    arma::vec expectedBravaisLatticeHash = {5.88755, 3.38755};
    arma::vec expectedMotifHash = {0, 5.40211};
    arma::vec expectedBravaisVectorsHash = {0, 1.22249, 3.72249, 3.38755, 5.88755};
    arma::vec expectedHamiltionianHash = {1.65625, -0.15, -0.15, 0.7125, 0.7125};
    
    REQUIRE(config.systemInfo.ndim == expectedDim);
    REQUIRE(config.systemInfo.filling == expectedFilling);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        REQUIRE(config.systemInfo.norbitals(i) == expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisLatticeHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedMotifHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisVectorsHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedHamiltionianHash(i), 1E-4));
    }        
}

TEST_CASE("CRYSTAL file parsing", "[CRYSTAL_parsing]"){

    std::string modelfile = "../models/DFT/hBN_base_HSE06.outp";
    xatu::CrystalDFTConfiguration config = xatu::CrystalDFTConfiguration(modelfile);

    int expectedFilling = 6;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {33, 36};
    uint expectedFockMatricesNumber = 20;

    arma::vec expectedBravaisLatticeHash = {3.3931, 5.18333};
    arma::vec expectedMotifHash = {0.818875, 5.09563};
    arma::vec expectedBravaisVectorsHash = {0, 1.21938, 5.9031, 3.72938, 3.3931, 0.163333, 
                                            5.18333, 0.275427, 4.62287, 0.382713, 10.0864,
                                            7.91271, 2.55643, 0.43876, 9.8062, 5.45876,
                                            4.7862, -0.673333, 9.366674, 0.494807};
    arma::vec expectedHamiltionianHash = {2468.69, 2108.44, 2104.88, 1793.87,
                                            2254.54, 2258.97, 1792.2, 1009.11, 1325.13,
                                            1325.32, 1009.1, 928.996, 1636.76, 856.064,
                                            856.041, 734.011, 1172, 1172.08, 734.016, 404.001};
    arma::vec expectedOverlapHash = {897.123, 952.327, 1197.62, 1042.23, 1303.16,
                                        1246.87, 980.796, 445.932, 623.911, 767.765,
                                        589.939, 543.999, 924.147, 508.966, 508.98,
                                        437.995, 687.007, 660.933, 413.993, 243};
    
    REQUIRE(config.systemInfo.ndim == expectedDim);
    REQUIRE(config.systemInfo.filling == expectedFilling);
    REQUIRE(config.systemInfo.hamiltonian.n_slices == expectedFockMatricesNumber);
    REQUIRE(config.systemInfo.overlap.n_slices == expectedFockMatricesNumber);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        REQUIRE(config.systemInfo.norbitals(i) == expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisLatticeHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedMotifHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedBravaisVectorsHash(i), 1E-4));
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedHamiltionianHash(i), 1E-2));
    } 
    for(uint i = 0; i < config.systemInfo.overlap.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.overlap.slice(i));
        REQUIRE_THAT(hash, Catch::Matchers::WithinAbs(expectedOverlapHash(i), 1E-2));
    }
}

TEST_CASE("Exciton file parsing", "[excitonfile_parsing]"){

    std::string excitonconfig = "../excitonconfig/hBN_spinless.txt";
    xatu::ExcitonConfiguration config = xatu::ExcitonConfiguration(excitonconfig);

    std::string expectedName = "hBN_N30";
    int expectedNcell = 30;
    int expectedNbands = 1;
    arma::vec expectedDielectric = {1, 1, 10};

    REQUIRE(config.excitonInfo.ncell == expectedNcell);
    REQUIRE(config.excitonInfo.nbands == expectedNbands);
    for (uint i = 0; i < expectedDielectric.n_elem; i++){
        REQUIRE(config.excitonInfo.eps(i) == expectedDielectric(i));
    }
    REQUIRE(xatu::array2hash(config.excitonInfo.Q) == 0);
}
    

    // xatu::Exciton bulkExciton = xatu::Exciton(config, ncell, nbands, nrmbands, parameters);
    // bulkExciton.setMode("realspace");

    // bulkExciton.brillouinZoneMesh(ncell);



    // bulkExciton.initializeHamiltonian();
    // bulkExciton.BShamiltonian();
    // auto results = bulkExciton.diagonalize("diag", nstates);

    // auto energies = xatu::detectDegeneracies(results.eigval, nstates, 6);
    
    // std::vector<std::vector<double>> expectedEnergies = {{5.335687, 2}, {6.073800, 1}, {6.164057, 2}, {6.172253, 1}, {6.351066, 2}};
    // for(int i = 0; i < energies.size(); i++){
    //     REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-5));
    //     REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    // }

    
