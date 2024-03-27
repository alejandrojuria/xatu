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

    std::cout << std::setw(40) << std::left << "Testing model file parsing... ";
    std::cout.setstate(std::ios_base::failbit);
    
    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    int expectedFilling = 1;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {1, 1};

    arma::vec expectedBravaisLatticeHash = {4.55422, 2.05422};
    arma::vec expectedMotifHash = {0, 3.90211};
    arma::vec expectedBravaisVectorsHash = {0, -0.110843, 2.38916, 2.05422, 4.55422};
    arma::vec expectedHamiltionianHash = {-1.34375, -0.9, -0.9, -0.0375, -0.0375};
    
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

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("CRYSTAL file parsing", "[CRYSTAL_parsing]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing CRYSTAL file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config = xatu::CRYSTALConfiguration(modelfile);

    int expectedFilling = 6;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {33, 36};
    uint expectedFockMatricesNumber = 20;

    arma::vec expectedBravaisLatticeHash = {2.05977, 4.51667};
    arma::vec expectedMotifHash = {0.068875, 2.84563};
    arma::vec expectedBravaisVectorsHash = {0, -0.113953, 4.56977, 2.39605, 2.05977,
                                            -0.503333, 4.51667, -0.39124, 3.9562, -0.95062,
                                            8.7531, 6.57938, 1.2231, -0.894573, 8.47287,
                                            4.12543, 3.45287, -1.34, 8.7, -0.838527};
    arma::vec expectedHamiltionianHash = {61.2007, 2.87756, -0.679479, 0.246344, -0.981652,
                                          0.447539, -1.41873, 0.32476, 0.408354, 0.596467,
                                          0.309267, 0.191303, 0.107262, 0.244018, 0.220732,
                                          0.165454, 0.248117, 0.32759, 0.169674, 0.0856007};
    arma::vec expectedOverlapHash = {46.302, -1.47274, -0.128825, 0.444089, 1.43755,
                                    -0.87285, 1.00188, 0.0256365, 0.0421825, -0.0733795,
                                    0.0628083, 0.113745, 0.340853, 0.0733232, 0.0869901,
                                    0.0867042, 0.151335, 0.0713888, 0.0795013, 0.0507182};
    
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

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("Exciton file parsing", "[excitonfile_parsing]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing exciton file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string excitonconfig = "./data/hBN_spinless.txt";
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

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}
    
TEST_CASE("TB hBN energies (full diagonalization)", "[tb-hBN-fulldiag]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (fulldiag)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (davidson)", "[tb-hBN-davidson]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (Davidson)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("davidson", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (Lanczos)", "[hBN-lanczos]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (Lanczos)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("sparse", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (reciprocal)", "[tb-hBN-reciprocal]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN energies (reciprocal)... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.setMode("reciprocalspace");
    exciton.setReciprocalVectors(5);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{6.234291, 1}, 
                                                         {6.236636, 1},
                                                         {6.731819, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN reciprocal w.f.", "[tb-hBN-kwf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN reciprocal w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.086145105;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN real-space w.f.", "[tb-hBN-rswf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN real-space w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;

    double radius = arma::norm(exciton.system->bravaisLattice.row(0)) * exciton.ncell;
    arma::mat cellCombinations = exciton.system->truncateSupercell(exciton.ncell, radius);
    arma::vec rswf = arma::zeros(cellCombinations.n_rows*exciton.system->motif.n_rows);

    // Compute probabilities
    for(int n = 0; n < nstates; n++){
        int it = 0;
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.system->motif.n_rows; atomIndex++){
            rswf(it) += results->realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
        }
    }

    double rswfHash = xatu::array2hash(rswf);
    double expectedRSwfHash = 3.4185866144;
    REQUIRE_THAT(rswfHash, Catch::Matchers::WithinAbs(expectedRSwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN absorption", "[tb-hBN-kubo]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN conductivity... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 47.4140063784;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("TB hBN energies (spinful)", "[tb-hBN-spinful]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN spinful... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 12;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 8}, 
                                                         {6.074062, 4}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}
    
TEST_CASE("TB hBN spin", "[tb-hBN-spin]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing TB hBN spin... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 16;
    int nstates = 4;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});
    
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::mat expectedSpin = arma::vec{-1, 0, 0, 1};  

    for(uint i = 0; i < nstates; i++){
        arma::cx_vec spin = results->spinX(i);
        REQUIRE(spin.n_elem == 3);
        REQUIRE_THAT(std::real(spin(0)), Catch::Matchers::WithinAbs(expectedSpin(i), 1E-2));
        REQUIRE_THAT(std::abs(spin(1)), Catch::Matchers::WithinAbs(0.5, 1E-2));
        REQUIRE_THAT(std::abs(spin(2)), Catch::Matchers::WithinAbs(0.5, 1E-2));
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("DFT hBN", "[dft-hBN]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing DFT hBN... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    bool triangular = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config = xatu::CRYSTALConfiguration(modelfile, 100);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.system->setAU(true);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{4.442317, 1}, 
                                                         {4.442427, 1}};
                                                         
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    // Check reciprocal w.f.
    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.0864895707;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    // Check realspace w.f.
    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;

    double radius = arma::norm(exciton.system->bravaisLattice.row(0)) * exciton.ncell;
    arma::mat cellCombinations = exciton.system->truncateSupercell(exciton.ncell, 2);
    arma::vec rswf = arma::zeros(cellCombinations.n_rows*exciton.system->motif.n_rows);

    // Compute probabilities
    for(int n = 0; n < nstates; n++){
        int it = 0;
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for(unsigned int cellIndex = 0; cellIndex < cellCombinations.n_rows; cellIndex++){
        arma::rowvec cell = cellCombinations.row(cellIndex);
        for (unsigned int atomIndex = 0; atomIndex < exciton.system->motif.n_rows; atomIndex++){
            rswf(it) += results->realSpaceWavefunction(statecoefs, atomIndex, holeIndex, cell, holeCell);
            it++;
        }
        }
    }

    double rswfHash = xatu::array2hash(rswf);
    double expectedRSwfHash = 83.2242560463;
    REQUIRE_THAT(rswfHash, Catch::Matchers::WithinAbs(expectedRSwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 energies", "[MoS2-energies]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 energies... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.768783, 2}, 
                                                         {1.780562, 2}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 reciprocal w.f.", "[MoS2-kwf]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 reciprocal w.f... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    int nbandsCombinations = exciton.conductionBands.n_elem * exciton.valenceBands.n_elem;
    arma::cx_vec kwf = arma::zeros<arma::cx_vec>(exciton.system->kpoints.n_rows);
    for (int n = 0; n < nstates; n++){
        arma::cx_vec statecoefs = results->eigvec.col(n);
        for (int i = 0; i < exciton.system->kpoints.n_rows; i++){
        double coef = 0;
        for(int nband = 0; nband < nbandsCombinations; nband++){
            coef += abs(statecoefs(nbandsCombinations*i + nband))*
                    abs(statecoefs(nbandsCombinations*i + nband));
        };
        coef /= arma::norm(exciton.system->kpoints.row(1) - exciton.system->kpoints.row(0)); // L2 norm instead of l2
        kwf(i) += coef;
        };
    }

    double kwfHash = xatu::array2hash(kwf);
    double expectedKwfHash = 1.1814790902;
    REQUIRE_THAT(kwfHash, Catch::Matchers::WithinAbs(expectedKwfHash, 1E-5));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 spin", "[MoS2-spin]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 spin... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.system->reducedBrillouinZoneMesh(ncell, factor);
    exciton.system->shiftBZ({0.6628, -1.1480, 0});

    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::mat expectedSpin = arma::mat{{-9.915e-01, -5.000e-01, -4.915e-01},
                                       {-6.966e-04, -5.000e-01,  4.993e-01},
                                       { 8.297e-03,  4.998e-01, -4.915e-01},
                                       { 9.991e-01,  4.998e-01,  4.993e-01}}; 
                              
    for(uint i = 0; i < nstates; i++){
        arma::cx_vec spin = results->spinX(i);
        for(uint j = 0; j < 3; j++){
            double spinValue = real(spin(j));
            REQUIRE_THAT(spinValue, Catch::Matchers::WithinAbs(expectedSpin(i, j), 1E-4));
        }
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 exciton bands", "[MoS2-Q]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 exciton bands... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);
    
    arma::vec Q_values = {-0.1, -0.05, 0.0, 0.05, 0.1};
    arma::rowvec Q = {0., 0., 0.};
    std::vector<std::vector<double>> energies;

    for (const auto& q : Q_values){
        Q(1) = q;
        xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55}, Q);

        exciton.brillouinZoneMesh(ncell);
        exciton.initializeHamiltonian();
        exciton.BShamiltonian();

        auto results = exciton.diagonalize("diag", nstates);
        auto resultingEnergies = xatu::detectDegeneracies(results->eigval, nstates, 6);

        energies.push_back(resultingEnergies[0]);
    }
    
    std::vector<std::vector<double>> expectedEnergies = {{1.810464, 2}, 
                                                         {1.780106, 2},
                                                         {1.768783, 2},
                                                         {1.780106, 2},
                                                         {1.810464, 2}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-5));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 exchange", "[MoS2-exchange]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 with exchange... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    arma::rowvec Q = {0., 0.1, 0.};

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55}, Q);
    exciton.setExchange(true);
    
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();

    auto results = exciton.diagonalize("diag", nstates);
    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.810464, 2}, 
                                                         {1.825740, 1},
                                                         {1.855544, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-5));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 reduced BZ", "[MoS2-reducedBZ]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 with reduced BZ... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.system->reducedBrillouinZoneMesh(ncell, factor);
    exciton.system->shiftBZ({0.6628, -1.1480, 0});

    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.795378, 1}, 
                                                         {1.809810, 1},
                                                         {1.935880, 1},
                                                         {1.950567, 1}};
    for(uint i = 0; i < energies.size(); i++){
        REQUIRE_THAT(energies[i][0], Catch::Matchers::WithinAbs(expectedEnergies[i][0], 1E-4));
        REQUIRE(energies[i][1] == expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST_CASE("MoS2 absorption", "[MoS2-kubo]"){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing MoS2 conductivity... ";
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config = xatu::SystemConfiguration(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 37.2792093358;
    REQUIRE_THAT(cum_norm_vme_ex, Catch::Matchers::WithinAbs(expectedTotalOscillator, 1E-7));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

