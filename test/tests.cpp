#include <gtest/gtest.h>
#include <iostream>
#include <armadillo>
#include <stdlib.h>
#include <string>
#include <vector>

#include <xatu.hpp>

#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif


TEST(FileParsing, ModelfileParsing){

    std::cout << std::setw(50) << std::left << "Testing model file parsing... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);
    
    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    int expectedFilling = 1;
    int expectedDim = 2;
    arma::urowvec expectedOrbitals = {1, 1};

    arma::vec expectedBravaisLatticeHash = {4.55422, 2.05422};
    arma::vec expectedMotifHash = {0, 3.90211};
    arma::vec expectedBravaisVectorsHash = {0, -0.110843, 2.38916, 2.05422, 4.55422};
    arma::vec expectedHamiltionianHash = {-1.34375, -0.9, -0.9, -0.0375, -0.0375};
    
    EXPECT_EQ(config.systemInfo.ndim, expectedDim);
    EXPECT_EQ(config.systemInfo.filling, expectedFilling);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        EXPECT_EQ(config.systemInfo.norbitals(i), expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        EXPECT_NEAR(hash, expectedBravaisLatticeHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        EXPECT_NEAR(hash, expectedMotifHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        EXPECT_NEAR(hash, expectedBravaisVectorsHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        EXPECT_NEAR(hash, expectedHamiltionianHash(i), 1E-4);
    }
}

TEST(FileParsing, CRYSTALfileParsing){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing CRYSTAL file parsing... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config(modelfile);

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
    
    EXPECT_EQ(config.systemInfo.ndim, expectedDim);
    EXPECT_EQ(config.systemInfo.filling, expectedFilling);
    EXPECT_EQ(config.systemInfo.hamiltonian.n_slices, expectedFockMatricesNumber);
    EXPECT_EQ(config.systemInfo.overlap.n_slices, expectedFockMatricesNumber);
    for(uint i = 0; i < config.systemInfo.norbitals.n_cols; i++){
        EXPECT_EQ(config.systemInfo.norbitals(i), expectedOrbitals(i));
    }

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        EXPECT_NEAR(hash, expectedBravaisLatticeHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        EXPECT_NEAR(hash, expectedMotifHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        EXPECT_NEAR(hash, expectedBravaisVectorsHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        EXPECT_NEAR(hash, expectedHamiltionianHash(i), 1E-2);
    } 
    for(uint i = 0; i < config.systemInfo.overlap.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.overlap.slice(i));
        EXPECT_NEAR(hash, expectedOverlapHash(i), 1E-2);
    }
}

TEST(FileParsing, W90fileParsing){
    // Suppress output during test
    std::cout << std::setw(50) << std::left << "Testing Wannier90 file parsing... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    // Load Wannier90 file
    std::string modelfile = "../examples/material_models/wannier/hBN_tb.dat";
    std::cout << "got file: " << modelfile << std::endl;
    xatu::Wannier90Configuration config(modelfile, 4);

    int expectedDim = 2;                         // Dimension of the system
    int expectedFilling = 4;                     // Number of filled states
    int expectedmSize = 6;                       // Number of unit cells
    int expectednFock = 83;                       // Number of Fock matrices
    
    // Hashes for real-space lattice vectors (flattened 3x3 matrix)
    arma::vec expectedBravaisLatticeHash = {2.425222, 4.080734};
    
    arma::vec expectedBravaisVectorsHash = {-2.971610, -2.456147, -1.940685, -3.584157, -3.068694, -2.553232, -2.037769, -1.522307,
    -1.340178,  2.407223, -4.014574, -3.165779, -2.650316, -2.134854, -1.619392, -1.103929,
    -0.921800,  2.825601,  6.239668, -1.910645, -2.565272, -1.716476, -1.201014, -0.685551,
    -0.503422,  3.243979,  6.658046, 10.072114,  0.181244, -0.140049, -0.461343, -1.115969,
    -0.267174, -0.085044,  3.662356,  7.076424, 10.490492,  1.951840,  1.630546,  1.309253,
     0.987960,  0.000000,  4.080734,  7.494802, 10.908869, 14.322937,  3.722436,  3.401142,
     3.079849,  2.425222,  5.335868,  7.579846, 11.327247, 14.741315, 18.155382,  5.814325,
     5.493031,  5.171738,  4.517112,  7.427757, 10.005069, 12.582381, 14.826359, 18.573760,
     7.584920,  7.263627,  6.609001,  9.519646, 12.096958, 14.674270, 17.251582, 19.828894,
    22.072872,  9.355516,  8.700890, 11.611535, 14.188847, 16.766159, 19.343471, 21.920783,
    13.703424, 16.280736, 18.858048
};

    // Hashes for motif positions (flattened matrix)
    arma::vec expectedMotifHash = {2.566478, 5.266827, 5.661396, 10.201080, 13.255683, 14.551132};

    // arma::vec expectedDegenHash = {1, 1, 2, 1, 1};

    // Hashes for Fock matrices (real/imag parts handled separately)
    arma::vec expectedHamiltonianHash = {0.978634, 0.990920, 1.046875, 1.054984, 0.996068, 0.968801, 1.101170, 0.961007,
        0.968796, 1.043775, 1.032634, 1.017193, 0.940727, 1.172638, 0.890277, 1.110180,
        1.088959, 1.027729, 0.962491, 0.970735, 1.036763, 1.089453, 1.079106, -0.531977,
        0.946280, 0.965448, 1.031659, 1.011340, 1.016294, 1.009594, 0.999775, 0.951863,
        0.458796, 0.579032, 1.125541, 1.059755, 0.985202, 0.986659, 1.042266, 1.042844,
        2.030694, -11.551247, 2.313774, 1.014736, 1.049511, 0.994099, 0.997162, 1.026182,
        1.131210, 2.914045, 6.239713, 0.937042, 1.013736, 0.995272, 1.019086, 0.998125,
        1.032483, 0.973986, 1.026194, 1.357262, 1.224249, 1.032389, 1.027929, 0.987408,
        0.982629, 1.003563, 1.071731, 1.044472, 0.880948, 1.118926, 0.963565, 0.998474,
        1.029212, 1.046875, 0.982783, 0.964047, 1.057069, 0.999934, 0.973145, 1.045742,
        1.043775, 1.005556, 0.966527
    };
    
    // Hashes for Rhop matrices (if used)
    arma::vec expectedRhopHash = {0.997114, 0.990306, 0.997694, 0.991232, 1.005100, 1.004297, 0.989814, 0.977408,
    1.005649, 0.983461, 1.006747, 1.007933, 1.012504, 0.994354, 0.991596, 1.013520,
    1.007218, 1.005119, 0.996115, 0.995221, 1.012035, 1.011999, 1.006950, 0.997909,
    1.000703, 0.989825, 1.002780, 1.002820, 0.985856, 1.007615, 0.986269, 0.993846,
    1.007167, 1.002167, 0.989714, 0.993715, 1.014779, 1.004536, 0.999528, 1.020440,
    1.012040, 1.008106, 0.994788, 1.035806, 0.989133, 0.991763, 0.972419, 0.991202,
    1.022495, 1.038578, 1.006076, 1.005136, 1.000803, 0.992482, 0.990968, 1.014487,
    1.003521, 1.004732, 0.985285, 1.003408, 0.988708, 0.953025, 1.004636, 1.021256,
    1.001702, 0.998779, 1.022129, 0.971661, 1.021564, 1.088059, 1.024697, 1.026738,
    1.006775, 0.963720, 0.993811, 1.026752, 1.029950, 1.005488, 0.957017, 1.021292,
    1.008024, 0.997984, 1.019589, 0.994519, 1.053920, 1.020668, 1.007176, 0.985443,
    1.008066, 0.993171, 1.007998, 0.976088, 1.000136, 0.992240, 0.992790, 1.035847,
    1.024869, 1.066353, 1.095552, 1.019990, 1.021736, 1.087952, 0.970075, 0.995309,
    1.004839, 1.015675, 1.005783, 0.995738, 0.980129, 0.997494, 1.004768, 0.983212,
    0.996325, 1.002814, 1.005364, 1.001616, 1.002397, 0.954128, 1.003217, 1.018674,
    1.085609, 1.078255, 1.012354, 4.548316, 2.873857, 1.197828, 1.115675, 1.006616,
    0.964120, 0.984689, 0.978937, 1.024818, 1.025584, 1.003606, 1.005235, 1.017286,
    0.988046, 0.995220, 0.996773, 0.999681, 0.996953, 1.013851, 0.991732, 1.009270,
    0.957047, 1.006585, 0.995265, 0.992805, 0.997985, 1.040597, 0.972494, 0.999376,
    1.115256, 1.006251, 0.996555, 1.018850, 1.021810, 0.972551, 0.998658, 0.984656,
    1.009052, 1.005229, 1.040141, 1.028630, 1.004017, 1.007448, 1.002345, 1.001194,
    0.991527, 1.005624, 1.007665, 1.025050, 1.006177, 1.011064, 0.970850, 0.992055,
    0.981971, 1.049555, 1.007129, 1.029757, 1.030265, 0.961084, 0.987074, 1.007111,
    0.998259, 1.021026, 0.994860, 1.013045, 1.006976, 0.999025, 1.012405, 0.992590,
    1.009446, 0.996579, 0.992565, 1.007522, 0.998227, 1.001843, 1.009624, 1.028827,
    1.022804, 1.022722, 1.002362, 0.977078, 0.983739, 1.012481, 0.999343, 1.005116,
    1.003181, 1.022401, 0.997381, 0.990140, 0.993093, 1.004442, 1.016081, 0.998645,
    0.993441, 1.001894, 1.008756, 0.989814, 0.977408, 1.005649, 1.006771, 1.007816,
    0.996179, 0.981866, 1.020622, 1.012108, 1.056906, 1.033585, 1.020818, 1.011852,
    0.997104, 1.001710, 1.015061, 0.997085, 1.002824, 0.989569, 1.023114, 1.005145,
    1.002820, 0.985856, 1.007615, 1.010556, 1.008685, 0.990653, 1.017227, 1.020576,
    1.002143
    };

    // --- Validation ---
    EXPECT_EQ(config.ndim, expectedDim);
    EXPECT_EQ(config.filling, expectedFilling);
    EXPECT_EQ(config.mSize, expectedmSize);
    EXPECT_EQ(config.nFock, expectednFock);

    for(uint i = 0; i < config.systemInfo.bravaisLattice.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisLattice.row(i));
        EXPECT_NEAR(hash, expectedBravaisLatticeHash(i), 1E-4);
    }

    // double degenHash = xatu::array2hash(arma::vectorise(config.degeneracies));
    // EXPECT_NEAR(degenHash, expectedDegenHash), 1E-4);

    for(uint i = 0; i < config.systemInfo.motif.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.motif.row(i));
        EXPECT_NEAR(hash, expectedMotifHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.bravaisVectors.n_rows; i++){
        double hash = xatu::array2hash(config.systemInfo.bravaisVectors.row(i));
        EXPECT_NEAR(hash, expectedBravaisVectorsHash(i), 1E-4);
    }
    for(uint i = 0; i < config.systemInfo.hamiltonian.n_slices; i++){
        double hash = xatu::array2hash(config.systemInfo.hamiltonian.slice(i));
        EXPECT_NEAR(hash, expectedHamiltonianHash(i), 1E-2);
    } 

    uint hashIndex = 0;
    for(uint i = 0; i < config.Rhop.n_elem; i++) {
        for(uint s = 0; s < config.Rhop(i).n_slices; s++) {
            arma::cx_mat slice = config.Rhop(i).slice(s);
            double sliceHash = xatu::array2hash(slice);
            EXPECT_NEAR(sliceHash, expectedRhopHash(hashIndex), 1E-2);
            hashIndex++;
        }
    }
}

TEST(FileParsing, ExcitonfileParsing){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing exciton file parsing... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    std::string excitonconfig = "./data/hBN_spinless.txt";
    xatu::ExcitonConfiguration config(excitonconfig);

    std::string expectedName = "hBN_N30";
    int expectedNcell = 30;
    int expectedNbands = 1;
    arma::vec expectedDielectric = {1, 1, 10};

    EXPECT_EQ(config.excitonInfo.ncell, expectedNcell);
    EXPECT_EQ(config.excitonInfo.nbands, expectedNbands);
    for (uint i = 0; i < expectedDielectric.n_elem; i++){
        EXPECT_EQ(config.excitonInfo.eps(i), expectedDielectric(i));
    }
    EXPECT_EQ(xatu::array2hash(config.excitonInfo.Q), 0);
}

TEST(FileParsing, ExcitonfileParsingAnisotropic){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing anisotropic exciton file parsing... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    std::string excitonconfig = "./data/hBN_anisotropic.txt";
    xatu::ExcitonConfiguration config(excitonconfig);

    std::string expectedName = "hBN_ani";
    int expectedNcell = 30;
    int expectedNbands = 1;
    arma::vec expectedDielectric = {1, 1, 10, 25, 35};

    EXPECT_EQ(config.excitonInfo.ncell, expectedNcell);
    EXPECT_EQ(config.excitonInfo.nbands, expectedNbands);
    for (uint i = 0; i < expectedDielectric.n_elem; i++){
        EXPECT_EQ(config.excitonInfo.eps(i), expectedDielectric(i));
    }
    EXPECT_EQ(xatu::array2hash(config.excitonInfo.Q), 0);

}

TEST(FileParsing, ScreeningParsing){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing screening file parsing... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string screeningfile = "../examples/screeningconfig/hBN_DFT_screening.txt";
    xatu::ScreeningConfiguration config = xatu::ScreeningConfiguration(screeningfile);

    std::string expectedfunction = "dielectric";
    int expectedNcell_aux = 10;
    int expected_valence_bands = 6;
    int expected_conduction_bands = 63;
    bool expectedspin = false;
    bool expectedisotropic = true;
    double expectedGcutoff = 3.0;
    arma::ivec expectedvectors = arma::ivec({0, 0});
    arma::rowvec expectedmomentum = {0.2, 0., 0.};

    EXPECT_EQ(config.screeningInfo.ncell_aux, expectedNcell_aux);
    EXPECT_EQ(config.screeningInfo.nvbands, expected_valence_bands);
    EXPECT_EQ(config.screeningInfo.ncbands, expected_conduction_bands);
    EXPECT_EQ(config.screeningInfo.spin, expectedspin);
    EXPECT_EQ(config.screeningInfo.isotropic, expectedisotropic);
    EXPECT_EQ(config.screeningInfo.Gcutoff, expectedGcutoff);
    EXPECT_EQ(config.screeningInfo.Gs(0), expectedvectors(0));
    EXPECT_EQ(config.screeningInfo.Gs(1), expectedvectors(1));
    EXPECT_EQ(config.screeningInfo.q(0), expectedmomentum(0));
    EXPECT_EQ(config.screeningInfo.q(1), expectedmomentum(1));
    EXPECT_EQ(config.screeningInfo.q(2), expectedmomentum(2));

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}
    
TEST(hBNTest, FullDiagonalization){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN energies (fulldiag)... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(hBNTest, DavidsonDiagonalization){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN energies (Davidson)... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("davidson", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(hBNTest, LanczosDiagonalization){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN energies (Lanczos)... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("sparse", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 2}, 
                                                         {6.074062, 1}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(hBNTest, ReciprocalSpaceMethod){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN energies (reciprocal)... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.setMode("reciprocalspace");
    exciton.setGcutoff(3.0);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{6.200045075, 1}, 
                                                         {6.300352, 1},
                                                         {6.832558, 1}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(hBNTest, ReciprocalWavefunction){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN reciprocal w.f... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

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
    EXPECT_NEAR(kwfHash, expectedKwfHash, 1E-5);

}

TEST(hBNTest, RealSpaceWavefunction){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN real-space w.f... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

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
    EXPECT_NEAR(rswfHash, expectedRSwfHash, 1E-5);

}

TEST(hBNTest, Conductivity){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN conductivity... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 47.4140064;
    EXPECT_NEAR(cum_norm_vme_ex, expectedTotalOscillator, 1E-6);

}

TEST(hBNTest, SpinfulSpectrumCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN spinful... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 12;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{5.335690, 8}, 
                                                         {6.074062, 4}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}
    
TEST(hBNTest, SpinProjectionCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN spin... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 16;
    int nstates = 4;

    std::string modelfile = "../examples/material_models/hBN_spinful.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1, 1, 10});
    
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::mat expectedSpin = arma::vec{-1, 0, 0, 1};  

    for(uint i = 0; i < nstates; i++){
        arma::cx_vec spin = results->spinX(i);
        EXPECT_EQ(spin.n_elem, 3);
        EXPECT_NEAR(std::real(spin(0)), expectedSpin(i), 1E-2);
        EXPECT_NEAR(std::abs(spin(1)), 0.5, 1E-2);
        EXPECT_NEAR(std::abs(spin(2)), 0.5, 1E-2);
    }

}

TEST(hBNTest, DFTCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing DFT hBN... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;
    bool triangular = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config(modelfile, 100);

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
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
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
    EXPECT_NEAR(kwfHash, expectedKwfHash, 1E-5);

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
    EXPECT_NEAR(rswfHash, expectedRSwfHash, 1E-5);

}

TEST(hBNTest, WannierCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing Wannierized hBN... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 30;
    int nstates = 2;
    bool triangular = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    std::string modelfile = "../examples/material_models/wannier/hBN_tb.dat";
    xatu::Wannier90Configuration config(modelfile, 4);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10});
    exciton.system->setAU(true);

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{2.833347, 1}, 
                                                         {2.834596, 1}};
                                                         
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }

    // check oscillator strength for absorption
    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 103.0325881;
    EXPECT_NEAR(cum_norm_vme_ex, expectedTotalOscillator, 1E-6);

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
    double expectedKwfHash = 1.057666;
    EXPECT_NEAR(kwfHash, expectedKwfHash, 1E-5);

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
    double expectedRSwfHash = 1.003540;
    EXPECT_NEAR(rswfHash, expectedRSwfHash, 1E-5);

}

TEST(hBNTest, AnisotropicSpectrumCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN-ani energies (fulldiag)... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 3;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10, 25, 35});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{6.428535, 1}, 
                                                         {6.456359, 1},
                                                         {6.732466, 1}};

    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(hBNTest, AnisotropicConductivity){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing TB hBN anisotropic conductivity... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 20;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/hBN.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 1, 0, {1, 1, 10, 25, 35});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 47.4140064;
    EXPECT_NEAR(cum_norm_vme_ex, expectedTotalOscillator, 1E-6);

}

TEST(MoS2Test, FullDiagonalization){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 energies... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{1.768783, 2}, 
                                                         {1.780562, 2}};
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(MoS2Test, ReciprocalWavefunction){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 reciprocal w.f... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

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
    EXPECT_NEAR(kwfHash, expectedKwfHash, 1E-5);

}

TEST(MoS2Test, SpinProjectionCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 spin... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

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
            EXPECT_NEAR(spinValue, expectedSpin(i, j), 1E-4);
        }
    }

}

TEST(MoS2Test, ExcitonBands){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 exciton bands... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);
    
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
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-5);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(MoS2Test, ExchangeCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 with exchange... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    arma::rowvec Q = {0., 0.1, 0.};

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

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
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-5);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(MoS2Test, ReducedBZCalculation){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 with reduced BZ... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 4;
    int factor = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

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
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }
}

TEST(MoS2Test, Conductivity){

    std::cout.clear();
    std::cout << std::setw(50) << std::left << "Testing MoS2 conductivity... " << std::endl;
    std::cout.setstate(std::ios_base::failbit);

    int ncell = 12;
    int nstates = 2;

    std::string modelfile = "../examples/material_models/MoS2.model";    
    xatu::SystemConfiguration config(modelfile);

    xatu::ExcitonTB exciton = xatu::ExcitonTB(config, ncell, 2, 0, {1., 4., 13.55});

    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();
    exciton.BShamiltonian();
    auto results = exciton.diagonalize("diag", nstates);

    arma::cx_mat vme_ex = results->excitonOscillatorStrength();
    arma::mat norm_vme_ex = arma::square(arma::abs(vme_ex));
    double cum_norm_vme_ex = arma::accu(norm_vme_ex);

    double expectedTotalOscillator = 37.2822837;
    EXPECT_NEAR(cum_norm_vme_ex, expectedTotalOscillator, 1E-6);

}

TEST(ScreeningTest, hBNDFTScreening){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing DFT hBN screening (polarizability, direct dielectric matrix, and inverse dielectric matrix)... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    xatu::CRYSTALConfiguration config = xatu::CRYSTALConfiguration(modelfile, 100);

    std::unique_ptr<xatu::SystemConfiguration> systemConfig = std::unique_ptr<xatu::CRYSTALConfiguration>(new xatu::CRYSTALConfiguration(modelfile, 100));

    int ncell = 10;
    double Gcutoff = 3.0;
    double d = 3.33;
    double percentage = 0.5;
    int nstates = 2;
    uint nvbands = 6;
    uint ncbands = 63;
    uint ncell_aux = 10;
    bool spinfull = false;
    bool isotropic = true;

    xatu::ExcitonTB exciton = xatu::ExcitonTB(*systemConfig, ncell, 1, 0, {1., 1., 10.}, {0., 0., 0.}, ncell_aux, nvbands, ncbands, Gcutoff, Gcutoff, d, spinfull, isotropic);
    
    exciton.system->setAU(true);
    exciton.brillouinZoneMesh(ncell);
    exciton.initializeHamiltonian();

    std::complex<double> expectedepsilon = 2.001461924;

    std::complex<double> epsilon = exciton.computesingleDielectricFunctionMatrixElement();

    EXPECT_NEAR(std::real(epsilon), std::real(expectedepsilon), 1E-4);

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;
}

TEST(ScreeningTest, hBNDFTExciton){

    std::cout.clear();
    std::cout << std::setw(40) << std::left << "Testing DFT hBN exciton calculation with screening... ";
    std::cout.setstate(std::ios_base::failbit);

    std::string modelfile = "../examples/material_models/DFT/hBN_base_HSE06.outp";
    std::unique_ptr<xatu::SystemConfiguration> systemConfig = std::unique_ptr<xatu::CRYSTALConfiguration>(new xatu::CRYSTALConfiguration(modelfile, 100));
    std::string dielectric_matrix_path = "../data/hBN_DFT_HSE06_Q2D_BZ_N20_Gc_5_1_invepsilon.dat";

    int ncell = 20;
    double Gcutoff = 5.1;
    double d = 3.33;
    double percentage = 0.5;
    int nstates = 2;
    uint nvbands = 6;
    uint ncbands = 63;
    uint ncell_aux = 10;
    bool spinfull = false;
    bool isotropic = true;
    int holeIndex = 1;
    arma::rowvec holeCell = {0, 0, 0};

    xatu::ExcitonTB exciton = xatu::ExcitonTB(*systemConfig, ncell, 1, 0, {1., 1., 10.}, {0., 0., 0.}, ncell_aux, nvbands, ncbands, Gcutoff, Gcutoff, d, spinfull, isotropic);
    
    exciton.system->setAU(true);
    exciton.setPercentage(percentage);

    exciton.brillouinZoneMesh(ncell);
    exciton.readInverseDielectricMatrix(dielectric_matrix_path);
    exciton.initializeHamiltonian();

    exciton.BShamiltonian();

    auto results = exciton.diagonalize("diag", nstates);

    auto energies = xatu::detectDegeneracies(results->eigval, nstates, 6);
    
    std::vector<std::vector<double>> expectedEnergies = {{3.696992, 1}, 
                                                         {3.730510, 1}};
                                                         
    for(uint i = 0; i < energies.size(); i++){
        EXPECT_NEAR(energies[i][0], expectedEnergies[i][0], 1E-4);
        EXPECT_EQ(energies[i][1], expectedEnergies[i][1]);
    }

    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;

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
    double expectedKwfHash = 1.086489570724429;
    if(std::abs(kwfHash-expectedKwfHash) < 1E-5){
        std::cout << "\033[1;32m Reciprocal w.f. correct \033[0m" << std::endl;
    }
    std::cout << std::setw(40) << "kwfHash: " << kwfHash << std::endl;

    // Check realspace w.f.
    arma::rowvec holePosition = exciton.system->motif.row(holeIndex).subvec(0, 2) + holeCell;
    std::cout << "Hole position: " << holePosition << std::endl;
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
    double expectedRSwfHash = 91.4518564389955063;
    if(std::abs(rswfHash-expectedRSwfHash) < 1E-5){
        std::cout << "\033[1;32m Real space w.f. correct \033[0m" << std::endl;
    }
    std::cout << std::setw(40) << "rswfHash: " << rswfHash << std::endl; 
    std::cout.clear();
    std::cout << std::setw(40) << "\033[1;32m Success \033[0m" << std::endl;

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}