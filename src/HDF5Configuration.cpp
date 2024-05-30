#include "xatu/HDF5Configuration.hpp"

namespace xatu {

/**
 * File constructor for HDF5Configuration. 
 * The provided filename corresponds to an HDF5 file, which is then 
 * parsed with armadillo to extract the content.
 * @param filename Name of the HDF5 file to be parsed.
*/
HDF5Configuration::HDF5Configuration(const std::string& filename) : ConfigurationBase{filename} {
    #ifdef ARMA_USE_HDF5
        this->filename = filename;
        parseContent();
        checkContentCoherence();
    #else
        throw std::runtime_error("Xatu was not compiled with HDF5 support");
    #endif
}

/**
 * Method responsible of extracting the information from the HDF5 file.
*/
void HDF5Configuration::parseContent(){

    #ifdef ARMA_USE_HDF5
    if (!systemInfo.bravaisLattice.load(arma::hdf5_name(filename, "bravaisLattice"))){
        throw std::runtime_error("Error reading HDF5 file: 'bravaisLattice' not found");
    }
    systemInfo.ndim = systemInfo.bravaisLattice.n_rows;
    
    if (!systemInfo.motif.load(arma::hdf5_name(filename, "motif"))){
        throw std::runtime_error("Error reading HDF5 file: 'motif' not found");
    }
    if (!systemInfo.norbitals.load(arma::hdf5_name(filename, "norbitals"))){
        throw std::runtime_error("Error reading HDF5 file: 'norbitals' not found");
    };

    if (!systemInfo.hamiltonian.load(arma::hdf5_name(filename, "hamiltonian"))){
        throw std::runtime_error("Error reading HDF5 file: 'hamiltonian' not found");
    }
    if (!systemInfo.bravaisVectors.load(arma::hdf5_name(filename, "bravaisVectors"))){
        throw std::runtime_error("Error reading HDF5 file: 'unitCellVectors' not found");
    }
    for(unsigned int i = 0; i < systemInfo.bravaisVectors.n_rows; i++){
        arma::rowvec cartersian_vector = arma::zeros<arma::rowvec>(3);
        for (int j = 0; j < systemInfo.ndim; j++){
            cartersian_vector += systemInfo.bravaisVectors(i, j) * systemInfo.bravaisLattice.row(j);
        }
        systemInfo.bravaisVectors.row(i) = cartersian_vector;
    }

    if (!systemInfo.overlap.load(arma::hdf5_name(filename, "overlap"))){
        throw std::runtime_error("Error reading HDF5 file: 'overlap' not found");
    };

    arma::vec filling_vec;
    if (!filling_vec.load(arma::hdf5_name(filename, "filling"))){
        throw std::runtime_error("Error reading HDF5 file: 'filling' not found");
    }
    systemInfo.filling = filling_vec(0);

    arma::cx_cube hamiltonian_imag;
    if (hamiltonian_imag.load(arma::hdf5_name(filename, "hamiltonian_imag"))){
        systemInfo.hamiltonian += hamiltonian_imag;
    }

    #endif
}

}