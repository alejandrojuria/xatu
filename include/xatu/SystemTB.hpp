#pragma once
#include "xatu/utils.hpp"
#include "xatu/System.hpp"

#ifndef constants
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#define auToEV 27.211386
#endif

namespace xatu {

/**
 * The SystemTB class contains all information regarding the system where we want to compute
 * the exciton spectrum in the TB mode. It is defined as a sub-class of System
*/
class SystemTB : public System {

    protected:
        bool isTriangular_ = false;
        bool isCRYSTAL_ = false;
        arma::mat motif_;
        int natoms_;
        arma::urowvec orbitalsPerSpecies_;
        double a_, c_;
        int factor_ = 1;
        arma::mat meshBZ_;
        arma::mat inverseReciprocalMatrix_;
        arma::cx_cube hamiltonianMatrices_Zeeman_;

    public:  // Const references to attributes (read-only)
        // Returns true if the system matrices (hamiltonian and overlap) are triangular
        const bool& isTriangular = isTriangular_;
        // Returns true if the system has been parsed from a CRYSTAL file (and, in particular, has the Hamiltonian in atomic units)
        const bool& isCRYSTAL = isCRYSTAL_;
        // Matrix containing the positions of the atoms of the motif by columns: (4,natoms). Each row has the format {x,y,z; species},
        // where x,y,z are in Angstrom and the last element is an index (0,1,...) labelling the chemical species
        const arma::mat& motif = motif_;
        // Number of atoms in the unit cell
        const int& natoms = natoms_;
        // Vector where each entry is a chemical species, and the value is its number of orbitals
        const arma::urowvec& orbitalsPerSpecies = orbitalsPerSpecies_;
        // Lattice parameters
        const double& a = a_;
        // Lattice parameters
        const double& c = c_;
        // Reduction factor of Brillouin zone mesh
        const int& factor = factor_;
        // Mesh of the Brillouin zone, stored by (3D) columns
        const arma::mat& meshBZ = meshBZ_;
        // Matrix that maps any kpoint back to the original Monkhorst-Pack mesh
        const arma::mat& inverseReciprocalMatrix = inverseReciprocalMatrix_;
        // Hamiltonian with an on-site Zeeman term, initialized only with SystemTB::addZeeman
        const arma::cx_cube& hamiltonianMatrices_Zeeman = hamiltonianMatrices_Zeeman_;
    

    protected:
        SystemTB() = default;
    public:
        SystemTB(const ConfigurationSystem&);  
        SystemTB(const SystemTB&) = default;

    public:
        /* Setters */
        void setCRYSTAL(const bool isCRYSTAL);
        void setFilling(const int filling);

        /* Bloch Hamiltonian */
        arma::cx_mat hamiltonian(const arma::colvec& k) const override;
        arma::cx_mat overlap(const arma::colvec& k) const override;
        using System::solveBands;
        void solveBands(const arma::colvec&, arma::vec&, arma::cx_mat&) const override;

        /* Gauge */
        arma::cx_vec latticeToAtomicGauge(const arma::cx_vec&, const arma::colvec&);
        arma::cx_vec atomicToLatticeGauge(const arma::cx_vec&, const arma::colvec&);

        /* Observables */
        double expectedSpinZValue(const arma::cx_vec&);
        // double expectedSpinYValue(const arma::cx_vec&);  // commented out until implemented
        // double expectedSpinXValue(const arma::cx_vec&);  // commented out until implemented   
        arma::cx_vec velocity(const arma::colvec&, int fBand, int sBand) const;  

        /* Modifiers */
        void addZeeman(const double amplitude);  

        /* Brillouin zone meshing, supercells & utilities */
        void brillouinZoneMesh(int);
        void reducedBrillouinZoneMesh(const int, const int);
        void shiftBZ(const arma::colvec&);
        void calculateInverseReciprocalMatrix();
        int findEquivalentPointBZ(const arma::colvec&, int);
        arma::mat truncateSupercell(const int, const double);
        arma::mat truncateReciprocalSupercell(const int, const double);
        void extractLatticeParameters();
        arma::rowvec rotateC3(const arma::rowvec&);
    
    private:
        void orthogonalize_hamiltonian(const arma::colvec&, arma::cx_mat&) const;

};

}
