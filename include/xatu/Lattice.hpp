#pragma once
#include <armadillo>
#include <string>
#include <iostream>
#include "xatu/SystemConfiguration.hpp"


#ifndef constants
#define PI 3.141592653589793
#define ec 1.6021766E-19
#define eps0 8.8541878E-12
#endif

namespace xatu {

/// @brief The Lattice class is designed to hold all information regarding both
/// the Bravais lattice and the reciprocal lattice of the system.
class Lattice {
    
    //// Attributes
    protected:
        int ndim_, natoms_, nk_, ncells_;
        int factor_ = 1;
        double a_, c_, unitCellArea_;
        arma::mat bravaisLattice_, motif_, unitCellList_;
        arma::mat reciprocalLattice_, kpoints_, meshBZ_;
        std::map<std::string, int> atomToIndex;
        arma::mat inverseReciprocalMatrix;

    // Const references to attributes (read-only)
    public:
        // Returns system dimension
        const int& ndim = ndim_;
        // Returns number of atoms in the motif
        const int& natoms = natoms_;
        // Returns Bravais lattice vectors of system
        const arma::mat& bravaisLattice = bravaisLattice_;
        // Returns matrix containing the positions of the atoms of the motif
        // by rows. Each row has the format {x,y,z; species}, where the 
        // last elements is an index representing the chemical species.
        // const arma::mat& motif = motif_;
        const arma::mat& motif = motif_;
        // Returns list of Bravais vectors whose atom connect with the unit cell at the origin
        const arma::mat& unitCellList = unitCellList_;
        // Returns reciprocal lattice vectors
        const arma::mat& reciprocalLattice = reciprocalLattice_;
        // Returns kpoints
        const arma::mat& kpoints = kpoints_;
        // Number of k points
        const int& nk = nk_;
        // Lattice parameters
        const double& a = a_;
        // Lattice parameters
        const double& c = c_;
        // Unit cell area
        const double& unitCellArea = unitCellArea_;
        // Number of unit cells connected to the origin
        const int& ncells = ncells_;

    //// Methods
    protected:
        Lattice(){}; // Protected so that Lattice can not be initialized (abstract)
    public:
        Lattice(const Lattice&); // Copy constructor
        ~Lattice(){};

        int getDimension();
        int getNumAtoms();
        arma::mat getBravaisLattice();
        arma::mat getMotif();
        arma::mat getUnitCellList();
        arma::mat getReciprocalLattice();

        /* Mesh generation routines */
        void brillouinZoneMesh(int);
        arma::mat brillouinZoneMeshCrystalCoordinates(int);
        void reducedBrillouinZoneMesh(int, int);
        void shiftBZ(const arma::rowvec&);
        void preserveC3();
        void brillouinZoneC3Mesh(int);
        arma::mat wignerSeitzSupercell(int);
        arma::mat truncateSupercell(int, double);
        arma::mat truncateReciprocalSupercell(int, double);
        arma::mat generateCombinations(int n, int ndim, bool centered = false);
        arma::mat supercellCutoff(int);

        /* Crystal operations */
        arma::cx_mat inversionOperator(const arma::cx_vec&);
        arma::rowvec rotateC3(const arma::rowvec&);
        int findEquivalentPointBZ(const arma::rowvec&, int);
        void calculateInverseReciprocalMatrix();


    protected:
        void initializeLatticeAttributes(const SystemConfiguration&);
        void extractLatticeParameters();
        void computeUnitCellArea();
        void calculateReciprocalLattice();
        bool isInsideWsCell(const arma::rowvec&, const arma::mat&, 
                            const arma::rowvec&);
        
};

}
