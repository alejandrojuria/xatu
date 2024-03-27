#include <armadillo>

#include "xatu.hpp"

// TO DO, NOTA BENE: This file is currently unused. It is aimed at extracting the calculation
// of the interaction matrix elements out of the exciton class to be available for
// a wider class of calculations that involve V.


namespace xatu {

// ---------------------------- Potentials ----------------------------

/** 
 * Calculate value of interaction potential (Keldysh). Units are eV.
 * @details If the distance is zero, then the interaction is renormalized to be V(a) since
 * V(0) is infinite, where a is the lattice parameter. Also, for r > cutoff the interaction is taken to be zero.
 * @param r Distance at which we evaluate the potential.
 * @param r0 Screening distance.
 * @param eps_s Dielectric constant of substrate.
 * @param eps_m Dielectric constant of embedding medium (Air = 1).
 * @param cutoff Cutoff distance in potential, V(r > Rc) = 0.
 * @param a Regularization distance, V(0) = V(a).
 * @return Value of Keldysh potential, V(r).
 */
double keldysh(double r, double r0, double eps_s, double eps_m, double cutoff, double a){

    double eps_bar = (eps_m + eps_s)/2;
    double SH0;
    double R = abs(r)/r0;
    double potential_value;
    if(r == 0){
        STVH0(a/r0, &SH0);
        potential_value = ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(a/r0));
    }
    else if (r > cutoff){
        potential_value = 0.0;
    }
    else{
        STVH0(R, &SH0);
        potential_value = ec/(8E-10*eps0*eps_bar*r0)*(SH0 - y0(R));
    };

    return potential_value;

}

/** 
 * Purpose: Compute Struve function H0(x).
 * Source: http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mstvh0_cpp.txt 
 * @param X x --- Argument of H0(x) ( x ò 0 )
 * @param SH0 SH0 --- H0(x). The return value is written to the direction of the pointer.
*/
void STVH0(double X, double *SH0) {
    double A0,BY0,P0,Q0,R,S,T,T2,TA0;
	int K, KM;

        S=1.0;
        R=1.0;
        if (X <= 20.0) {
           A0=2.0*X/PI;
           for (K=1; K<61; K++) {
              R=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
    e15:       *SH0=A0*S;
        }
        else {
           KM=int(0.5*(X+1.0));
           if (X >= 50.0) KM=25;
           for (K=1; K<=KM; K++) {
              R=-R*pow((2.0*K-1.0)/X,2);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
    e25:       T=4.0/X;
           T2=T*T;
           P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
           Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
           TA0=X-0.25*PI;
           BY0=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
           *SH0=2.0/(PI*X)*S+BY0;
        }
}

/**
 * Evaluates the Fourier transform of the Keldysh potential, which is an analytical expression.
 * @param q kpoint where we evaluate the FT.
 * @param r0 Screening distance.
 * @param eps_s Dielectric constant of substrate.
 * @param eps_m Dielectric constant of embedding medium.
 * @param unitCellArea Area of unit cell.
 * @param totalCells Number of unit cells of the system.
 * @return Fourier transform of the potential at q, FT[V](q).
 */
double keldyshFT(const arma::rowvec& q, double r0, double eps_s, double eps_m, double unitCellArea, int totalCells){

    double potential = 0;
    double eps_bar = (eps_m + eps_s)/2;

    double qnorm = arma::norm(q);
    double eps = 1E-8;
    if (qnorm < eps){
        potential = 0;
    }
    else{
        potential = 1/(qnorm*(1 + r0*qnorm));
    }
    
    potential = potential*ec*1E10/(2*eps0*eps_bar*unitCellArea*totalCells);
    return potential;
}

/**
 * Coulomb potential in real space.
 * @param r Distance at which we evaluate the potential.
 * @param regularization Regularization distance to remove divergence at r=0.
 * @return Value of Coulomb potential, V(r).
*/
double coulomb(double r, double regularization){
    return (r > 1E-8) ? ec/(4E-10*PI*eps0*r) : ec*1E10/(4*PI*eps0*regularization);    
}

// ----------------------- Interaction matrix element initialization -----------------------

/**
 * Routine to compute the lattice Fourier transform with the potential displaced by some
 * vectors of the motif.
 * @param firstAtom Vector of first atom.
 * @param secondAtom Vector of second atom.
 * @param k kpoint where we evaluate the FT.
 * @param cells Matrix with the unit cells over which we sum to compute the lattice FT.
 * @param totalCells Number of unit cells of the system.
 * @return Motif lattice Fourier transform of the Keldysh potential at k.
 */
std::complex<double> motifFourierTransform(const arma::rowvec& firstAtom, const arma::rowvec& secondAtom, 
                                           const arma::rowvec& k, const arma::mat& cells, int totalCells,
                                           double r0, double eps_s, double eps_m, double cutoff, double a){

    std::complex<double> imag(0,1);
    std::complex<double> Vk = 0.0;

    for(int n = 0; n < cells.n_rows; n++){
        arma::rowvec cell = cells.row(n);
        double module = arma::norm(cell + firstAtom - secondAtom);
        Vk += keldysh(module, r0, eps_s, eps_m, cutoff, a)*std::exp(imag*arma::dot(k, cell));
    }
    Vk /= pow(totalCells, 1);

    return Vk;
}

/**
 * Method to extend the motif Fourier transform matrix to match the dimension of the
 * one-particle basis. 
 * @param motifFT Matrix storing the motif Fourier transform to be extended.
 * @return Extended matrix.
 */
arma::cx_mat extendMotifFT(const arma::cx_mat& motifFT, int basisdim, const arma::mat& motif, const arma::urowvec& orbitals){
    arma::cx_mat extendedMFT = arma::zeros<arma::cx_mat>(basisdim, basisdim);
    int rowIterator = 0;
    int colIterator = 0;
    for(unsigned int atom_index_r = 0; atom_index_r < motif.n_rows; atom_index_r++){
        int species_r = motif.row(atom_index_r)(3);
        int norbitals_r = orbitals(species_r);
        colIterator = 0;
        for(unsigned int atom_index_c = 0; atom_index_c < motif.n_rows; atom_index_c++){
            int species_c = motif.row(atom_index_c)(3);
            int norbitals_c = orbitals(species_c);
            extendedMFT.submat(rowIterator, colIterator, 
                               rowIterator + norbitals_r - 1, colIterator + norbitals_c - 1) = 
                          motifFT(atom_index_r, atom_index_c) * arma::ones(norbitals_r, norbitals_c);
            colIterator += norbitals_c;
        }
        rowIterator += norbitals_r;
    }

    return extendedMFT;
}

/** 
 * Real space implementation of interaction term, valid for both direct and exchange.
 * To compute the direct term, the expected order is (ck,v'k',c'k',vk).
 * For the exchange term, the order is (ck,v'k',vk,c'k').
 * @param coefsK1 First eigenstate vector.
 * @param coefsK2 Second eigenstate vector.
 * @param coefsK3 Third eigenstate vector.
 * @param coefsK4 Fourth eigenstate vector.
 * @param motifFT Motif Fourier transform.
 * @return Interaction term.
 */
std::complex<double> exactInteractionTermMFT(const arma::cx_vec& coefsK1, 
                                     const arma::cx_vec& coefsK2,
                                     const arma::cx_vec& coefsK3, 
                                     const arma::cx_vec& coefsK4,
                                     const arma::cx_mat& motifFT,
                                     const arma::mat& motif,
                                     const arma::urowvec& orbitals){
    
    arma::cx_vec firstCoefArray = arma::conj(coefsK1) % coefsK3;
    arma::cx_vec secondCoefArray = arma::conj(coefsK2) % coefsK4;
    std::complex<double> term = arma::dot(firstCoefArray, extendMotifFT(motifFT, coefsK1.n_elem, motif, orbitals) * secondCoefArray);
        
    return term;
};

}