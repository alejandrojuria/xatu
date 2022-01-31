/* This file includes several routines to implement the computation 
of the Z2 invariant of a time-reversal topological insulator in an 
approximate way using the spin Bott index. */

#include <armadillo>
#include <complex>
#include <math.h>
#include <stdlib.h>

#include "zigzag.hpp"
#include "libinvariant_dense.hpp"

using namespace arma;
using namespace std::chrono;

/* Routine to calculate the matrix representation of the ground-state 
projected spin P_z=PS_zP. Here, we compute the representation using
the tight-binding eigenfunctions, i.e. Bloch states |n,k>.
Input: Eigenvector coefficient matrix from diagonalization of the TB 
hamiltonian (matrix dimTB x dimTB)
Output: complex matrix PspinMatrix dimTB x dimTB */
cx_mat PspinMatrix(const cx_mat& eigvectors){
	int dimTB = (int)eigvectors.n_cols;
    int numAtoms = dimTB/8;

	cx_mat spinMat = arma::zeros<cx_mat>(dimTB, dimTB);
	cx_vec spinEigvalues = {1./2, -1./2, 1./2, 1./2, 1./2, -1./2, -1./2, -1./2};
	cx_vec spinVector = arma::kron(arma::ones(dimTB/8, 1), spinEigvalues);

	for(int i = 0; i < 5*numAtoms; i++){
		for(int j = i; j < 5*numAtoms; j++){

			cx_vec effectiveCoefs = eigvectors.col(j) % spinVector;
			
			if(i == j){
				spinMat(i,j) = arma::cdot(eigvectors.col(i), effectiveCoefs)/2.;
			}
			else{
			    spinMat(i,j) = arma::cdot(eigvectors.col(i), effectiveCoefs);
			};

		}
	}
	spinMat = spinMat + spinMat.t();

	return spinMat;
};

/* Routine to compute all the necessary matrices for the computation
of the Bott index. Here we calculate the tight-binding coefficient matrix 
(eigenvalues) for ALL n,k (not just n), and for all unit cells 
(not only R=0). We do the same for the spin sector projectors P_+ and P_-.
Note that for those two, their matrix representation is not square.
This is due to the fact that we are considering only occupied states, 
and of top of that we separate them in two sectors. 
Input: vec kpoints, int N (finite direction), int nOrb (num. orbitals). 
Output: void. Actually overwrite predefined matrices. */

void calculateFullTBMatrices(const vec& kpoints, int N, int nOrb){

    int nk = (int)kpoints.n_elem;
    int dimTB = 2*(N+1)*nOrb;
    int dimFull = dimTB*nk;
    int occStates = dimTB*5/8;    

    fullTB_matrix = arma::zeros<cx_mat>(dimFull, dimFull);
    PpTB_matrix = arma::zeros<cx_mat>(dimFull, dimFull);
    PmTB_matrix = arma::zeros<cx_mat>(dimFull, dimFull);

    vec eigval;
    cx_mat eigvec;
    std::complex<double> imag(0,1);
    int it = 0;

    for(int i = 0; i < nk; i++){

		cx_mat h = hamiltonian(kpoints(i), H0, Ha, Hsoc);
		arma::eig_sym(eigval, eigvec, h);

        // Full TB coefficient matrix calculation
        for(int ncell = 0; ncell < nk; ncell++){
            std::complex<double> expKR = exp(imag*kpoints(i)*a*(double)ncell);

            cx_mat coefs = eigvec * expKR;

            fullTB_matrix.submat(ncell*dimTB, i*dimTB, (ncell+1)*dimTB - 1, (i + 1)*dimTB - 1) 
                      = coefs;
        };

        // Pz sectors TB matrix calculation
        cx_mat PzMatrix = PspinMatrix(eigvec);
        arma::eig_sym(eigval, eigvec, PzMatrix);
        cx_vec coefVec;

        for(int j = 0; j < dimTB; j++){

            coefVec = arma::zeros<cx_vec>(dimFull, 1);
            coefVec.subvec(i*dimTB, (i + 1)*dimTB - 1) = eigvec.col(j);

            if(eigval(j) > 0.001){
                PpTB_matrix.col(it) = coefVec;
                it++;
            }
            else if (eigval(j) < -0.001){
                PmTB_matrix.col(it) = coefVec;
                it++;
            };
        };
    };

    return;
};

/* Routine to calculate the matrix representation of exp(i2*PI*X) in the
localized orbital basis. 
Input: int posVar = {0,1}. Denotes the computation either for x (0) or y (1).
mat motiv contains all the positions of the atoms inside the unit cell.
int nOrbitals, int Ncell.
Output: cx_mat positionMatrix, corresponding to X or Y, in the localized 
orbital basis. 
Tha basis ordering is: we fix unit cell and atom, and run over orbitals. 
When finished, change atom inside unit cell (motif) and repeat, until we 
ran over the whole motif. Then switch to the next unit cell and repeat. */
cx_mat positionMatrix(int posVar, const mat& motif, 
                      int nOrbitals, int Ncell, int N){

    int dim = (int)motif.n_rows*(2*Ncell+1)*nOrbitals;
    int dimTB = motif.n_rows*nOrbitals;
    cx_mat positionMatrix = zeros<cx_mat>(dim, dim);
    std::complex<double> imagNum(0,1);
    double Lx = motif.row(motif.n_rows - 1)(0);
    double Ly = motif.row(1)(1) + 2*Ncell*a;

    for(int i = 0; i < dim; i++){

        int actualCell = (int)i/(motif.n_rows*nOrbitals);
        int actualAtom = (int)(i/nOrbitals) - actualCell*motif.n_rows;

        double posValue = motif.row(actualAtom)(posVar) 
                        + posVar*actualCell*a;

        posValue /= (posVar*Ly + abs((posVar - 1))*Lx);

        std::complex<double> expPos = std::exp(imagNum*2.*PI*posValue);
        
        positionMatrix(i,i) = expPos;
    };

    return positionMatrix;
};

/* Routine to compute the spin sector projection P+- matrix in the localized 
orbital representation. -- See theory for details of implementation --.
Input: int sign (0 for +, 1 for -). sp_cx_mat tbStack: sparse matrix
with all the eigenvectors for all n,k
sp_cx_mat Psector: sparse matrix with the eigenvectors of the corresponding
sector (+ or -) for all n,k.
Output: sp_cx_mat projectorMatrix (sector projector representation in 
localized orbital basis) */
cx_mat projectorMatrix(int sign, const arma::cx_mat& tbStack,
                       const arma::cx_mat& Psector, int Ncell){
    
    cx_mat Pprod = Psector*Psector.t();
    cx_mat projectorMatrix = tbStack*Pprod*tbStack.t()/(2.*Ncell + 1);

    return projectorMatrix;
};

double sectorBottIndex(const cx_mat& eX, const cx_mat& eY, 
                       const cx_mat& Psector){

    cx_mat U = Psector*eX*Psector;
    cout << "U computed" << endl;
    cx_mat V = Psector*eY*Psector;
    cout << "V computed" << endl;

    cx_mat totalMat = V*U*V.t()*U.t();
    cout << "VUV.tU.t computed" << endl;

    cx_vec eigval;
    arma::cx_mat eigvec;
    eig_gen(eigval, eigvec, arma::cx_mat(totalMat));
    cout << eigval.subvec(0, 10) << endl;
    std::complex<double> logeig = arma::sum(arma::log(eigval));

    double bottIndex = imag(logeig)/(2*PI);

    return bottIndex;
};