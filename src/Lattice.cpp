#include "xatu/Lattice.hpp"

namespace xatu {

/**
 * Default constructor to initialize Lattice attributes from a ConfigurationSystem object.
 * @param SystemConfig ConfigurationSystem object. 
 */
Lattice::Lattice(const ConfigurationSystem& SystemConfig){

    this->ndim_   = SystemConfig.ndim;	
	this->Rbasis_ = SystemConfig.Rbasis;
	computeUnitCellVolume();
    calculateGbasis();

}

/**
 * Method to generate a kronecker-like list of integer combinations, to be used with direct or reciprocal lattice vectors.
 * Each row contains the ndim coefficients of a different point. Matrix dimension: (n_{1}*..*n_{ndim}, ndim)
 * @param n1 Number of cells in the direction of the first vector.
 * @param n2 Number of cells in the direction of the second vector. Irrelevant for 1D.
 * @param n3 Number of cells in the direction of the third vector. Irrelevant for 1D and 2D.
 * @param centered If true, the combinations are centered at zero.
 * 				   If false, all combinations have positive coefficients.
 * @returns arma::mat List of cell combinations.
*/
arma::mat Lattice::generateCombinations(const int32_t n1, const int32_t n2, const int32_t n3, const bool centered){

    if(n1 < 0 || n2 < 0 || n3 < 0){
        throw std::invalid_argument("ERROR generateCombinations: number of points per direction cannot be negative");
    }
    arma::icolvec nvec = {1,n1,n2,n3};
    nvec = nvec.subvec(0, ndim);
    nvec.insert_rows(ndim + 1, 1);
    nvec(ndim + 1) = 1;   // nvec = {1,n_{1},..,n_{ndim},1}

    int64_t ncombinations = arma::prod(nvec);
    arma::mat combinations(ncombinations, ndim);

    for(int d = 0; d < ndim; d++){
        int dshift = centered ? (int)nvec(d+1)/2 : 0;
        arma::colvec dvalues = arma::regspace<arma::colvec>(0, nvec(d+1) - 1) - dshift;
        arma::colvec com_aux = arma::repelem( dvalues, arma::prod( nvec.subvec(0,d) ), 1 ); 
        combinations.col(d)  = arma::repmat( com_aux , arma::prod( nvec.subvec(d+2, nvec.n_elem - 1) ), 1 );
    }

    return combinations;

}

arma::mat Lattice::generateCombinations(const int32_t n1, const int32_t n2, const bool centered){

    return generateCombinations(n1,n2,n2,centered);

}

arma::mat Lattice::generateCombinations(const int32_t n1, const bool centered){

    return generateCombinations(n1,n1,n1,centered);

}

/* --------------------------- Reciprocal Lattice methods --------------------------- */

/**
 * Compute the reciprocal lattice vectors {G_1,..,G_ndim} and return them by columns in arma::mat (3,ndim). 
 * The units (Angstrom^-1 or atomic units) are preserved from those of the input Bravais vectors (Angstrom or atomic units).
 * @return void.
 */
void Lattice::calculateGbasis(){
    
    arma::mat Gbasis = arma::zeros<arma::mat>(3,ndim);
    arma::colvec R1 = Rbasis_.col(0);
    if(ndim == 1){
        Gbasis(0,0) = 2*PI/R1(0);
    } 
    else if(ndim == 2){
        arma::colvec R2 = Rbasis_.col(1);
        arma::mat Rot2D = {{0,-1,0},{1,0,0},{0,0,1}};
        arma::mat RotFac = (2*PI/arma::dot(R1,Rot2D*R2))*Rot2D;
        Gbasis.col(0) = RotFac*R2;
        Gbasis.col(1) = -RotFac*R1;
    }
    else if(ndim == 3){
        arma::colvec R2 = Rbasis_.col(1);
        arma::colvec R3 = Rbasis_.col(2);
        double volFac = 2*PI/std::abs(arma::det(Rbasis_));
        Gbasis.col(0) = volFac*arma::cross(R2,R3);
        Gbasis.col(1) = volFac*arma::cross(R3,R1);
        Gbasis.col(2) = volFac*arma::cross(R1,R2);
    }
    else {
        throw std::invalid_argument("ERROR calculateGbasis: lattice dimensionality is not 1, 2 or 3.");
    }

    this->Gbasis_ = Gbasis;

}

/**
 * Compute a Monkhorst-Pack grid in the interval [0 G1)x...x[0 Gn_dim), and return the k-points by columns in arma::mat (3,nk). 
 * The units (Angstrom^-1 or atomic units) are preserved from those of the input basis of reciprocal vector.
 * @param shrink1,shrink2,shrink3 Number of sampled k-points along each reciprocal lattice vectors G1,G2,G3 (resp.). 
 *        Only the first ndim are taken into accout.
 * @param containsGamma True (false) for a grid containing Gamma (displaced by half the corresponding step in each Gi, respectively).
 */
arma::mat Lattice::gridMonkhorstPack(const int32_t shrink1, const int32_t shrink2, const int32_t shrink3, const bool containsGamma){

    arma::mat combs = ( generateCombinations(shrink1, shrink2, shrink3, false) ).t();
    arma::colvec shrink_vec = {(double)shrink1,(double)shrink2,(double)shrink3};
    shrink_vec = shrink_vec.subvec(0, ndim - 1);

    uint32_t nk = combs.n_cols;
    if(!containsGamma){
        combs += 0.5; 
    }
    combs /= arma::repmat(shrink_vec,1,nk);
    arma::mat Klist = Gbasis*combs;

    return Klist;

}

arma::mat Lattice::gridMonkhorstPack(const int32_t shrink1, const int32_t shrink2, const bool containsGamma){

    return gridMonkhorstPack(shrink1, shrink2, shrink2, containsGamma);

}

arma::mat Lattice::gridMonkhorstPack(const int32_t shrink1, const bool containsGamma){

    return gridMonkhorstPack(shrink1, shrink1, shrink1, containsGamma);

}

/* --------------------------- Direct Lattice methods --------------------------- */

/**
 * Method to compute the unit cell volume, area or length (depending on the lattice dimensionality).
 * @return void
 */
void Lattice::computeUnitCellVolume(){

	arma::mat Rbasis_red = Rbasis_.submat(0, 0, ndim - 1, ndim - 1);
	this->unitCellVolume_ = std::abs( arma::det( Rbasis_red ) );

}

/**
* Method to create the matrix of the first nR (at least) 3-component Bravais vectors, stored by columns and ordered by ascending norm.
* @details The number of returned vectors is at least nR because full stars are given. 
* It basically substitutes Rlist for the integrals when more R-vectors are requested than contained in the .outp.
* @param nR Minimum number of direct lattice vectors that will be listed.
* @param IntegralType String which will be printed along with norm(nR) and which indicates the type of integrals for which the 
*        list is generated
* @return arma::mat (3, nR' >= nR) matrix with the aforementioned Bravais vectors by columns.
*/
arma::mat Lattice::generateRlist(const uint32_t nR, const std::string& IntegralType){

    arma::rowvec norms_Ri = arma::sqrt(arma::sum(Rbasis_ % Rbasis_, 0));
    // Automatic correction accounting for possibly large differences of norms in the lattice vectors
    uint normRatio = std::ceil(0.5*arma::max(norms_Ri) / arma::min(norms_Ri)); 
    // Conservative estimate to make sure that none of the first n vectors is left out
    uint32_t RindmaxAux = std::ceil(3*normRatio*std::sqrt(nR));
	RindmaxAux += 1 - (RindmaxAux % 2);
    
	arma::mat combinations = generateCombinations(RindmaxAux, true);
	uint64_t ncombinations = combinations.n_rows;
	arma::mat generated_Rlist(3,ncombinations); 
    arma::rowvec generated_norms(ncombinations);

	#pragma omp parallel for 
	for(uint64_t ind_comb = 0; ind_comb < ncombinations; ind_comb++){
		arma::colvec comb = (combinations.row(ind_comb)).t();
		arma::colvec Rvec = Rbasis_*comb;

		generated_Rlist.col(ind_comb) = Rvec;
		generated_norms(ind_comb) = arma::norm(Rvec);
	}

    arma::urowvec indices = (arma::sort_index(generated_norms)).t();
    generated_norms = arma::sort(generated_norms);
    generated_Rlist = generated_Rlist.cols(indices); // Order the lattice vectors (columns of generated_Rlist) according to the norms 
    
    double requested_norm = generated_norms(nR - 1);
    uint32_t countr = nR;
    double current_norm = generated_norms(countr);
    std::cout << "Requested direct lattice vectors maximum norm: " << requested_norm << " Angstrom (" + IntegralType + ")" << std::endl;

    while(current_norm - requested_norm < 1e-3){ // Complete the current star of direct lattice vectors
        countr++;
        current_norm = generated_norms(countr);
    }
    return generated_Rlist.cols(0, countr - 1);

}

/**
 * Returns a map where each entry is the index of the direct lattice vector in the input generated_Rlist (generalization of Rlist for
 * an arbitrary number of direct lattice vectors) opposite to the lattice vector whose index is the corresponding map's key. 
 * @param generated_Rlist List of Bravais vectors (3,nR), as given by generateRlist. Not necessarily the attribute Rlist
 * @return std::map The list whose n-th entry is -R_{n}, where R_{n} = Rlist(n).
 */
std::map<uint32_t,uint32_t> Lattice::generateRlistOpposite(const arma::mat& generated_Rlist){

    uint32_t nRlist = generated_Rlist.n_cols; 
    std::map<uint32_t,uint32_t> RlistOpposites;
    for(uint32_t RindOpp = 0; RindOpp < nRlist; RindOpp++){
        for(uint32_t Rind = 0; Rind < nRlist; Rind++){
            arma::colvec Rsum = generated_Rlist.col(Rind) + generated_Rlist.col(RindOpp);
            if( Rsum.is_zero(0.001) ){
                RlistOpposites[RindOpp] = Rind;
            }
        }

    }
    return RlistOpposites;

}


}

