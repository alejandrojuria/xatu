#include "xatu/IntegralsBase.hpp"


namespace xatu {

/**
 * Standard constructor. Should only be used once, after that use the copy constructor.
 * @details Constructor which takes in a GTFConfiguration object, i.e.
 * to init the parameters needed for all integrals from the DFT configuration files.
 * @param configuration GTFConfiguration object.
 * @param integrals_directory Location where the integrals files generated in this family of classes will be stored. The folder must exist in the first place!
 */
IntegralsBase::IntegralsBase(const GTFConfiguration& GTFconfig, const std::string& integrals_directory) {

    this->IntFiles_Dir = integrals_directory;
	initializeBasesAttributes(GTFconfig);
    buildOrbitalsInfo(GTFconfig);

    gfun( maxL );
    triangIndfun(dimMat_AUX);
    FAC12fun(GTFconfig, maxL, true );
    FAC12fun(GTFconfig, maxL, false );
    FAC3fun( true );
    FAC3fun( false );

};

/**
 * Copy constructor. It should be used after the first use of the standard constructor.
 */
IntegralsBase::IntegralsBase(const IntegralsBase& IntBase) = default;

/**
 * Method to initialize IntegralsBase attributes from GTFConfiguration object.
 */
void IntegralsBase::initializeBasesAttributes(const GTFConfiguration& GTFconfig){

    ndim_                    = GTFconfig.ndim;	
    nspecies_                = GTFconfig.nspecies;
    motif_                   = GTFconfig.motif;
    RlistAU_                 = ANG2AU*(GTFconfig.Rlist);   //conversion from Angstrom to atomic units
    RlistOpposites_          = GTFconfig.RlistOpposites;
    R_unpaired_              = GTFconfig.R_unpaired;
    
    natoms_ = motif.n_rows;
	ncells_ = RlistAU.n_cols;

}


/**
 * Method to build the arma::mat orbitals_info attributes. 
 */
void IntegralsBase::buildOrbitalsInfo(const GTFConfiguration& GTFconfig){

    std::vector<int> maxL_spec_SCF, maxL_spec_AUX;
    for(int atom = 0; atom < natoms; atom++){
        int spec = motif.at(atom,3);
        std::vector<int> L_spec_SCF = GTFconfig.L_all_species_SCF[spec];
        std::vector<int> L_spec_AUX = GTFconfig.L_all_species_AUX[spec];
        maxL_spec_SCF.push_back( *std::max_element(L_spec_SCF.begin(),L_spec_SCF.end()) );
        maxL_spec_AUX.push_back( *std::max_element(L_spec_AUX.begin(),L_spec_AUX.end()) );
    }
    int maxL_SCF  = *std::max_element(maxL_spec_SCF.begin(),maxL_spec_SCF.end()); //maximum L q.num. among all species in SCF basis
    int maxL_AUX  = *std::max_element(maxL_spec_AUX.begin(),maxL_spec_AUX.end()); //maximum L q.num. among all species in AUX basis
    int maxL = std::max(maxL_SCF,maxL_AUX);                                       //maximum L q.num. among both basis sets

    std::vector<std::vector<int>>    orbitals_info_int_SCF,  orbitals_info_int_AUX;
    std::vector<std::vector<double>> orbitals_info_real_SCF, orbitals_info_real_AUX;

    for(int atom = 0; atom < natoms; atom++){
        int spec = motif(atom,3);
        std::vector<double> coords_atom {ANG2AU*motif(atom,0), ANG2AU*motif(atom,1), ANG2AU*motif(atom,2)};
        int nsh_spec_SCF = GTFconfig.nshells_all_species_SCF[spec];
        int nsh_spec_AUX = GTFconfig.nshells_all_species_AUX[spec];
        std::vector<int> L_spec_SCF  = GTFconfig.L_all_species_SCF[spec];
        std::vector<int> L_spec_AUX  = GTFconfig.L_all_species_AUX[spec];
        std::vector<int> nG_spec_SCF = GTFconfig.nG_all_species_SCF[spec];
        std::vector<int> nG_spec_AUX = GTFconfig.nG_all_species_AUX[spec];
        std::vector<std::vector<double>> shells_spec_SCF = GTFconfig.shells_all_species_SCF[spec];
        std::vector<std::vector<double>> shells_spec_AUX = GTFconfig.shells_all_species_AUX[spec];

        for(int shl = 0; shl < nsh_spec_SCF; shl++){
            int L_shl  = L_spec_SCF[shl];
            int nG_shl = nG_spec_SCF[shl];
            std::vector<double> contGs_in_shell_SCF      = shells_spec_SCF[shl];
            std::vector<double> orbitals_info_real_SCF_pre = coords_atom;
            orbitals_info_real_SCF_pre.insert(orbitals_info_real_SCF_pre.end(), contGs_in_shell_SCF.begin(), contGs_in_shell_SCF.end());
            for(int m = -L_shl; m <= L_shl; m++){
                orbitals_info_int_SCF.insert(orbitals_info_int_SCF.end(), {spec,shl,L_shl,m,nG_shl});
                orbitals_info_real_SCF.push_back( orbitals_info_real_SCF_pre );
            }
        }

        for(int shl = 0; shl < nsh_spec_AUX; shl++){
            int L_shl  = L_spec_AUX[shl];
            int nG_shl = nG_spec_AUX[shl];
            std::vector<double> contGs_in_shell_AUX      = shells_spec_AUX[shl];
            std::vector<double> orbitals_info_real_AUX_pre = coords_atom;
            orbitals_info_real_AUX_pre.insert(orbitals_info_real_AUX_pre.end(), contGs_in_shell_AUX.begin(), contGs_in_shell_AUX.end());
            for(int m = -L_shl; m <= L_shl; m++){
                orbitals_info_int_AUX.insert(orbitals_info_int_AUX.end(), {spec,shl,L_shl,m,nG_shl});
                orbitals_info_real_AUX.push_back( orbitals_info_real_AUX_pre );
            }
        }
        
    }

    this->orbitals_info_int_SCF  = orbitals_info_int_SCF;
    this->orbitals_info_real_SCF = orbitals_info_real_SCF;
    this->orbitals_info_int_AUX  = orbitals_info_int_AUX;
    this->orbitals_info_real_AUX = orbitals_info_real_AUX;
    this->dimMat_SCF = ( orbitals_info_int_SCF.size() );
    this->dimMat_AUX = ( orbitals_info_int_AUX.size() );
    this->maxL = maxL;

}

/**
* Method to create the matrix of the first nR (at least) Bravais vectors, stored by columns and ordered by ascending norm.
* The number of returned vectors is at least nR because full stars are given. bravaisLattice contains R1,R2,R3 by columns.
* It basically substitutes RlistAU for the integrals when more R-vectors are requested than contained in the .outp.
*/
arma::mat IntegralsBase::generateRlist(const arma::mat& bravaisLattice, const int nR){

    arma::rowvec norms_Ri = arma::sqrt(arma::sum(bravaisLattice % bravaisLattice));
    // Automatic correction accounting for possibly large differences of norms in the lattice vectors
    int normRatio = std::ceil(0.5*arma::max(norms_Ri.cols(0,ndim-1))/arma::min(norms_Ri.cols(0,ndim-1))); 
    // Conservative estimate to make sure that none of the first n vectors is left out
    double RindmaxAux = std::ceil(3*normRatio*std::sqrt(nR));
    arma::rowvec Rindmax;
    if(ndim == 3){
        Rindmax = {RindmaxAux, RindmaxAux, RindmaxAux};
    } else if(ndim == 2){
        Rindmax = {RindmaxAux, RindmaxAux, 0}; 
    } else if(ndim == 1){
        Rindmax = {RindmaxAux, 0, 0}; 
    } else {
        throw std::logic_error("IntegralsBase::generateRlist error: attempted to construct list of Bravais vectors for a finite system");
    }

    arma::mat generated_Rlist = {}; 
    arma::colvec generated_norms = {};     
    for(int Rind1 = -Rindmax(0); Rind1 <= Rindmax(0); Rind1++){
        for(int Rind2 = -Rindmax(1); Rind2 <= Rindmax(1); Rind2++){
            for(int Rind3 = -Rindmax(2); Rind3 <= Rindmax(2); Rind3++){
                arma::irowvec Rindvec = {Rind1, Rind2, Rind3};
                arma::rowvec Rvec = Rindvec.cols(0,ndim-1)*bravaisLattice;
                generated_Rlist = arma::join_vert(generated_Rlist, Rvec);

                arma::colvec norm_aux(1,arma::fill::value(arma::norm(Rvec)));
                generated_norms = arma::join_vert(generated_norms, norm_aux);
            }       
        }
    }

    arma::ucolvec indices = arma::sort_index(generated_norms);
    generated_norms = arma::sort(generated_norms);
    generated_Rlist = generated_Rlist.rows(indices); // Order the lattice vectors (stored as rows still) according to the norms 

    double requested_norm = generated_norms(nR - 1);
    int countr = nR;
    double current_norm = generated_norms(countr);
    std::cout << "Requested direct lattice vectors maximum norm: " << requested_norm << " Angstrom (IntegralsBase::generateRlist)" << std::endl;

    while(current_norm - requested_norm < 1e-3){
        countr++;
        current_norm = generated_norms(countr);
    }
    return (generated_Rlist.rows(0, countr - 1)).t();

}


/**
 * Method to build the matrix of the normalization prefactor FAC1(m,l)->FAC1[l][m]. 
 */
std::vector<std::vector<double>> IntegralsBase::FAC1fun(const int maxL) {

    std::vector<std::vector<double>> FAC1;
    FAC1.reserve( maxL+1 );
    double FAC1l;
    for(int l = 0; l <= maxL; l++){
        std::vector<double> FAC1_pre;
        FAC1l = std::pow(2,l)*std::pow(PI,-1.5)/doubleFactorial(2*l-1);
        for(int m = -l; m <= l; m++){
            FAC1_pre.push_back( std::sqrt( FAC1l*(2-(m==0))*factorial(l - std::abs(m))/factorial(l + std::abs(m)) ) );
        }
        FAC1.push_back( FAC1_pre );
    }
    if(maxL >= 2){
    FAC1[2][0] *= std::sqrt(3);
    FAC1[2][1] *= 3;
    FAC1[2][2] *= std::sqrt(3);
    FAC1[2][3] *= 1.5;
    FAC1[2][4] *= 6;
    if(maxL >= 3){
    FAC1[3][0] *= 3*std::sqrt(10);
    FAC1[3][1] *= 1.5*std::sqrt(10);
    FAC1[3][2] *= 1.5;
    FAC1[3][3] *= 0.5*std::sqrt(15);
    FAC1[3][4] *= 3*std::sqrt(10);
    FAC1[3][5] *= 2.5*std::sqrt(6);
    FAC1[3][6] *= 15;
    if(maxL == 4){
    FAC1[4][0] *= 3*std::sqrt(35);
    FAC1[4][1] *= 15*std::sqrt(7);
    FAC1[4][2] *= 7.5*std::sqrt(2);
    FAC1[4][3] *= 1.25*std::sqrt(2);
    FAC1[4][4] *= 0.5*std::sqrt(5);
    FAC1[4][5] *= 2.5*std::sqrt(7);
    FAC1[4][6] *= 7.5*std::sqrt(14);
    FAC1[4][7] *= 26.25*std::sqrt(2);
    FAC1[4][8] *= 420;
    }
    }
    }
    return FAC1;

}

/**
 * Method to build the vector of the normalization prefactor FAC2(shell,l)->FAC2[shell]. 
 * The shells of only one atom per species are included. 
 */
std::vector<double> IntegralsBase::FAC2fun(const GTFConfiguration& GTFconfig, const bool basis_id){

    std::vector<double> FAC2;
    for(int spec = 0; spec < nspecies; spec++){
        int nsh_spec;
        std::vector<int> L_spec, nG_spec;
        std::vector<std::vector<double>> shells_spec;
        if(basis_id){ //basis_id == true => SCF basis; basis_id == false => auxiliary basis
            nsh_spec    = GTFconfig.nshells_all_species_SCF[spec];
            L_spec      = GTFconfig.L_all_species_SCF[spec];
            nG_spec     = GTFconfig.nG_all_species_SCF[spec];
            shells_spec = GTFconfig.shells_all_species_SCF[spec];
        } else {
            nsh_spec    = GTFconfig.nshells_all_species_AUX[spec];
            L_spec      = GTFconfig.L_all_species_AUX[spec];
            nG_spec     = GTFconfig.nG_all_species_AUX[spec];
            shells_spec = GTFconfig.shells_all_species_AUX[spec];
        }
        for(int shl = 0; shl < nsh_spec; shl++){
            int L_shl  = L_spec[shl];
            int nG_shl = nG_spec[shl];
            std::vector<double> contGs_in_shell = shells_spec[shl];

            double FAC2_elem = 0;
            for(int gaussC1 = 0; gaussC1 < nG_shl; gaussC1++){
                double alpha1 = contGs_in_shell[2*gaussC1];
                double d1     = contGs_in_shell[2*gaussC1 + 1];
                for(int gaussC2 = 0; gaussC2 < nG_shl; gaussC2++){
                    double alpha2 = contGs_in_shell[2*gaussC2];
                    double d2     = contGs_in_shell[2*gaussC2 + 1];
                    FAC2_elem += d1*d2*std::pow( std::sqrt(alpha1*alpha2)/(alpha1+alpha2), L_shl + 1.5 );
                }
            }
            FAC2.push_back( std::pow(FAC2_elem,-0.5) );
   
        }
    }
    return FAC2;

}

/**
 * Method to build the vector FAC12[orb] attribute: FAC12[orb] = FAC1[l(shell)][m(orb)] * FAC2[shell]. 
 */
void IntegralsBase::FAC12fun(const GTFConfiguration& GTFconfig, const int maxL, const bool basis_id){

    std::vector<std::vector<double>> FAC1 = FAC1fun( maxL );
    std::vector<double> FAC2 = FAC2fun(GTFconfig, basis_id);

    std::vector<double> FAC12;
    if(basis_id){ //basis_id == true => SCF basis; basis_id == false => auxiliary basis
        FAC12.reserve( dimMat_SCF );
        for(uint32_t orb = 0; orb < dimMat_SCF; orb++){
            std::vector<int> orbital = orbitals_info_int_SCF[orb]; 
            int L_orb = orbital[2];
            int m_orb = orbital[3];
            int shl_spec_orb = orbital[1];
            for(int spec = 0; spec < orbital[0]; spec++){
                shl_spec_orb += GTFconfig.nshells_all_species_SCF[spec];
            }
            FAC12.push_back( (FAC1[L_orb][m_orb + L_orb])*(FAC2[shl_spec_orb]) );
        }
        this->FAC12_SCF = FAC12;

    } else {
        FAC12.reserve( dimMat_AUX );
        for(uint32_t orb = 0; orb < dimMat_AUX; orb++){
            std::vector<int> orbital = orbitals_info_int_AUX[orb]; 
            int L_orb = orbital[2];
            int m_orb = orbital[3];
            int shl_spec_orb = orbital[1];
            for(int spec = 0; spec < orbital[0]; spec++){
                shl_spec_orb += GTFconfig.nshells_all_species_AUX[spec];
            }
            FAC12.push_back( (FAC1[L_orb][m_orb + L_orb])*(FAC2[shl_spec_orb]) );
        }
        this->FAC12_AUX = FAC12;
    }

}

/**
 * Method to build the vector of vectors FAC3[orb][gaussian] attribute. 
 */
void IntegralsBase::FAC3fun(const bool basis_id){

    std::vector<std::vector<double>> FAC3;
    if(basis_id){ //basis_id == true => SCF basis; basis_id == false => auxiliary basis
        FAC3.reserve( dimMat_SCF );
        for(uint32_t orb = 0; orb < dimMat_SCF; orb++){
            std::vector<double> FAC3_pre;
            std::vector<int> orbital          = orbitals_info_int_SCF[orb]; 
            std::vector<double> orbital_coefs = orbitals_info_real_SCF[orb]; 
            int L_orb  = orbital[2];
            int nG_orb = orbital[4];
            FAC3_pre.reserve( nG_orb );
            for(int nG = 0; nG < nG_orb; nG++){
                FAC3_pre.push_back( orbital_coefs[2*nG+4]*std::pow(orbital_coefs[2*nG+3],0.5*L_orb+0.75) );
            }
            FAC3.push_back( FAC3_pre );
        }
        this->FAC3_SCF = FAC3;

    } else {
        FAC3.reserve( dimMat_AUX );
        for(uint32_t orb = 0; orb < dimMat_AUX; orb++){
            std::vector<double> FAC3_pre;
            std::vector<int> orbital          = orbitals_info_int_AUX[orb]; 
            std::vector<double> orbital_coefs = orbitals_info_real_AUX[orb]; 
            int L_orb  = orbital[2];
            int nG_orb = orbital[4];
            FAC3_pre.reserve( nG_orb );
            for(int nG = 0; nG < nG_orb; nG++){
                FAC3_pre.push_back( orbital_coefs[2*nG+4]*std::pow(orbital_coefs[2*nG+3],0.5*L_orb+0.75) );
            }
            FAC3.push_back( FAC3_pre );
        }
        this->FAC3_AUX = FAC3;
    }

}

/**
 * Method to build the unordered_map method containing the g^{l,m}_{i,j,k} expansion coefficients.
 */
void IntegralsBase::gfun(const int maxL){

    std::unordered_map<int,std::vector<int>> g_coefs;
    // l = 0
    g_coefs[0] = {1,0,0,0,1}; 
    // l = 1
    g_coefs[1] = {1,1,0,0,1};
    g_coefs[2] = {1,0,1,0,1};
    g_coefs[3] = {1,0,0,1,1};
    if(maxL >= 2){
    // l = 2
    g_coefs[4] = {3,2,0,0,-1,0,2,0,-1,0,0,2,2};
    g_coefs[5] = {1,1,0,1,1};
    g_coefs[6] = {1,0,1,1,1};
    g_coefs[7] = {2,2,0,0,1,0,2,0,-1};
    g_coefs[8] = {1,1,1,0,1};
    if(maxL >= 3){
    // l = 3
    g_coefs[9]  = {3,2,0,1,-3,0,2,1,-3,0,0,3,2};
    g_coefs[10] = {3,3,0,0,-1,1,2,0,-1,1,0,2,4};
    g_coefs[11] = {3,2,1,0,-1,0,3,0,-1,0,1,2,4};
    g_coefs[12] = {2,2,0,1,1,0,2,1,-1};
    g_coefs[13] = {1,1,1,1,1};
    g_coefs[14] = {2,3,0,0,1,1,2,0,-3};
    g_coefs[15] = {2,2,1,0,3,0,3,0,-1};
    if(maxL >= 4){
    // l = 4
    g_coefs[16] = {6,4,0,0,3,0,4,0,3,0,0,4,8,2,0,2,-24,0,2,2,-24,2,2,0,6};
    g_coefs[17] = {3,3,0,1,-3,1,2,1,-3,1,0,3,4};
    g_coefs[18] = {3,2,1,1,-3,0,3,1,-3,0,1,3,4};
    g_coefs[19] = {4,4,0,0,-1,0,4,0,1,2,0,2,6,0,2,2,-6};
    g_coefs[20] = {3,3,1,0,-1,1,3,0,-1,1,1,2,6};
    g_coefs[21] = {2,3,0,1,1,1,2,1,-3};
    g_coefs[22] = {2,2,1,1,3,0,3,1,-1};
    g_coefs[23] = {3,4,0,0,1,0,4,0,1,2,2,0,-6};
    g_coefs[24] = {2,3,1,0,1,1,3,0,-1};
    }
    }
    }
    this->g_coefs = g_coefs;

}

/**
 * Method to build the unordered_map method containing the inverse of the (bijective) function: s(i,j) = j + i(i+1)/2.
 */
void IntegralsBase::triangIndfun(const uint32_t dimMat_AUX){
    
    std::unordered_map<uint64_t,std::array<uint32_t,2>> triangInd_to_rowcol;
    uint64_t countr = 0;
    for(uint32_t i = 0; i < dimMat_AUX; i++){
        for(uint32_t j = 0; j <= i; j++){
            triangInd_to_rowcol[countr] = {i,j};
            countr++;
        }
    }
    this->triangInd_to_rowcol = triangInd_to_rowcol;

}

/**
 * Method to compute the E^{i,i'}_{t} coefficients, for i,i'<=4. Returns the vector for 0 <= t <= (i+i').
 */
arma::colvec IntegralsBase::Efun(const int index, const double p, const double PA, const double PB){

    switch(index)
    {
    case 0:  {// (i,j) = (0,0)
        return arma::colvec {1.0};
    }
    case 1:  {// (i,j) = (1,0)
        double facp = 0.5/p;
        return arma::colvec {PA, facp};
    }
    case 2:  {// (i,j) = (1,1)
        double facp = 0.5/p;
        return arma::colvec {PA*PB + facp,  (PA + PB)*facp,  facp*facp};
    }
    case 3:  {// (i,j) = (2,0)
        double facp = 0.5/p;
        return arma::colvec {PA*PA + facp,  PA*2*facp,  facp*facp};
    }
    case 4:  {// (i,j) = (2,1)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        return arma::colvec {PA*2*facp + PB*(PA*PA + facp),  (PA*2*p*(PA + 2*PB) + 3)*facp_to2, (2*PA + PB)*facp_to2,  facp_to2*facp};
    }
    case 5:  {// (i,j) = (2,2)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return arma::colvec {(4*PAPAp*PBPBp + 2*PAPAp + 2*PBPBp + 8*PAPBp + 3)*facp_to2,  (PA + PB)*(PAPBp*4 + 6)*facp_to2,
            (PAPAp + PBPBp + 4*PAPBp + 3)*2*facp_to3,  (PA + PB)*2*facp_to3, facp_to3*facp}; 
    }
    case 6:  {// (i,j) = (3,0)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double PAPAp = PA*PA*p;
        return arma::colvec {PA*(PAPAp*2 + 3)*facp,  (PAPAp*6 + 3)*facp_to2,  3*PA*facp_to2,  facp_to2*facp};
    }
    case 7:  {// (i,j) = (3,1)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        return arma::colvec {(4*PAPAp*PAPBp + 6*PAPAp + 6*PAPBp + 3)*facp_to2,  (PA*(PAPAp*2 + PAPBp*6 + 9) + 3*PB)*facp_to2,
            (PAPAp + PAPBp + 1)*6*facp_to3,  (3*PA + PB)*facp_to3,  facp_to3*facp};
    }
    case 8:  {// (i,j) = (3,2)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        double facp_to4 = facp_to3*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return arma::colvec {(PAPAp*PA*(PBPBp*4 + 2) + PAPBp*6*(2*PA + PB) + 9*PA + 6*PB)*facp_to2,
            (PAPAp*(8*PAPBp + 12*PBPBp + 18) + 36*PAPBp + 6*PBPBp + 15)*facp_to3,
            (PAPAp*(PA + 6*PB) + 3*PB*(PAPBp + 2) + 9*PA)*2*facp_to3,  (3*PAPAp + 6*PAPBp + PBPBp + 5)*2*facp_to4,
            (3*PA + 2*PB)*facp_to4, facp_to4*facp};
    }
    case 9:  {// (i,j) = (3,3)
        double facp = 0.5/p;
        double facp_to3 = facp*facp*facp;
        double facp_to4 = facp_to3*facp;
        double facp_to5 = facp_to4*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        double PBPBmas3PAPB = PBPBp + 3*PAPBp;
        return arma::colvec {(PAPBp*(PAPAp*(PBPBp*8 + 12) + 12*PBPBmas3PAPB + 54) + 18*(PAPAp + PBPBp) + 15)*facp_to3,
        (PAPAp*(PBPBp*12 + 6) + 48*PAPBp + 6*PBPBp + 45)*(PA + PB)*facp_to3, 
        (PAPBp*12*(PAPAp + PBPBmas3PAPB + 9) + 36*(PAPAp + PBPBp) + 45)*facp_to4, 
        (PAPAp + PBPBp + 8*PAPBp + 15)*(PA + PB)*2*facp_to4,  ((PAPAp + PBPBmas3PAPB)*6 + 15)*facp_to5,  
        (PA + PB)*3*facp_to5,  facp_to5*facp};
    }
    case 10: {// (i,j) = (4,0)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        double PAPAp = PA*PA*p;
        return arma::colvec {(PAPAp*4*(PAPAp + 3) + 3)*facp_to2,  (PAPAp*8 + 12)*PA*facp,  (PAPAp*12 + 6)*facp_to3,
            PA*4*facp_to3,  facp_to3*facp};
    }
    case 11: {// (i,j) = (4,1)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        double facp_to4 = facp_to3*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        return arma::colvec {(PAPAp*4*((PAPBp + 2)*PA + 3*PB) + 12*PA + 3*PB)*facp_to2, 
            (PAPAp*4*(PAPAp + PAPBp*4 + 9) + PAPBp*24 + 15)*facp_to3,  (PAPAp*(8*PA + 12*PB) + 24*PA + 6*PB)*facp_to3,
            (12*PAPAp + 8*PAPBp + 10)*facp_to4,  (4*PA + PB)*facp_to4,  facp_to4*facp};
    }
    case 12: {// (i,j) = (4,2)
        double facp = 0.5/p;
        double facp_to3 = facp*facp*facp;
        double facp_to4 = facp_to3*facp;
        double facp_to5 = facp_to4*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return arma::colvec {(PAPAp*(PAPAp*(4 + PBPBp*8) + 32*PAPBp + 24*PBPBp + 36) + 48*PAPBp + 6*PBPBp + 15)*facp_to3,  
            (PAPAp*8*(3*PA + 9*PB + PAPBp*(PA + 2*PB)) + PB*(PAPBp*24 + 30) + 60*PA)*facp_to3,  
            (PAPAp*4*(PAPBp*8 + PBPBp*6 + PAPAp + 18) + 96*PAPBp + 12*PBPBp + 45)*facp_to4,  
            (PAPAp*8*(PA + 3*PB) + PA*8*(PBPBp + 5) + 20*PB)*facp_to4,  (12*PAPAp + 16*PAPBp + 2*PBPBp + 15)*facp_to5,
            (2*PA + PB)*2*facp_to5,  facp_to5*facp};
    }
    case 13: {// (i,j) = (4,3)
        double facp = 0.5/p;
        double facp_to3 = facp*facp*facp;
        double facp_to4 = facp_to3*facp;
        double facp_to5 = facp_to4*facp;
        double facp_to6 = facp_to5*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p; 
        return arma::colvec {(PAPAp*(PAPAp*PB*(8*PBPBp + 12) + 24*PBPBp*(2*PA + PB) + 24*PA + 108*PB) + PBPBp*(72*PA + 6*PB) + 60*PA + 45*PB)*facp_to3,
        (PAPAp*(PAPAp*(12 + PBPBp*24) + PAPBp*(PBPBp*32 + 144) + (PBPBp*216 + 180)) + PAPBp*(PBPBp*48 + 360) + PBPBp*90 + 105)*facp_to4,
        (PAPAp*(12*PA*(PAPBp + PBPBp*4 + 4) + 24*PB*(PBPBp + 9)) + 12*PBPBp*(12*PA + PB) + 180*PA + 135*PB)*facp_to4,
        (PAPAp*4*(PAPAp + 12*PAPBp + 18*PBPBp + 30) + 16*PAPBp*(PBPBp + 15) + 60*PBPBp + 105)*facp_to5,
        (PAPAp*(8*PA + 36*PB) + PBPBp*2*(12*PA + PB) + 60*PA + 45*PB)*facp_to5,  (6*(2*PAPAp + PBPBp + 4*PAPBp) + 21)*facp_to6,
        (4*PA + 3*PB)*facp_to6,  facp_to6*facp};
    }
    case 14: {// (i,j) = (4,4)
        double facp = 0.5/p;
        double facp_to4 = facp*facp*facp*facp;
        double facp_to5 = facp_to4*facp;
        double facp_to6 = facp_to5*facp;
        double facp_to7 = facp_to6*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        double PAmasPB = PA + PB;
        return arma::colvec {(PAPAp*(PAPAp*(PBPBp*16*(PBPBp + 3) + 12) + PAPBp*(PBPBp*128 + 192) + PBPBp*48*(PBPBp + 9) + 180)
        + PAPBp*(PBPBp*192 + 480) + PBPBp*12*(PBPBp + 15) + 105)*facp_to4,  (PAPAp*(PAPBp*(PBPBp*32 + 48) + (PBPBp*240 + 120))
        + PAPBp*(PBPBp*48 + 600) + PBPBp*120 + 420)*PAmasPB*facp_to4,  (PAPAp*(PAPAp*(PBPBp*48 + 24) + PAPBp*128*(PBPBp + 3) 
        + PBPBp*48*(PBPBp + 18) + 540) + PBPBp*(PBPBp*24 + 384*PAPBp + 540) + PAPBp*1440 + 420)*facp_to5,
        (PAPAp*16*(PAPBp + PBPBp*5 + 5) + PAPBp*16*(PBPBp + 25) + PBPBp*80 + 420)*PAmasPB*facp_to5,
        (PAPAp*(4*PAPAp + 64*PAPBp + PBPBp*144 + 180) + PAPBp*(PBPBp*64 + 480) + PBPBp*(4*PBPBp + 180) + 210)*facp_to6,
        (8*(PAPAp + PBPBp + 5*PAPBp) + 84)*PAmasPB*facp_to6,  (12*(PAPAp + PBPBp) + 32*PAPBp + 28)*facp_to7,
        PAmasPB*4*facp_to7,  facp_to7*facp};
    }
    default: {
        throw std::invalid_argument("IntegralsBase::Efun error: the E^{i,i'}_{t} coefficients are being evaluated for i and/or i' >= 5");
    }
    }
        
}

/**
 * Method to compute the factorial of a given integer. 
 */
int IntegralsBase::factorial(int a){
     
    return (a <= 1) ? 1 : a*factorial(a-1);

}

/**
 * Method to compute the double factorial of a given integer. Should not be used for integers <= -2. 
 */
int IntegralsBase::doubleFactorial(int a){
     
    return (a <= 1) ? 1 : a*doubleFactorial(a-2);

}

}
