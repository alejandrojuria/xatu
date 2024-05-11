#include "xatu/Overlap2Centers.hpp"
#include <chrono>

namespace xatu {

/**
 * Constructor that also initializes IntegralsBase from the DFT files. It computes the matrices in the auxiliary basis set.
 * Should only be used to test the 2-center overlap matrices in isolation. 
 * @param GTFConfig GTFConfiguration object.
 * @param o2Mat_name Name of the file where the 2-center overlap matrices will be stored as a cube (o2Mat_name.o2c).
 * @param comp True => Compute the integrals, False => Don't compute the integrals.
 */
Overlap2Centers::Overlap2Centers(const GTFConfiguration& GTFConfig, const uint32_t nR, const std::string& o2Mat_name, const bool comp) : IntegralsBase{GTFConfig} {

    if(comp){
        overlap2Cfun(nR, o2Mat_name);
    }

}

/**
 * Constructor that copies a pre-initialized IntegralsBase.
 * @param IntBase IntegralsBase object.
 * @param o2Mat_name Name of the file where the 2-center overlap matrices will be stored as a cube (o2Mat_name.o2c).
 * @param basis_id True => SCF basis, False => Auxiliary basis.
 */
Overlap2Centers::Overlap2Centers(const IntegralsBase& IntBase, const uint32_t nR, const std::string& o2Mat_name, const bool basis_id) : IntegralsBase{IntBase} {

    overlap2Cfun(nR, o2Mat_name, basis_id);

}

/**
 * Method to compute the overlap matrices in the auxiliary (if basis_id == false) or SCF (if basis_id == true) basis (<P,0|P',R> or <mu,0|mu',R>) 
 * for the first nR Bravais vectors R. These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
 * The resulting cube (third dimension spans the Bravais vectors) is saved in the o2Mat_name.o2c file.
 * The basis_id parameter determines the basis for which the integrals will be computed: true => SCF, false => AUX.
 */
void Overlap2Centers::overlap2Cfun(const uint32_t nR, const std::string& o2Mat_name, const bool basis_id){

arma::mat RlistAU = ANG2AU*generateRlist(bravaisLattice, nR);
uint32_t nR_star = RlistAU.n_cols;
std::map<uint32_t,uint32_t> RlistOpposites = generateRlistOpposite(RlistAU);

uint32_t dimMat {basis_id? dimMat_SCF : dimMat_AUX};
std::vector<std::vector<int>> orbitals_info_int {basis_id? orbitals_info_int_SCF : orbitals_info_int_AUX};
std::vector<std::vector<double>> orbitals_info_real {basis_id? orbitals_info_real_SCF : orbitals_info_real_AUX};
std::vector<double> FAC12 {basis_id? FAC12_SCF : FAC12_AUX};
std::vector<std::vector<double>> FAC3 {basis_id? FAC3_SCF : FAC3_AUX};

std::string basis_string {basis_id? "SCF" : "AUX"};
std::cout << "Computing " << nR_star << " " << dimMat << "x" << dimMat << " 2-center overlap matrices in the " << basis_string  << " basis..." << std::endl;
auto begin = std::chrono::high_resolution_clock::now();  

    uint64_t nelem_triang = 0.5*dimMat*(dimMat + 1);
    uint64_t total_elem = nelem_triang*nR_star;
    arma::cube overlap2Matrices {arma::zeros<arma::cube>(dimMat,dimMat,nR_star)};

    #pragma omp parallel for 
    for(uint64_t s = 0; s < total_elem; s++){ //Spans the lower triangle of all the nR_star matrices <P,0|P',R>
        uint64_t sind {s % nelem_triang};    //Index for the corresponding entry in the overlap matrix, irrespective of the specific R
        uint32_t Rind {s / nelem_triang};         //Position in RlistAU (i.e. 0 for R=0) of the corresponding Bravais vector 
        std::array<uint32_t,2> orb_braket {triangInd_to_rowcol.at(sind)}; 
        // arma::colvec R {RlistAU.col(Rind)};  //Bravais vector (a.u.) corresponding to the "s" matrix element
        uint32_t RindOpp    {RlistOpposites.at(Rind)};  //Position in RlistAU (i.e. 0 for R=0) of the opposite of the corresponding Bravais vector 

        uint32_t orb_bra {orb_braket[0]};   //Orbital number (<dimMat) of the bra corresponding to the index s 
        int L_bra  {orbitals_info_int[orb_bra][2]};
        int m_bra  {orbitals_info_int[orb_bra][3]};
        int nG_bra {orbitals_info_int[orb_bra][4]};
        arma::colvec coords_bra {orbitals_info_real[orb_bra][0], orbitals_info_real[orb_bra][1], orbitals_info_real[orb_bra][2]};  //Position (a.u.) of bra atom
        std::vector<int> g_coefs_bra   {g_coefs.at( L_bra*(L_bra + 1) + m_bra )};

        uint32_t orb_ket {orb_braket[1]};   //Orbital number (<dimMat) of the ket corresponding to the index s. orb_ket <= orb_bra (lower triangle)
        int L_ket  {orbitals_info_int[orb_ket][2]};
        int m_ket  {orbitals_info_int[orb_ket][3]};
        int nG_ket {orbitals_info_int[orb_ket][4]};
        arma::colvec coords_ket {RlistAU.col(Rind) + arma::colvec{orbitals_info_real[orb_ket][0], orbitals_info_real[orb_ket][1], orbitals_info_real[orb_ket][2]} };  //Position (a.u.) of ket atom
        std::vector<int> g_coefs_ket   {g_coefs.at( L_ket*(L_ket + 1) + m_ket )};

        double norm_braket {arma::dot(coords_bra - coords_ket, coords_bra - coords_ket)};
        double FAC12_braket = FAC12[orb_bra]*FAC12[orb_ket];

        for(int gaussC_bra = 0; gaussC_bra < nG_bra; gaussC_bra++){ //Iterate over the contracted Gaussians in the bra orbital
            double exponent_bra {orbitals_info_real[orb_bra][2*gaussC_bra + 3]};
            //double d_bra {orbitals_info_real[orb_bra][2*gaussC_bra + 4]};

            for(int gaussC_ket = 0; gaussC_ket < nG_ket; gaussC_ket++){ //Iterate over the contracted Gaussians in the ket orbital
                double exponent_ket {orbitals_info_real[orb_ket][2*gaussC_ket + 3]};
                //double d_ket {orbitals_info_real[orb_ket][2*gaussC_ket + 4]};

                double p {exponent_bra + exponent_ket};  //Exponent coefficient of the Hermite Gaussian
                arma::colvec P {(exponent_bra*coords_bra + exponent_ket*coords_ket)/p};  //Center of the Hermite Gaussian
                double PAx {P(0) - coords_bra(0)}; 
                double PAy {P(1) - coords_bra(1)}; 
                double PAz {P(2) - coords_bra(2)}; 
                double PBx {P(0) - coords_ket(0)}; 
                double PBy {P(1) - coords_ket(1)}; 
                double PBz {P(2) - coords_ket(2)}; 

                double overlap_g_pre {0.};
                std::vector<int>::iterator g_itr_bra {g_coefs_bra.begin()};
                for(int numg_bra = 0; numg_bra < g_coefs_bra[0]; numg_bra++){ //Iterate over the summands of the corresponding spherical harmonic in the bra orbital
                    int i_bra {*(++g_itr_bra)};
                    int j_bra {*(++g_itr_bra)};
                    int k_bra {*(++g_itr_bra)};
                    int g_bra {*(++g_itr_bra)};
                    int Ei_bra {i_bra*(i_bra + 1)/2};
                    int Ej_bra {j_bra*(j_bra + 1)/2};
                    int Ek_bra {k_bra*(k_bra + 1)/2};
                    
                    std::vector<int>::iterator g_itr_ket {g_coefs_ket.begin()};
                    for(int numg_ket = 0; numg_ket < g_coefs_ket[0]; numg_ket++){ //Iterate over the summands of the corresponding spherical harmonic in the ket orbital
                        int i_ket {*(++g_itr_ket)};
                        int j_ket {*(++g_itr_ket)};
                        int k_ket {*(++g_itr_ket)};
                        int g_ket {*(++g_itr_ket)};

                        double Eii0 {(i_bra >= i_ket)? Efunt0(i_ket + Ei_bra, p, PAx, PBx) : Efunt0(i_bra + i_ket*(i_ket + 1)/2, p, PBx, PAx)};
                        double Ejj0 {(j_bra >= j_ket)? Efunt0(j_ket + Ej_bra, p, PAy, PBy) : Efunt0(j_bra + j_ket*(j_ket + 1)/2, p, PBy, PAy)};
                        double Ekk0 {(k_bra >= k_ket)? Efunt0(k_ket + Ek_bra, p, PAz, PBz) : Efunt0(k_bra + k_ket*(k_ket + 1)/2, p, PBz, PAz)};

                        overlap_g_pre += g_bra*g_ket*Eii0*Ejj0*Ekk0;

                    }
                }
                overlap_g_pre *= FAC3[orb_bra][gaussC_bra]*FAC3[orb_ket][gaussC_ket]*std::pow(p,-1.5)*std::exp(-exponent_bra*exponent_ket*norm_braket/p);
                overlap2Matrices(orb_bra,orb_ket,Rind) += overlap_g_pre;
            }
        }
        overlap2Matrices(orb_bra,orb_ket,Rind) *= FAC12_braket*std::pow(PI,1.5);
        if(orb_bra > orb_ket){
            overlap2Matrices(orb_ket,orb_bra,RindOpp) = overlap2Matrices(orb_bra,orb_ket,Rind);
        }

    }

auto end = std::chrono::high_resolution_clock::now(); 
auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 

overlap2Matrices.save(IntFiles_Dir + o2Mat_name + ".o2c",arma::arma_ascii);  //save matrices to file

std::cout << "Done! Elapsed wall-clock time: " << std::to_string( elapsed.count() * 1e-3 ) << " seconds. Matrices stored in the file: " << 
    IntFiles_Dir + o2Mat_name << ".o2c." << std::endl;

double trace_overlap0 {arma::trace(overlap2Matrices.slice(0))};
if(std::abs(trace_overlap0-dimMat) >= 0.1){
    std::cerr << "WARNING! There is a deviation of " << 100*abs(trace_overlap0-dimMat)/dimMat << 
        "% in the trace of the 2-center overlap matrix in the reference cell. Please, contact the devs!" << std::endl;
}

}

/**
 * Method to compute the E^{i,i'}_{0} coefficients, for i,i'<=4. Returns only the t=0 component.
 */
double Overlap2Centers::Efunt0(const int index, const double p, const double PA, const double PB){

    switch(index)
    {
    case 0:  {// (i,j) = (0,0)
        return 1.0;
    }
    case 1:  {// (i,j) = (1,0)
        return PA;
    }
    case 2:  {// (i,j) = (1,1)
        return (PA*PB + 0.5/p);
    }
    case 3:  {// (i,j) = (2,0)
        return (PA*PA + 0.5/p);
    }
    case 4:  {// (i,j) = (2,1)
        double facp = 0.5/p;
        return (PA*2*facp + PB*(PA*PA + facp)); 
    }
    case 5:  {// (i,j) = (2,2)
        double facp = 0.5/p;
        return (facp*(3*facp + PB*(4*PA + PB)) + PA*PA*(facp + PB*PB));  
    }
    case 6:  {// (i,j) = (3,0)
        return (PA*(PA*PA*p + 1.5)/p);
    }
    case 7:  {// (i,j) = (3,1)
        double facp = 0.5/p;
        return ((PA*p*(6*(PA + PB) + PB*PA*PA*p*4) + 3)*facp*facp);
    }
    case 8:  {// (i,j) = (3,2)
        double facp = 0.5/p;
        return ((PA*p*(PA*PA*(PB*PB*p*4 + 2) + 6*PB*(2*PA + PB)) + 9*PA + 6*PB)*facp*facp);
    }
    case 9:  {// (i,j) = (3,3)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return ((PAPBp*(PAPAp*(PBPBp*8 + 12) + 12*(PBPBp + 3*PAPBp) + 54) + 18*(PAPAp + PBPBp) + 15)*facp*facp*facp);
    }
    case 10: {// (i,j) = (4,0)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p; 
        return ((PAPAp*4*(PAPAp + 3) + 3)*facp*facp);
    }
    case 11: {// (i,j) = (4,1)
        double facp = 0.5/p;
        return ((PA*PA*p*4*((PA*PB*p + 2)*PA + 3*PB) + 12*PA + 3*PB)*facp*facp); 
    }
    case 12: {// (i,j) = (4,2)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PAPAp*(PAPAp*(4 + PB*PB*p*8) + PB*p*(32*PA + 24*PB) + 36) + PB*6*p*(8*PA + PB) + 15)*facp*facp*facp);
    }
    case 13: {// (i,j) = (4,3)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*PB*(8*PBPBp + 12) + 24*PBPBp*(2*PA + PB) + 24*PA + 108*PB) + PBPBp*(72*PA + 6*PB) + 60*PA + 45*PB)*facp*facp*facp);
    }
    case 14: {// (i,j) = (4,4)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return ((PAPAp*(PAPAp*(PBPBp*16*(PBPBp + 3) + 12) + PAPBp*(PBPBp*128 + 192) + PBPBp*48*(PBPBp + 9) + 180)
        + PAPBp*(PBPBp*192 + 480) + PBPBp*12*(PBPBp + 15) + 105)*facp_to2*facp_to2);
    }
    default: {
        throw std::invalid_argument("Overlap2Centers::Efunt0 error: the E^{i,i'}_{0} coefficients are being evaluated for i and/or i' >= 5");
    }
    }
        
}


}