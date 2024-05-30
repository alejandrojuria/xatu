#include "xatu/IntegralsOverlap2MixC.hpp"

namespace xatu {

/**
 * Constructor that copies a pre-initialized IntegralsBase.
 * @param IntBase IntegralsBase object.
 * @param intName Name of the file where the mixed 2-center overlap matrices will be stored as a cube (o2MixMat_intName.o2mc).
 */
IntegralsOverlap2MixC::IntegralsOverlap2MixC(const IntegralsBase& IntBase, const uint32_t nR, const std::string& intName) : IntegralsBase{IntBase} {

    overlap2MixCfun(nR, intName);
    
}

/**
 * Method to compute the overlap matrices in the mixed SCF and auxiliary basis sets (<P,0|mu,R>) for the first nR Bravais vectors R. 
 * These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
 * The resulting cube (third dimension spans the Bravais vectors) is saved in the o2MixMat_intName.o2mc.
 * @param nR Minimum number of direct lattice vectors for which the 2-center overlap integrals will be computed.
 * @param intName Name of the file where the mixed 2-center overlap matrices will be stored as a cube (o2MixMat_name.o2mc).
 * @return void. Matrices and the corresponding list of lattice vectors are stored instead.
 */
void IntegralsOverlap2MixC::overlap2MixCfun(const uint32_t nR, const std::string& intName){

const double PIpow = std::pow(PI,1.5);
arma::mat RlistAU = ANG2AU*generateRlist(nR, "Overlap2MixC");
uint32_t nR_star = RlistAU.n_cols;

std::cout << "Computing " << nR_star << " " << dimMat_AUX << "x" << dimMat_SCF << " 2-center mixed overlap matrices..." << std::flush;

// Start the calculation
auto begin = std::chrono::high_resolution_clock::now();  

    uint64_t nelem_rectang = dimMat_AUX*dimMat_SCF;
    uint64_t total_elem = nelem_rectang*nR_star;
    arma::cube overlap2MixMatrices {arma::zeros<arma::cube>(dimMat_AUX,dimMat_SCF,nR_star)};

    #pragma omp parallel for
    for(uint64_t s = 0; s < total_elem; s++){ //Spans all the elements of all the nR_star matrices <P,0|mu,R>
        uint64_t sind {s % nelem_rectang};   //Index for the corresponding entry in the overlap matrix, irrespective of the specific R
        uint32_t Rind {static_cast<uint32_t>(s / nelem_rectang)};   //Position in RlistAU (i.e. 0 for R=0) of the corresponding Bravais vector 
        uint32_t orb_bra {static_cast<uint32_t>(sind % dimMat_AUX)};  //Orbital number (<dimMat_AUX) of the bra corresponding to the index s
        uint32_t orb_ket {static_cast<uint32_t>(sind / dimMat_AUX)};  //Orbital number (<dimMat_SCF) of the ket corresponding to the index s
        // arma::colvec R {RlistAU.col(Rind)};  //Bravais vector (a.u.) corresponding to the "s" matrix element

        int L_bra  {orbitals_info_int_AUX[orb_bra][2]};
        int m_bra  {orbitals_info_int_AUX[orb_bra][3]};
        int nG_bra {orbitals_info_int_AUX[orb_bra][4]};
        arma::colvec coords_bra {orbitals_info_real_AUX[orb_bra][0], orbitals_info_real_AUX[orb_bra][1], orbitals_info_real_AUX[orb_bra][2]};  //Position (a.u.) of bra atom
        std::vector<int> g_coefs_bra   {g_coefs.at( L_bra*(L_bra + 1) + m_bra )};

        int L_ket  {orbitals_info_int_SCF[orb_ket][2]};
        int m_ket  {orbitals_info_int_SCF[orb_ket][3]};
        int nG_ket {orbitals_info_int_SCF[orb_ket][4]};
        arma::colvec coords_ket {RlistAU.col(Rind) + arma::colvec{orbitals_info_real_SCF[orb_ket][0], orbitals_info_real_SCF[orb_ket][1], orbitals_info_real_SCF[orb_ket][2]} };  //Position (a.u.) of ket atom
        std::vector<int> g_coefs_ket   {g_coefs.at( L_ket*(L_ket + 1) + m_ket )};

        arma::colvec coords_braket {coords_bra - coords_ket};
        double norm_braket {arma::dot(coords_braket,coords_braket)};
        double FAC12_braket = FAC12_AUX[orb_bra]*FAC12_SCF[orb_ket];

        for(int gaussC_bra = 0; gaussC_bra < nG_bra; gaussC_bra++){ //Iterate over the contracted Gaussians in the bra orbital
            double exponent_bra {orbitals_info_real_AUX[orb_bra][2*gaussC_bra + 3]};
            //double d_bra {orbitals_info_real_AUX[orb_bra][2*gaussC_bra + 4]};

            for(int gaussC_ket = 0; gaussC_ket < nG_ket; gaussC_ket++){ //Iterate over the contracted Gaussians in the ket orbital
                double exponent_ket {orbitals_info_real_SCF[orb_ket][2*gaussC_ket + 3]};
                //double d_ket {orbitals_info_real_SCF[orb_ket][2*gaussC_ket + 4]};

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
                overlap_g_pre *= FAC3_AUX[orb_bra][gaussC_bra]*FAC3_SCF[orb_ket][gaussC_ket]*std::pow(p,-1.5)*std::exp(-exponent_bra*exponent_ket*norm_braket/p);
                overlap2MixMatrices(orb_bra,orb_ket,Rind) += overlap_g_pre;
            }
        }
        overlap2MixMatrices(orb_bra,orb_ket,Rind) *= FAC12_braket*PIpow;

    }

auto end = std::chrono::high_resolution_clock::now(); 
auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 

// Store the matrices and the list of direct lattice vectors
overlap2MixMatrices.save(IntFiles_Dir + "o2MixMat_" + intName + ".o2mc", arma::arma_ascii); 
RlistAU.save(IntFiles_Dir + "RlistAU_" + intName + ".o2mc", arma::arma_ascii);

std::cout << "Done! Elapsed wall-clock time: " << std::to_string( elapsed.count() * 1e-3 ) << " seconds." << std::endl;
std::cout << "Matrices stored in the file: " << IntFiles_Dir + "o2MixMat_" + intName + ".o2mc" << 
    " , and list of Bravais vectors in " << IntFiles_Dir + "RlistAU_" + intName + ".o2mc" << std::endl;

}

}