#include "xatu/Overlap3Centers.hpp"

namespace xatu {

/**
 * Constructor that also initializes IntegralsBase from the DFT files. 
 * Should only be used to test the 3-center overlap matrices in isolation.
 * @param GTFConfig GTFConfiguration object.
 * @param tol Tolerance for retaining the entries of the 3-center overlap matrices. These must be > 10^-tol, in absolute value.
 * @param o3Mat_name Name of the file where the 3-center overlap matrices will be stored as a cube (o3Mat_name.o3c).
 * @param o2Mat_name Name of the file where the 2-center overlap matrices in the AUX basis will be stored as a cube (o2Mat_name.o2c).
 */
Overlap3Centers::Overlap3Centers(const GTFConfiguration& GTFConfig, const int tol, const uint32_t nR2, const std::string& o3Mat_name, const std::string& o2Mat_name) 
    : Overlap2Centers{GTFConfig, nR2, o2Mat_name, false}, IntegralsBase{GTFConfig} {

    overlap3Cfun(nR2, tol, o3Mat_name);

}

/**
 * Constructor that copies a pre-initialized IntegralsBase and also computes the 2-center overlap in the AUX basis.
 * @param IntBase IntegralsBase object.
 * @param tol Tolerance for retaining the entries of the 3-center overlap matrices. These must be > 10^-tol, in absolute value.
 * @param nR2 Number of R and R' that will be considered for the integrals. By default it's IntegralsBase.ncells, and cannot be larger than it.
 * @param o3Mat_name Name of the file where the 3-center overlap matrices will be stored as a cube (o3Mat_name.o3c).
 * @param o2Mat_name Name of the file where the 2-center overlap matrices in the AUX basis will be stored as a cube (o2Mat_name.o2c).
 */
Overlap3Centers::Overlap3Centers(const IntegralsBase& IntBase, const int tol, const uint32_t nR2, const std::string& o3Mat_name, const std::string& o2Mat_name) 
    : Overlap2Centers{IntBase, nR2, o2Mat_name, false}, IntegralsBase{IntBase} {

    overlap3Cfun(nR2, tol, o3Mat_name);

}

/**
 * Method to compute the rectangular overlap matrices <P,0|mu,R;mu',R'> in the mixed SCF and auxiliary basis sets for the first nR Bravais
 * vectors R and R' (nR2^2 pairs of vectors). These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
 * The resulting cube (third dimension spans auxiliary basis {P}, while the Bravais lattice vectors R and R' vary by rows and columns
 * (respectively) once the block for {mu,mu'} has been completed, i.e., each slice is kron({mu,mu'},{R,R'}) ) is saved in the o3Mat_name.o2c file.
 */
void Overlap3Centers::overlap3Cfun(const uint32_t nR2, const int tol, const std::string& o3Mat_name){

// openmp directive to define std::vector insert as a reduction. SEE https://stackoverflow.com/questions/18669296/c-openmp-parallel-for-loop-alternatives-to-stdvector
#pragma omp declare reduction (merge : std::vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction (merge2 : std::vector<std::array<uint32_t,5>> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

arma::mat RlistAU = ANG2AU*generateRlist(bravaisLattice, nR2);
uint32_t nR2_star = RlistAU.n_cols;

double etol = std::pow(10.,-tol);
uint64_t dim_Slice = dimMat_SCF*nR2_star;
uint64_t nelems_slice = dim_Slice*dim_Slice;
uint64_t total_elem = nelems_slice*dimMat_AUX;
std::cout << "Computing " << dimMat_AUX << " " << dim_Slice << "x" << dim_Slice << " 3-center overlap matrices..." << std::endl;
auto begin = std::chrono::high_resolution_clock::now();  

    std::vector<double> overlap3Matrices_val; //Each entry is a 3-center overlap matrix element above the defined tolerance 
    std::vector<std::array<uint32_t,5>> overlap3Matrices_ind; //The first dimension is in bijection with the entries of overlap3Matrices_val. The second dimension contains the corresponding P,mu,mu',R,R'

    #pragma omp parallel for reduction(merge: overlap3Matrices_val) reduction(merge2: overlap3Matrices_ind)
    for(uint64_t s = 0; s < total_elem; s++){ //Spans all the matrix elements of the cube
        uint64_t sind {s % nelems_slice};   //Index for the entry within the slice matrix (single index for {mu,R,mu',R'}) corresponding to loop index s. Iterates first over columns (mu,R)
        
        uint32_t Pind {s / nelems_slice};     //Slice of the cube (or AUX basis component P, <dimMat_AUX) corresponding to loop index s 
        uint64_t muRind1 {sind % dim_Slice};  //Row in the slice matrix (single index for {mu,R})
        uint32_t muind1 {muRind1 % dimMat_SCF};    //Orbital number mu (<dimMat_SCF) corresponding to loop index s
        uint32_t Rind1  {muRind1 / dimMat_SCF};    //Direct lattice vector index (<nR2_star) corresponding to loop index s
        uint64_t muRind2 {sind / dim_Slice};  //Column in the slice matrix (single index for {mu',R'}) 
        uint32_t muind2 {muRind2 % dimMat_SCF};    //Orbital number mu' (<dimMat_SCF) corresponding to loop index s 
        uint32_t Rind2  {muRind2 / dimMat_SCF};    //Direct lattice vector R' (<nR2_star) corresponding to loop index s

        int L_P  {orbitals_info_int_AUX[Pind][2]};
        int m_P  {orbitals_info_int_AUX[Pind][3]};
        int nG_P {orbitals_info_int_AUX[Pind][4]};
        arma::colvec coords_P {orbitals_info_real_AUX[Pind][0], orbitals_info_real_AUX[Pind][1], orbitals_info_real_AUX[Pind][2]};  //Position (a.u.) of P atom
        std::vector<int> g_coefs_P {g_coefs.at( L_P*(L_P + 1) + m_P )};

        int L_mu1  {orbitals_info_int_SCF[muind1][2]};
        int m_mu1  {orbitals_info_int_SCF[muind1][3]};
        int nG_mu1 {orbitals_info_int_SCF[muind1][4]};
        // arma::colvec R1 {RlistAU.col(Rind1)};  //Direct lattice vector R corresponding to loop index s
        arma::colvec coords_mu1 {RlistAU.col(Rind1) + arma::colvec{orbitals_info_real_SCF[muind1][0], orbitals_info_real_SCF[muind1][1], orbitals_info_real_SCF[muind1][2]} };  //Position (a.u.) of mu atom
        std::vector<int> g_coefs_mu1 {g_coefs.at( L_mu1*(L_mu1 + 1) + m_mu1 )};

        int L_mu2  {orbitals_info_int_SCF[muind2][2]};
        int m_mu2  {orbitals_info_int_SCF[muind2][3]};
        int nG_mu2 {orbitals_info_int_SCF[muind2][4]};
        // arma::colvec R2 {RlistAU.col(Rind2)};  //Direct lattice vector R corresponding to loop index s
        arma::colvec coords_mu2 {RlistAU.col(Rind2) + arma::colvec{orbitals_info_real_SCF[muind2][0], orbitals_info_real_SCF[muind2][1], orbitals_info_real_SCF[muind2][2]} };  //Position (a.u.) of mu' atom
        std::vector<int> g_coefs_mu2 {g_coefs.at( L_mu2*(L_mu2 + 1) + m_mu2 )};

        double norm_mu1mu2 {arma::dot(coords_mu1 - coords_mu2,coords_mu1 - coords_mu2)};
        double overlap3_g_pre0 {0.};
        for(int gaussC_mu1 = 0; gaussC_mu1 < nG_mu1; gaussC_mu1++){ //Iterate over the contracted Gaussians in the mu orbital 
            double exponent_mu1 {orbitals_info_real_SCF[muind1][2*gaussC_mu1 + 3]};

            for(int gaussC_mu2 = 0; gaussC_mu2 < nG_mu2; gaussC_mu2++){ //Iterate over the contracted Gaussians in the mu' orbital
                double exponent_mu2 {orbitals_info_real_SCF[muind2][2*gaussC_mu2 + 3]};
                double p {exponent_mu1 + exponent_mu2};  //Exponent coefficient of the Hermite Gaussian between mu,mu'
                arma::colvec P {(exponent_mu1*coords_mu1 + exponent_mu2*coords_mu2)/p};  //Center of the Hermite Gaussian between mu,mu'
                double PAx {P(0) - coords_mu1(0)}; 
                double PAy {P(1) - coords_mu1(1)}; 
                double PAz {P(2) - coords_mu1(2)}; 
                double PBx {P(0) - coords_mu2(0)}; 
                double PBy {P(1) - coords_mu2(1)}; 
                double PBz {P(2) - coords_mu2(2)};

                double norm_PP_3 {arma::dot(P - coords_P, P - coords_P)};
                for(int gaussC_P = 0; gaussC_P < nG_P; gaussC_P++){ //Iterate over the contracted Gaussians in the P orbital
                    double exponent_P {orbitals_info_real_AUX[Pind][2*gaussC_P + 3]};

                    double p_3 {exponent_P + p};  //Exponent coefficient of the triple Hermite Gaussian (between P and resulting of the Hermite Gaussian mu,mu')
                    arma::colvec P_3 {(p*P + exponent_P*coords_P)/p_3};  //Center of the triple Hermite Gaussian
                    double PAx_3 {P_3(0) - P(0)};
                    double PAy_3 {P_3(1) - P(1)};
                    double PAz_3 {P_3(2) - P(2)};
                    double PBx_3 {P_3(0) - coords_P(0)};
                    double PBy_3 {P_3(1) - coords_P(1)};
                    double PBz_3 {P_3(2) - coords_P(2)};

                    double overlap3_g_pre1 {0.};
                    std::vector<int>::iterator g_itr_mu1 {g_coefs_mu1.begin()};
                    for(int numg_mu1 = 0; numg_mu1 < g_coefs_mu1[0]; numg_mu1++){ //Iterate over the summands of the corresponding spherical harmonic in the mu orbital
                        int i_mu1 {*(++g_itr_mu1)};
                        int j_mu1 {*(++g_itr_mu1)};
                        int k_mu1 {*(++g_itr_mu1)};
                        int g_mu1 {*(++g_itr_mu1)};
                        int Ei_mu1 {i_mu1*(i_mu1 + 1)/2};
                        int Ej_mu1 {j_mu1*(j_mu1 + 1)/2};
                        int Ek_mu1 {k_mu1*(k_mu1 + 1)/2};

                        std::vector<int>::iterator g_itr_mu2 {g_coefs_mu2.begin()};
                        for(int numg_mu2 = 0; numg_mu2 < g_coefs_mu2[0]; numg_mu2++){ //Iterate over the summands of the corresponding spherical harmonic in the mu' orbital
                            int i_mu2 {*(++g_itr_mu2)};
                            int j_mu2 {*(++g_itr_mu2)};
                            int k_mu2 {*(++g_itr_mu2)};
                            int g_mu2 {*(++g_itr_mu2)};
                            int Ei_mu2 {i_mu2*(i_mu2 + 1)/2};
                            int Ej_mu2 {j_mu2*(j_mu2 + 1)/2};
                            int Ek_mu2 {k_mu2*(k_mu2 + 1)/2};

                            arma::colvec Eii_vec {(i_mu1 >= i_mu2)? Efun(i_mu2 + Ei_mu1, p, PAx, PBx) : Efun(i_mu1 + Ei_mu2, p, PBx, PAx)};
                            arma::colvec Ejj_vec {(j_mu1 >= j_mu2)? Efun(j_mu2 + Ej_mu1, p, PAy, PBy) : Efun(j_mu1 + Ej_mu2, p, PBy, PAy)};
                            arma::colvec Ekk_vec {(k_mu1 >= k_mu2)? Efun(k_mu2 + Ek_mu1, p, PAz, PBz) : Efun(k_mu1 + Ek_mu2, p, PBz, PAz)};

                            std::vector<int>::iterator g_itr_P {g_coefs_P.begin()};
                            for(int numg_P = 0; numg_P < g_coefs_P[0]; numg_P++){ //Iterate over the summands of the corresponding spherical harmonic in the P orbital
                                int i_P {*(++g_itr_P)};
                                int j_P {*(++g_itr_P)};
                                int k_P {*(++g_itr_P)};
                                int g_P {*(++g_itr_P)};

                                // X contribution
                                double sum_t {0.};
                                for(int t = 0; t <= (i_mu1 + i_mu2); t++){
                                    std::vector<double> Dt_vec {Dfun(t, p)};
                                    double sum_ipp {0.};
                                    double E3ii_0;
                                    for(int i_pp = t; i_pp >= 0; i_pp -= 2){
                                        if(i_pp <= 4){
                                            E3ii_0 = (i_pp >= i_P)? Efunt0(i_P + i_pp*(i_pp+1)/2, p_3, PAx_3, PBx_3) : Efunt0(i_pp + i_P*(i_P+1)/2, p_3, PBx_3, PAx_3);
                                        } else {
                                            E3ii_0 = (i_pp >= i_P)? EfunTriplet0(5*i_pp + i_P - 25, p_3, PAx_3, PBx_3) : EfunTriplet0(5*i_P + i_pp - 25, p_3, PBx_3, PAx_3);
                                        }
                                        sum_ipp += Dt_vec[i_pp/2]*E3ii_0;
                                    }
                                    sum_t += Eii_vec(t)*sum_ipp;
                                }   
                                // Y contribution
                                double sum_u {0.};
                                for(int u = 0; u <= (j_mu1 + j_mu2); u++){
                                    std::vector<double> Du_vec {Dfun(u, p)};
                                    double sum_jpp {0.};
                                    double E3jj_0;
                                    for(int j_pp = u; j_pp >= 0; j_pp -= 2){
                                        if(j_pp <= 4){
                                            E3jj_0 = (j_pp >= j_P)? Efunt0(j_P + j_pp*(j_pp+1)/2, p_3, PAy_3, PBy_3) : Efunt0(j_pp + j_P*(j_P+1)/2, p_3, PBy_3, PAy_3);
                                        } else {
                                            E3jj_0 = (j_pp >= j_P)? EfunTriplet0(5*j_pp + j_P - 25, p_3, PAy_3, PBy_3) : EfunTriplet0(5*j_P + j_pp - 25, p_3, PBy_3, PAy_3);
                                        }
                                        sum_jpp += Du_vec[j_pp/2]*E3jj_0;
                                    }
                                    sum_u += Ejj_vec(u)*sum_jpp;
                                }
                                // Z contribution
                                double sum_v {0.};
                                for(int v = 0; v <= (k_mu1 + k_mu2); v++){
                                    std::vector<double> Dv_vec {Dfun(v, p)};
                                    double sum_kpp {0.};
                                    double E3kk_0;
                                    for(int k_pp = v; k_pp >= 0; k_pp -= 2){
                                        if(k_pp <= 4){
                                            E3kk_0 = (k_pp >= k_P)? Efunt0(k_P + k_pp*(k_pp+1)/2, p_3, PAz_3, PBz_3) : Efunt0(k_pp + k_P*(k_P+1)/2, p_3, PBz_3, PAz_3);
                                        } else {
                                            E3kk_0 = (k_pp >= k_P)? EfunTriplet0(5*k_pp + k_P - 25, p_3, PAz_3, PBz_3) : EfunTriplet0(5*k_P + k_pp - 25, p_3, PBz_3, PAz_3);
                                        }
                                        sum_kpp += Dv_vec[k_pp/2]*E3kk_0;
                                    }
                                    sum_v += Ekk_vec(v)*sum_kpp;
                                }

                                overlap3_g_pre1 += g_mu1*g_mu2*g_P*sum_t*sum_u*sum_v;

                            }
                        }
                    }
                    overlap3_g_pre1 *= FAC3_SCF[muind1][gaussC_mu1]*FAC3_SCF[muind2][gaussC_mu2]*FAC3_AUX[Pind][gaussC_P]*std::pow(p_3,-1.5)*std::exp(-(exponent_mu1*exponent_mu2*norm_mu1mu2/p) -(p*exponent_P*norm_PP_3/p_3));
                    overlap3_g_pre0 += overlap3_g_pre1;
                }
            }
        }

        overlap3_g_pre0 *= FAC12_AUX[Pind]*FAC12_SCF[muind1]*FAC12_SCF[muind2]*std::pow(PI,1.5);
        if(std::abs(overlap3_g_pre0) > etol){
            overlap3Matrices_val.push_back(overlap3_g_pre0);
            overlap3Matrices_ind.push_back({Pind,muind1,muind2,Rind1,Rind2});
        }

    }

    auto end = std::chrono::high_resolution_clock::now(); 
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 

uint64_t n_entries = overlap3Matrices_val.size();
std::ofstream output_file(IntFiles_Dir + o3Mat_name + ".o3c");
output_file << "Entry, P, mu, mu', R, R'" << std::endl;
output_file << n_entries << std::endl;
output_file << "Matrix density: " << ((double)n_entries/total_elem)*100 << " %" << std::endl;
for(uint64_t ent = 0; ent < n_entries; ent++){
    output_file << overlap3Matrices_val[ent] << ", " << overlap3Matrices_ind[ent][0] << ", " << overlap3Matrices_ind[ent][1] << ", " <<
        overlap3Matrices_ind[ent][2] << ", " << overlap3Matrices_ind[ent][3] << ", " << overlap3Matrices_ind[ent][4] << std::endl;
}

// overlap3Matrices.save(IntFiles_Dir + o3Mat_name + ".o3c",arma::arma_ascii);  //save matrices to file

std::cout << "Done! Elapsed wall-clock time: " << std::to_string( elapsed.count() * 1e-3 ) << " seconds. Values above 10^-" << std::to_string(tol) <<
    " stored in the file: " << IntFiles_Dir + o3Mat_name << ".o3c." << std::endl;

}

/**
 * Method to compute the E^{i,i'}_{0} coefficients, for (5 <= i <= 8, i' <= 4). Returns only the t=0 component.
 */
double Overlap3Centers::EfunTriplet0(const int index, const double p, const double PA, const double PB){

    switch(index)
    {
    case 0:  { // (i,j) = (5,0)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return (PA*(PAPAp*4*(PAPAp + 5) + 15)*facp*facp);
    } 
    case 1:  { // (i,j) = (5,1)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        return ((PAPAp*(PAPBp*8*(PAPAp + 5) + PAPAp*20 + 60) + 30*PAPBp + 15)*facp*facp*facp);
    }
    case 2:  { // (i,j) = (5,2)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        return ((PAPAp*((8*PAPBp*PAPBp + 4*PAPAp + 60)*PA + 40*PAPBp*(PA + PB) + 120*PB) + 30*PB*(PAPBp + 1) + 75*PA)*facp*facp*facp);
    }
    case 3:  { // (i,j) = (5,3)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPBp*(16*PBPBp*(PAPAp + 5) + 24*PAPAp + 120*PAPBp + 360) + 60*PAPAp + 360*PBPBp + 300) + PAPBp*(60*PBPBp + 450) + 90*PBPBp + 105)*facp2*facp2);
    }
    case 4:  { // (i,j) = (5,4)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(16*PAPBp*(PB*(PAPAp*(PBPBp + 3) + 10*PAPBp + 5*PBPBp + 45) + 15*PA) + 12*PA*(PAPAp + 25) + PB*(480*PBPBp + 1200)) + 60*PBPBp*(PA*(PBPBp + 15) + 2*PB) + 525*PA + 420*PB)*facp2*facp2);
    }
    case 5:  { // (i,j) = (6,0)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PAPAp*(PAPAp*(8*PAPAp + 60) + 90) + 15)*facp*facp*facp);
    } 
    case 6:  { // (i,j) = (6,1)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PAPAp*(PAPAp*(PB*(8*PAPAp + 60) + 24*PA) + 120*PA + 90*PB) + 90*PA + 15*PB)*facp*facp*facp);
    }
    case 7:  { // (i,j) = (6,2)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(8*PAPAp*(2*PBPBp + 1) + 96*PAPBp + 120*PBPBp + 180) + 480*PAPBp + 180*PBPBp + 450) + 30*PBPBp + 360*PAPBp + 105)*facp2*facp2);
    }
    case 8:  { // (i,j) = (6,3)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*PB*(16*PBPBp + 24) + 72*PA*(2*PBPBp + 1) + PB*(120*PBPBp + 540)) + 180*PBPBp*(PB + 4*PA) + 600*PA + 1350*PB) + PBPBp*(540*PA + 30*PB) + 630*PA + 315*PB)*facp2*facp2);
    }
    case 9:  { // (i,j) = (6,4)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*(32*PBPBp*(PBPBp + 3) + 24) + PBPBp*(384*PAPBp + 240*PBPBp + 2160) + 576*PAPBp + 900) + PBPBp*(1920*PAPBp + 360*PBPBp + 5400) + 4800*PAPBp + 3150) + PBPBp*(1440*PAPBp + 60*PBPBp + 1260) + 5040*PAPBp + 945)*facp2*facp2*facp);
    }
    case 10: { // (i,j) = (7,0)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PA*(PAPAp*(8*PAPAp*PAPAp + 84*PAPAp + 210) + 105))*facp*facp*facp);
    } 
    case 11: { // (i,j) = (7,1)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*(16*PAPBp + 56) + 168*PAPBp + 420) + 420*PAPBp + 630) + 210*PAPBp + 105)*facp2*facp2);
    }
    case 12: { // (i,j) = (7,2)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PA*PAPAp*(PAPAp*(PAPAp*(16*PBPBp + 8) + 112*PAPBp + 168*PBPBp + 252) + 840*PAPBp + 420*PBPBp + 1050) + PAPBp*(1260*PA + 210*PB + 735) + 210*PB)*facp2*facp2);
    }
    case 13: { // (i,j) = (7,3)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*(PAPBp*(32*PBPBp + 48) + 336*PBPBp + 168) + PAPBp*(336*PBPBp + 1512) + 2520*PBPBp + 2100) + PAPBp*(840*PBPBp + 6300) + 3780*PBPBp + 4410) + PAPBp*(420*PBPBp + 4410) + 630*PBPBp + 945)*facp2*facp2*facp);
    }
    case 14: { // (i,j) = (7,4)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PA*(PAPAp*(PAPAp*(32*PBPBp*(PBPBp + 3) + 24) + 336*PBPBp*(PBPBp + 9) + 1260) + 840*PBPBp*(PBPBp + 15) + 7350) + PB*(PAPAp*(PAPAp*(448*PBPBp + 672) + 3360*PBPBp + 8400) + 5040*PBPBp + 17640)) + PA*(420*PBPBp*(PBPBp + 21) + 6615) + PB*(840*PBPBp + 3780))*facp2*facp2*facp);
    }
    case 15: { // (i,j) = (8,0)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        return ((PAPAp*(PAPAp*(16*PAPAp*(PAPAp + 14) + 840) + 840) + 105)*facp2*facp2);
    }
    case 16: { // (i,j) = (8,1)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        return ((PAPAp*(PA*(PAPAp*(64*PAPAp + 672) + 1680) + PB*(PAPAp*(16*PAPAp*(PAPAp + 14) + 840) + 840)) + 840*PA + 105*PB)*facp2*facp2);
    }
    case 17: { // (i,j) = (8,2)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*(PAPAp*(32*PBPBp + 16) + 256*PAPBp + 448*PBPBp + 672) + 2688*PAPBp + 1680*PBPBp + 4200) + 6720*PAPBp + 1680*PBPBp + 5880) + 3360*PAPBp + 210*PBPBp + 945)*facp2*facp2*facp);
    }
    case 18: { // (i,j) = (8,3)
        double facp = 0.5/p;
        double facp2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PB*(PAPAp*(PAPAp*(32*PBPBp*(PAPAp + 14) + 48*PAPAp + 2016) + 1680*PBPBp + 12600) + 1680*PBPBp + 17640) + PA*(PAPAp*(PAPAp*(384*PBPBp + 192) + 4032*PBPBp + 3360) + 10080*PBPBp + 11760)) + PB*(210*PBPBp + 2835) + PA*(5040*PBPBp + 7560))*facp2*facp2*facp);
    }
    case 19: { // (i,j) = (8,4)
        double facp = 0.5/p;
        double facp3 = facp*facp*facp;
        double PAPAp = PA*PA*p;
        double PAPBp = PA*PB*p;
        double PBPBp = PB*PB*p;
        return ((PAPAp*(PAPAp*(PAPAp*(PAPAp*(64*PBPBp*(PBPBp + 3) + 48) + PBPBp*(1024*PAPBp + 896*PBPBp + 8064) + 1536*PAPBp + 3360) + PBPBp*(10752*PAPBp + 3360*PBPBp + 50400) + 26880*PAPBp + 29400)+ 3360*PBPBp*(8*PAPBp + PBPBp + 21) + 94080*PAPBp + 52920) + PAPBp*(13440*PBPBp + 60480) + 420*PBPBp*(PBPBp + 27) + 10395)*facp3*facp3);
    }
    default: {
        throw std::invalid_argument("Overlap3Centers_Mat::EfunTriplet0 error: the E^{i,i'}_{0} coefficients are being evaluated for i >= 9 and/or i' >= 5");
    }
    }

}

// Coefficients D^{t}_{l}(p) in the expansion \Lambda_{t}(x,p,Px) = D^{t}_{l}(p) *(x-Px)^l *exp(-p(x-Px)^2).
// The argument t spans 0 <= t <= 8, and the entries of the returned vector are the nonzero D^{t}_{l} for the given t.  
// There are ceil((t+1)/2) of such entries, each corresponding to l = t, t-2, t-4, ... 1 (or 0)
std::vector<double> Overlap3Centers::Dfun(const int t, const double p){

    switch(t)
    {
    case 0:  {
        return std::vector<double> {1.0};
    } 
    case 1:  { 
        return std::vector<double> {2*p};
    }
    case 2:  { 
        double fac = -2*p;
        return std::vector<double> {fac, fac*fac};
    }
    case 3:  { 
        double fac = p*p;
        return std::vector<double> {-12*fac, 8*p*fac};
    }
    case 4:  { 
        double fac = 4*p*p;
        return std::vector<double> {3*fac, -12*p*fac, fac*fac};
    }
    case 5:  { 
        double fac = p*p*p;
        return std::vector<double> {120*fac, -160*p*fac, 32*p*p*fac};
    }
    case 6:  { 
        double fac = 8*p*p*p;
        return std::vector<double> {-15*fac, 90*p*fac, -60*p*p*fac, fac*fac};
    }
    case 7:  { 
        double fac = p*p*p*p;
        return std::vector<double> {-1680*fac, 3360*p*fac, -1344*p*p*fac, 128*p*p*p*fac};
    }
    case 8:  { 
        double fac = 16*p*p*p*p;
        return std::vector<double> {105*fac, -840*p*fac, 840*p*p*fac, -224*p*p*p*fac, fac*fac};
    }
    default: {
        throw std::invalid_argument("Overlap3Centers::Dfun error: the D^{t}_{l} coefficients are being evaluated for t not in [0,8]");
    }
    }

}

}