#include "xatu/Overlap2Centers.hpp"
#include <chrono>

namespace xatu {

/**
 * Constructor that also initializes IntegralsBase from the DFT files. 
 * Should only be used to test the 2-center overlap matrices in isolation from the 3-center ones.
 */
Overlap2Centers::Overlap2Centers(const GTFConfiguration& GTFConfig, const std::string& overlap2Matrices_filename, const bool comp) : IntegralsBase{GTFConfig} {

    if(comp){
        overlap2Cfun(ncells, overlap2Matrices_filename);
    }

}

/**
 * Constructor that copies a pre-initialized IntegralsBase.
 */
Overlap2Centers::Overlap2Centers(const IntegralsBase& IntBase, const std::string& overlap2Matrices_filename) : IntegralsBase{IntBase} {

    overlap2Cfun(ncells, overlap2Matrices_filename);

}

/**
 * Method to compute the overlap matrices in the auxiliary basis <P,0|P',R> for the first nR Bravais vectors R, 
 * where nR <= ncells (attribute of IntegralsBase and third argument of GTFConfiguration 's constructor).
 * The resulting cube (third dimension spans the Bravais vectors) is saved in the overlap2Matrices_filename.o2c file.
 */
void Overlap2Centers::overlap2Cfun(const int nR, const std::string& overlap2Matrices_filename){

std::cout << "Computing " << nR << " " << dimMat_AUX << "x" << dimMat_AUX << " 2-center overlap matrices..." << std::endl;
auto begin = std::chrono::high_resolution_clock::now();  

    int dimMat = dimMat_AUX;
    int nelem_triang = 0.5*dimMat*(dimMat + 1);
    arma::cube overlap2Matrices {arma::zeros<arma::cube>(dimMat,dimMat,nR)};

    #pragma omp parallel for
    for(long int s = 0; s < (nelem_triang*nR); s++){ //Spans the lower triangle of all the nR matrices <P,0|P',R>
        int sind {s % nelem_triang};   //Index for the corresponding entry in the overlap matrix, irrespective of the specific R
        int Rind {s / nelem_triang};   //Position in RlistAU (i.e. 0 for R=0) of the corresponding Bravais vector 
        std::array<int,2> orb_braket {triangInd_to_rowcol.at(sind)}; 
        arma::colvec R {RlistAU.col(Rind)};  //Bravais vector (a.u.) corresponding to the "s" matrix element
        int RindOpp    {RlistOpposites.at(Rind)};  //Position in RlistAU (i.e. 0 for R=0) of the opposite of the corresponding Bravais vector 

        int orb_bra {orb_braket[0]};   //Orbital number (<dimMat_AUX) of the bra corresponding to the index s 
        std::vector<int> orb_info_int_bra {orbitals_info_int_AUX[orb_bra]};
        int L_bra  {orb_info_int_bra[2]};
        int m_bra  {orb_info_int_bra[3]};
        int nG_bra {orb_info_int_bra[4]};
        std::vector<double> orb_info_real_bra {orbitals_info_real_AUX[orb_bra]};
        arma::colvec coords_bra {orb_info_real_bra[0], orb_info_real_bra[1], orb_info_real_bra[2]};  //Position (a.u.) of bra atom
        std::vector<double> FAC3_bra   {FAC3_AUX[orb_bra]};
        std::vector<int> g_coefs_bra   {g_coefs.at( L_bra*(L_bra + 1) + m_bra )};

        int orb_ket {orb_braket[1]};   //Orbital number (<dimMat_AUX) of the ket corresponding to the index s. orb_ket <= orb_bra (lower triangle)
        std::vector<int> orb_info_int_ket {orbitals_info_int_AUX[orb_ket]};
        int L_ket  {orb_info_int_ket[2]};
        int m_ket  {orb_info_int_ket[3]};
        int nG_ket {orb_info_int_ket[4]};
        std::vector<double> orb_info_real_ket {orbitals_info_real_AUX[orb_ket]};
        arma::colvec coords_ket {R + arma::colvec{orb_info_real_ket[0], orb_info_real_ket[1], orb_info_real_ket[2]} };  //Position (a.u.) of ket atom
        std::vector<double> FAC3_ket   {FAC3_AUX[orb_ket]};
        std::vector<int> g_coefs_ket   {g_coefs.at( L_ket*(L_ket + 1) + m_ket )};

        arma::colvec coords_braket {coords_bra - coords_ket};
        double norm_braket {arma::dot(coords_braket,coords_braket)};
        double FAC12_braket = FAC12_AUX[orb_bra]*FAC12_AUX[orb_ket];

        for(int gaussC_bra = 0; gaussC_bra < nG_bra; gaussC_bra++){
            double exponent_bra {orb_info_real_bra[2*gaussC_bra + 3]};
            //double d_bra {orb_info_real_bra[2*gaussC_bra + 4]};
            double FAC3_gaussC_bra {FAC3_bra[gaussC_bra]};

            for(int gaussC_ket = 0; gaussC_ket < nG_ket; gaussC_ket++){
                double exponent_ket {orb_info_real_ket[2*gaussC_ket + 3]};
                //double d_ket {orb_info_real_ket[2*gaussC_ket + 4]};
                double FAC3_gaussC_ket {FAC3_ket[gaussC_ket]};

                double p   {exponent_bra + exponent_ket};  //exponent coefficient of the Hermite Gaussian
                arma::colvec P {(exponent_bra*coords_bra + exponent_ket*coords_ket)/p};  //center of the Hermite Gaussian
                double PAx {P(0) - coords_bra(0)}; 
                double PAy {P(1) - coords_bra(1)}; 
                double PAz {P(2) - coords_bra(2)}; 
                double PBx {P(0) - coords_ket(0)}; 
                double PBy {P(1) - coords_ket(1)}; 
                double PBz {P(2) - coords_ket(2)}; 

                double overlap_g_pre {0};
                std::vector<int>::iterator g_itr_bra {g_coefs_bra.begin()};
                for(int numg_bra = 0; numg_bra < g_coefs_bra[0]; numg_bra++){
                    int i_bra {*(++g_itr_bra)};
                    int j_bra {*(++g_itr_bra)};
                    int k_bra {*(++g_itr_bra)};
                    int g_bra {*(++g_itr_bra)};
                    int Ei_bra {i_bra*(i_bra + 1)/2};
                    int Ej_bra {j_bra*(j_bra + 1)/2};
                    int Ek_bra {k_bra*(k_bra + 1)/2};
                    
                    std::vector<int>::iterator g_itr_ket {g_coefs_ket.begin()};
                    for(int numg_ket = 0; numg_ket < g_coefs_ket[0]; numg_ket++){
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
                overlap_g_pre *= FAC3_gaussC_bra*FAC3_gaussC_ket*std::pow(p,-1.5)*std::exp(-exponent_bra*exponent_ket*norm_braket/p);
                overlap2Matrices(orb_bra,orb_ket,Rind) += overlap_g_pre;
            }
        }
        overlap2Matrices(orb_bra,orb_ket,Rind) *= FAC12_braket*std::pow(PI,1.5);
        if((RindOpp >= 0) && (orb_bra > orb_ket)){
            overlap2Matrices(orb_ket,orb_bra,RindOpp) = overlap2Matrices(orb_bra,orb_ket,Rind);
        }

    }

    for(int Rind : R_unpaired){ //compute the upper triangles of the R-matrices corresponding to unpaired Bravais vectors 
        int nelem_triang_offd {nelem_triang - dimMat};

        #pragma omp parallel for
        for(int s = 0; s < nelem_triang_offd; s++){
            std::array<int,2> orb_braket {triangInd_to_rowcol.at(s)};
            arma::colvec R {RlistAU.col(Rind)};

            int orb_bra {orb_braket[1]};   
            std::vector<int> orb_info_int_bra {orbitals_info_int_AUX[orb_bra]};
            int L_bra  {orb_info_int_bra[2]};
            int m_bra  {orb_info_int_bra[3]};
            int nG_bra {orb_info_int_bra[4]};
            std::vector<double> orb_info_real_bra {orbitals_info_real_AUX[orb_bra]};
            arma::colvec coords_bra {orb_info_real_bra[0], orb_info_real_bra[1], orb_info_real_bra[2]};  
            std::vector<double> FAC3_bra   {FAC3_AUX[orb_bra]};
            std::vector<int> g_coefs_bra   {g_coefs.at( L_bra*(L_bra + 1) + m_bra )};

            int orb_ket {orb_braket[0] + 1};   
            std::vector<int> orb_info_int_ket {orbitals_info_int_AUX[orb_ket]};
            int L_ket  {orb_info_int_ket[2]};
            int m_ket  {orb_info_int_ket[3]};
            int nG_ket {orb_info_int_ket[4]};
            std::vector<double> orb_info_real_ket {orbitals_info_real_AUX[orb_ket]};
            arma::colvec coords_ket {R + arma::colvec{orb_info_real_ket[0], orb_info_real_ket[1], orb_info_real_ket[2]} }; 
            std::vector<double> FAC3_ket   {FAC3_AUX[orb_ket]};
            std::vector<int> g_coefs_ket   {g_coefs.at( L_ket*(L_ket + 1) + m_ket )};

            arma::colvec coords_braket {coords_bra - coords_ket};
            double norm_braket {arma::dot(coords_braket,coords_braket)};
            double FAC12_braket = FAC12_AUX[orb_bra]*FAC12_AUX[orb_ket];

            for(int gaussC_bra = 0; gaussC_bra < nG_bra; gaussC_bra++){
            double exponent_bra {orb_info_real_bra[2*gaussC_bra + 3]};
            //double d_bra {orb_info_real_bra[2*gaussC_bra + 4]};
            double FAC3_gaussC_bra {FAC3_bra[gaussC_bra]};

            for(int gaussC_ket = 0; gaussC_ket < nG_ket; gaussC_ket++){
                double exponent_ket {orb_info_real_ket[2*gaussC_ket + 3]};
                //double d_ket {orb_info_real_ket[2*gaussC_ket + 4]};
                double FAC3_gaussC_ket {FAC3_ket[gaussC_ket]};

                double p   {exponent_bra + exponent_ket};  //exponent coefficient of the Hermite Gaussian
                arma::colvec P {(exponent_bra*coords_bra + exponent_ket*coords_ket)/p};  //center of the Hermite Gaussian
                double PAx {P(0) - coords_bra(0)}; 
                double PAy {P(1) - coords_bra(1)}; 
                double PAz {P(2) - coords_bra(2)}; 
                double PBx {P(0) - coords_ket(0)}; 
                double PBy {P(1) - coords_ket(1)}; 
                double PBz {P(2) - coords_ket(2)}; 

                double overlap_g_pre {0};
                std::vector<int>::iterator g_itr_bra {g_coefs_bra.begin()};
                for(int numg_bra = 0; numg_bra < g_coefs_bra[0]; numg_bra++){
                    int i_bra {*(++g_itr_bra)};
                    int j_bra {*(++g_itr_bra)};
                    int k_bra {*(++g_itr_bra)};
                    int g_bra {*(++g_itr_bra)};
                    int Ei_bra {i_bra*(i_bra + 1)/2};
                    int Ej_bra {j_bra*(j_bra + 1)/2};
                    int Ek_bra {k_bra*(k_bra + 1)/2};
                    
                    std::vector<int>::iterator g_itr_ket {g_coefs_ket.begin()};
                    for(int numg_ket = 0; numg_ket < g_coefs_ket[0]; numg_ket++){
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
                overlap_g_pre *= FAC3_gaussC_bra*FAC3_gaussC_ket*std::pow(p,-1.5)*std::exp(-exponent_bra*exponent_ket*norm_braket/p);
                overlap2Matrices(orb_bra,orb_ket,Rind) += overlap_g_pre;
            }
        }
        overlap2Matrices(orb_bra,orb_ket,Rind) *= FAC12_braket*std::pow(PI,1.5);
        }
    }

auto end = std::chrono::high_resolution_clock::now(); 
auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 

overlap2Matrices.save(overlap2Matrices_filename + ".o2c",arma::arma_ascii);  //save matrices to file

std::cout << "Done! Elapsed wall-clock time: " << std::to_string( elapsed.count() * 1e-3 ) << " seconds. Matrices stored in the file: " << 
    overlap2Matrices_filename << ".o2c." << std::endl;

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

    if(index == 0){ // (i,j) = (0,0)
        return 1.0;
    } 
    else if(index == 1) { // (i,j) = (1,0)
        return PA;
    }
    else if(index == 2) { // (i,j) = (1,1)
        return (PA*PB + 0.5/p);
    }
    else if(index == 3) { // (i,j) = (2,0)
        return (PA*PA + 0.5/p);
    }
    else if(index == 4) { // (i,j) = (2,1)
        double facp = 0.5/p;
        return (PA*2*facp + PB*(PA*PA+facp)); 
    }
    else if(index == 5) { // (i,j) = (2,2)
        double facp = 0.5/p;
        return (facp*(3*facp + PB*(4*PA+PB)) + PA*PA*(facp+PB*PB));  
    }
    else if(index == 6) { // (i,j) = (3,0)
        return (0.5*PA*(PA*PA*p*2 + 3)/p);
    }
    else if(index == 7) { // (i,j) = (3,1)
        double facp = 0.5/p;
        return ((PA*2*p*(3*(PA+PB) + PB*PA*PA*p*2) + 3)*facp*facp);
    }
    else if(index == 8) { // (i,j) = (3,2)
        double facp = 0.5/p;
        return ((PA*2*p*(PA*PA*(PB*PB*p*2 + 1) + 3*PB*(2*PA+PB)) + 9*PA + 6*PB)*facp*facp);
    }
    else if(index == 9) { // (i,j) = (3,3)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return ((PAPBp*2*(PAPAp*2*(PBPBp*2 + 3) + PAPBp*18 + PBPBp*6 + 27) + 18*(PAPAp + PBPBp) + 15)*facp*facp*facp);
    }
    else if(index == 10) { // (i,j) = (4,0)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PAPAp*4*(PAPAp + 3) + 3)*facp*facp);
    }
    else if(index == 11) { // (i,j) = (4,1)
        double facp = 0.5/p;
        return ((PA*PA*p*4*((PA*PB*p + 2)*PA + 3*PB) + 12*PA + 3*PB)*facp*facp); 
    }
    else if(index == 12) { // (i,j) = (4,2)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        return ((PAPAp*4*(PAPAp*(1 + PB*PB*p*2) + PB*2*p*(4*PA + 3*PB) + 9) + PB*6*p*(8*PA + PB) + 15)*facp*facp*facp);
    }
    else if(index == 13) { // (i,j) = (4,3)
        double facp = 0.5/p;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return ((PAPBp*PAPBp*8*(3*PB + 6*PA + PAPBp*PA) + PA*12*(5 + 2*PAPAp + 6*PBPBp + PAPBp*(PAPAp + 9)) + 3*PB*(15 + PBPBp*2))*facp*facp*facp);
    }
    else if(index == 14) { // (i,j) = (4,4)
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double PAPAp = PA*PA*p;
        double PBPBp = PB*PB*p;
        double PAPBp = PA*PB*p;
        return ((PAPAp*PAPAp*4*(PBPBp*4*(PBPBp + 3) + 3) + PAPBp*PAPAp*64*(PBPBp*2 + 3) + PAPAp*12*(PBPBp*4*(PBPBp + 9) + 15)
        + PAPBp*96*(PBPBp*2 + 5) + PBPBp*12*(PBPBp + 15) + 105)*facp_to2*facp_to2);
    } 
    else {
        throw std::invalid_argument("Overlap2Centers::Efunt0 error: the E^{i,i'}_{0} coefficients are being evaluated for i and/or i' >= 5");
    }
        
}


}