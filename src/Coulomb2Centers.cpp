#include "xatu/Coulomb2Centers.hpp"

namespace xatu {

/**
 * Constructor that also initializes IntegralsBase from the DFT files. 
 * Should only be used to test the 3-center overlap matrices in isolation.
 * @param GTFConfig GTFConfiguration object.
 * @param tol Tolerance for retaining the entries of the 3-center overlap matrices. These must be > 10^-tol, in absolute value.
 * @param o3Mat_name Name of the file where the 3-center overlap matrices will be stored as a cube (o3Mat_name.o3c).
 * @param o2Mat_name Name of the file where the 2-center overlap matrices in the AUX basis will be stored as a cube (o2Mat_name.o2c).
 */
Coulomb2Centers::Coulomb2Centers(const GTFConfiguration& GTFConfig, const uint32_t nR, const std::string& C2Mat_name, const bool comp) : IntegralsBase{GTFConfig} {

    if(comp){
        Coulomb2Cfun(nR, C2Mat_name);
    }

}

/**
 * Constructor that copies a pre-initialized IntegralsBase and also computes the 2-center overlap in the AUX basis.
 * @param IntBase IntegralsBase object.
 * @param tol Tolerance for retaining the entries of the 3-center overlap matrices. These must be > 10^-tol, in absolute value.
 * @param nR2 Number of R and R' that will be considered for the integrals. By default it's IntegralsBase.ncells, and cannot be larger than it.
 * @param o3Mat_name Name of the file where the 3-center overlap matrices will be stored as a cube (o3Mat_name.o3c).
 * @param o2Mat_name Name of the file where the 2-center overlap matrices in the AUX basis will be stored as a cube (o2Mat_name.o2c).
 */
Coulomb2Centers::Coulomb2Centers(const IntegralsBase& IntBase, const uint32_t nR, const std::string& C2Mat_name) : IntegralsBase{IntBase} {

    Coulomb2Cfun(nR, C2Mat_name);

}

/**
 * Method to compute the Coulomb matrices in the auxiliary basis (<P,0|V_c|P',R>) for the first nR Bravais vectors R.
 * These first nR (at least, until the star of vectors is completed) are generated with IntegralsBase::generateRlist.
 * The resulting cube (third dimension spans the Bravais vectors) is saved in the C2Mat_name.C2c file.
 */
void Coulomb2Centers::Coulomb2Cfun(const uint32_t nR, const std::string& C2Mat_name){

arma::mat RlistAU = ANG2AU*generateRlist(bravaisLattice, nR);
uint32_t nR_star = RlistAU.n_cols;
std::map<uint32_t,uint32_t> RlistOpposites = generateRlistOpposite(RlistAU);

std::cout << "Computing " << nR_star << " " << dimMat_AUX << "x" << dimMat_AUX << " 2-center Coulomb matrices in the AUX basis..." << std::endl;
auto begin = std::chrono::high_resolution_clock::now();  

    uint64_t nelem_triang = 0.5*dimMat_AUX*(dimMat_AUX + 1);
    uint64_t total_elem = nelem_triang*nR_star;
    arma::cube Coulomb2Matrices {arma::zeros<arma::cube>(dimMat_AUX,dimMat_AUX,nR_star)};

    #pragma omp parallel for 
    for(uint64_t s = 0; s < total_elem; s++){ //Spans the lower triangle of all the nR_star matrices <P,0|P',R>
        uint64_t sind {s % nelem_triang};    //Index for the corresponding entry in the overlap matrix, irrespective of the specific R
        uint32_t Rind {s / nelem_triang};         //Position in RlistAU (i.e. 0 for R=0) of the corresponding Bravais vector 
        std::array<uint32_t,2> orb_braket {triangInd_to_rowcol.at(sind)}; 
        // arma::colvec R {RlistAU.col(Rind)};  //Bravais vector (a.u.) corresponding to the "s" matrix element
        uint32_t RindOpp    {RlistOpposites.at(Rind)};  //Position in RlistAU (i.e. 0 for R=0) of the opposite of the corresponding Bravais vector

        uint32_t orb_bra {orb_braket[0]};   //Orbital number (<dimMat) of the bra corresponding to the index s 
        int L_bra  {orbitals_info_int_AUX[orb_bra][2]};
        int m_bra  {orbitals_info_int_AUX[orb_bra][3]};
        int nG_bra {orbitals_info_int_AUX[orb_bra][4]};
        arma::colvec coords_bra {orbitals_info_real_AUX[orb_bra][0], orbitals_info_real_AUX[orb_bra][1], orbitals_info_real_AUX[orb_bra][2]};  //Position (a.u.) of bra atom
        std::vector<int> g_coefs_bra   {g_coefs.at( L_bra*(L_bra + 1) + m_bra )};

        uint32_t orb_ket {orb_braket[1]};   //Orbital number (<dimMat) of the ket corresponding to the index s. orb_ket <= orb_bra (lower triangle)
        int L_ket  {orbitals_info_int_AUX[orb_ket][2]};
        int m_ket  {orbitals_info_int_AUX[orb_ket][3]};
        int nG_ket {orbitals_info_int_AUX[orb_ket][4]};
        arma::colvec coords_ket {RlistAU.col(Rind) + arma::colvec{orbitals_info_real_AUX[orb_ket][0], orbitals_info_real_AUX[orb_ket][1], orbitals_info_real_AUX[orb_ket][2]} };  //Position (a.u.) of ket atom
        std::vector<int> g_coefs_ket   {g_coefs.at( L_ket*(L_ket + 1) + m_ket )};

        arma::colvec coords_braket {coords_bra - coords_ket};
        double FAC12_braket = FAC12_AUX[orb_bra]*FAC12_AUX[orb_ket]; 

        for(int gaussC_bra = 0; gaussC_bra < nG_bra; gaussC_bra++){ //Iterate over the contracted Gaussians in the bra orbital
            double exponent_bra {orbitals_info_real_AUX[orb_bra][2*gaussC_bra + 3]};
            //double d_bra {orbitals_info_real_AUX[orb_bra][2*gaussC_bra + 4]};

            for(int gaussC_ket = 0; gaussC_ket < nG_ket; gaussC_ket++){ //Iterate over the contracted Gaussians in the ket orbital
                double exponent_ket {orbitals_info_real_AUX[orb_ket][2*gaussC_ket + 3]};
                //double d_ket {orbitals_info_real_AUX[orb_ket][2*gaussC_ket + 4]};

                double p {exponent_bra + exponent_ket};  //Exponent coefficient of the Hermite Gaussian
                double mu {exponent_bra*exponent_ket/p};

                double Coulomb_g_pre {0.};
                std::vector<int>::iterator g_itr_bra {g_coefs_bra.begin()};
                for(int numg_bra = 0; numg_bra < g_coefs_bra[0]; numg_bra++){ //Iterate over the summands of the corresponding spherical harmonic in the bra orbital
                    int i_bra {*(++g_itr_bra)};
                    int j_bra {*(++g_itr_bra)};
                    int k_bra {*(++g_itr_bra)};
                    int g_bra {*(++g_itr_bra)};

                    arma::colvec Ei0vec_bra {Efun_single(i_bra, exponent_bra)};
                    arma::colvec Ej0vec_bra {Efun_single(j_bra, exponent_bra)};
                    arma::colvec Ek0vec_bra {Efun_single(k_bra, exponent_bra)};
                    
                    std::vector<int>::iterator g_itr_ket {g_coefs_ket.begin()};
                    for(int numg_ket = 0; numg_ket < g_coefs_ket[0]; numg_ket++){ //Iterate over the summands of the corresponding spherical harmonic in the ket orbital
                        int i_ket {*(++g_itr_ket)};
                        int j_ket {*(++g_itr_ket)};
                        int k_ket {*(++g_itr_ket)};
                        int g_ket {*(++g_itr_ket)};

                        arma::colvec Ei0vec_ket {Efun_single(i_ket, exponent_ket)};
                        arma::colvec Ej0vec_ket {Efun_single(j_ket, exponent_ket)};
                        arma::colvec Ek0vec_ket {Efun_single(k_ket, exponent_ket)};

                        for(int t_bra = 0; t_bra <= i_bra; t_bra++){
                            for(int u_bra = 0; u_bra <= j_bra; u_bra++){
                                for(int v_bra = 0; v_bra <= k_bra; v_bra++){
                                    double Eijk_bra = Ei0vec_bra(t_bra)*Ej0vec_bra(u_bra)*Ek0vec_bra(v_bra);

                                    for(int t_ket = 0; t_ket <= i_ket; t_ket++){
                                        for(int u_ket = 0; u_ket <= j_ket; u_ket++){
                                            for(int v_ket = 0; v_ket <= k_ket; v_ket++){
                                                int sign_ket {std::pow(-1, t_ket + u_ket + v_ket)};
                                                int t_tot = t_bra + t_ket;
                                                int u_tot = u_bra + u_ket;
                                                int v_tot = v_bra + v_ket;
                                                double Hermit;
                                                if(t_tot >= u_tot){
                                                    if(u_tot >= v_tot){ // (t,u,v)
                                                        Hermit = HermiteCoulomb(t_tot, u_tot, v_tot, mu, coords_braket(0), coords_braket(1), coords_braket(2));
                                                    }
                                                    else if(t_tot >= v_tot){ // (t,v,u)
                                                        Hermit = HermiteCoulomb(t_tot, v_tot, u_tot, mu, coords_braket(0), coords_braket(2), coords_braket(1));
                                                    }
                                                    else{ // (v,t,u)
                                                        Hermit = HermiteCoulomb(v_tot, t_tot, u_tot, mu, coords_braket(2), coords_braket(0), coords_braket(1));
                                                    }
                                                } 
                                                else if(u_tot >= v_tot){ 
                                                    if(t_tot >= v_tot){ // (u,t,v)
                                                        Hermit = HermiteCoulomb(u_tot, t_tot, v_tot, mu, coords_braket(1), coords_braket(0), coords_braket(2));
                                                    }
                                                    else{ // (u,v,t)
                                                        Hermit = HermiteCoulomb(u_tot, v_tot, t_tot, mu, coords_braket(1), coords_braket(2), coords_braket(0));
                                                    }
                                                }
                                                else{ // (v,u,t)
                                                    Hermit = HermiteCoulomb(v_tot, u_tot, t_tot, mu, coords_braket(2), coords_braket(1), coords_braket(0));
                                                }
                                                
                                                Coulomb_g_pre += g_bra*g_ket*sign_ket*Eijk_bra*Ei0vec_ket(t_ket)*Ej0vec_ket(u_ket)*Ek0vec_ket(v_ket)*Hermit;

                                            }
                                        }
                                    }

                                }
                            }
                        }


                    }
                }
                Coulomb_g_pre *= FAC3_AUX[orb_bra][gaussC_bra]*FAC3_AUX[orb_ket][gaussC_ket]*std::pow(p,-1.5)/mu;
                Coulomb2Matrices(orb_bra,orb_ket,Rind) += Coulomb_g_pre;

            }
        }
        Coulomb2Matrices(orb_bra,orb_ket,Rind) *= FAC12_braket*2*std::pow(PI,2.5);
        if((RindOpp >= 0) && (orb_bra > orb_ket)){
            Coulomb2Matrices(orb_ket,orb_bra,RindOpp) = Coulomb2Matrices(orb_bra,orb_ket,Rind);
        }

    }

    auto end = std::chrono::high_resolution_clock::now(); 
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin); 

    Coulomb2Matrices.save(IntFiles_Dir + C2Mat_name + ".C2c",arma::arma_ascii);  //save matrices to file

    std::cout << "Done! Elapsed wall-clock time: " << std::to_string( elapsed.count() * 1e-3 ) << " seconds. Matrices stored in the file: " << 
        IntFiles_Dir + C2Mat_name << ".C2c." << std::endl;

}

/**
 * Analogous to Efun in the parent class IntegralsBase, but restricted to i'=0 and setting PA=0. These are the 
 * expansion coefficients of a single GTF in Hermite Gaussians with the same exponent and center.
 */
arma::colvec Coulomb2Centers::Efun_single(const int i, const double p){

    switch(i)
    {
    case 0:  { 
        return arma::colvec {1.0};
    } 
    case 1:  { 
        return arma::colvec {0, 0.5/p};
    }
    case 2:  { 
        double facp = 0.5/p;
        return arma::colvec {facp,  0,  facp*facp};
    }
    case 3:  { 
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        return arma::colvec {0,  3*facp_to2,  0,  facp_to2*facp};
    }
    case 4:  { 
        double facp = 0.5/p;
        double facp_to2 = facp*facp;
        double facp_to3 = facp_to2*facp;
        return arma::colvec {3*facp_to2,  0,  6*facp_to3,  0,  facp_to3*facp};
    }
    default: {
        throw std::invalid_argument("Coulomb2Centers::Efun_single error: the E^{i,0}_{t} coefficients are being evaluated for i >= 5");
    }
    }

}

/**
 * Method to compute and return the Boys function F_{n}(arg) = \int_{0}^{1}t^{2n}exp(-arg*t^2)dt. It is computed 
 * with the lower incomplete Gamma function as: F_{n}(arg) = Gamma(n+0.5)*IncGamma(n+0.5,arg)/(2*arg^(n+0.5)), see (9.8.20)-Helgaker.
 */
double Coulomb2Centers::Boysfun(const int n, const double arg){

    if(arg < 3.1e-10){
        return (double)1/(2*n + 1);
    }
    else if(arg > 19.785){
        return doubleFactorial(2*n-1)*std::pow(2,-n-1)*std::sqrt(PI*std::pow(arg,-2*n-1));
    }
    else{
        double nfac = n + 0.5;
        return asa239::gammad(arg, nfac)*std::tgamma(nfac)*0.5*std::pow(arg,-nfac);
    }

}

/**
 * Method to compute and return the auxiliary Hermite Coulomb integral R^{n}_{0,0,0}(p,arg)=(-2p)^n *F_{n}(arg), see (9.9.14)-Helgaker.
 */
double Coulomb2Centers::Rn000(const int n, const double p, const double arg){

    return std::pow(-2*p,n)*Boysfun( n, arg );

}

/**
 * Method to compute and return the Hermite Coulomb integral R^{0}_{t,u,v}(r,(X,Y,Z)), see (9.9.9)-Helgaker.
 */
double Coulomb2Centers::HermiteCoulomb(const int t, const int u, const int v, const double p, const double X, const double Y, const double Z){

    double arg = p*(X*X + Y*Y + Z*Z);
    switch(v)
    {
    case 0:  {
        switch(u)
        {
        case 0:  {
            switch(t)
            {
            case 0:  {
                return Boysfun(0, arg);
            } 
            case 1:  {
                return X*Rn000(1,p,arg);
            }
            case 2:  {
                return Rn000(1,p,arg) + X*X*Rn000(2,p,arg);
            }
            case 3:  {
                return X*(3*Rn000(2,p,arg) + X*X*Rn000(3,p,arg));
            }
            case 4:  {
                return 3*Rn000(2,p,arg) + X*X*(6*Rn000(3,p,arg) + X*X*Rn000(4,p,arg));
            }
            case 5:  {
                double X2 = X*X;
                return X*(15*Rn000(3,p,arg) + X2*(10*Rn000(4,p,arg) + X2*Rn000(5,p,arg)));
            }
            case 6:  {
                double X2 = X*X;
                return 15*Rn000(3,p,arg) + X2*(45*Rn000(4,p,arg) + X2*(15*Rn000(5,p,arg) + X2*Rn000(6,p,arg)));
            } 
            case 7:  {
                double X2 = X*X;
                return X*(105*Rn000(4,p,arg) + X2*(105*Rn000(5,p,arg) + X2*(21*Rn000(6,p,arg) + X2*Rn000(7,p,arg))));
            }
            case 8:  {
                double X2 = X*X;
                return 105*Rn000(4,p,arg) + X2*(420*Rn000(5,p,arg) + X2*(210*Rn000(6,p,arg) + X2*(28*Rn000(7,p,arg) + X2*Rn000(8,p,arg))));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=0,v=0)");
            }
            }
        }
        case 1:  {
            switch(t)
            {
            case 1:  {
                return X*Y*Rn000(2,p,arg);
            }
            case 2:  {
                return Y*(Rn000(2,p,arg) + X*X*Rn000(3,p,arg));
            }
            case 3:  {
                return X*Y*(3*Rn000(3,p,arg) + X*X*Rn000(4,p,arg));
            }
            case 4:  {
                return Y*(3*Rn000(3,p,arg) + X*X*(6*Rn000(4,p,arg) + X*X*Rn000(5,p,arg)));
            }
            case 5:  {
                return X*Y*(15*Rn000(4,p,arg) + X*X*(10*Rn000(5,p,arg) + X*X*Rn000(6,p,arg)));
            }
            case 6:  {
                double X2 = X*X;
                return Y*(15*Rn000(4,p,arg) + X2*(45*Rn000(5,p,arg) + X2*(15*Rn000(6,p,arg) + X2*Rn000(7,p,arg))));
            }
            case 7:  {
                double X2 = X*X;
                return X*Y*(105*Rn000(5,p,arg) + X2*(105*Rn000(6,p,arg) + X2*(21*Rn000(7,p,arg) + X2*Rn000(8,p,arg))));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=1,v=0)");
            }
            }
        }
        case 2:  {
            switch(t)
            {
            case 2:  {
                return Rn000(2,p,arg) + (X*X + Y*Y)*Rn000(3,p,arg) + X*X*Y*Y*Rn000(4,p,arg);
            }
            case 3:  {
                return X*(3*Rn000(3,p,arg) + (X*X + 3*Y*Y)*Rn000(4,p,arg) + X*X*Y*Y*Rn000(5,p,arg));
            }
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return 3*Rn000(3,p,arg) + 3*(2*X2 + Y2)*Rn000(4,p,arg) + X2*((X2 + 6*Y2)*Rn000(5,p,arg) + X2*Y2*Rn000(6,p,arg));
            }
            case 5:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return X*(15*Rn000(4,p,arg) + (10*X2 + 15*Y2)*Rn000(5,p,arg) + X2*((X2 + 10*Y2)*Rn000(6,p,arg) + X2*Y2*Rn000(7,p,arg)));
            }
            case 6:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return 15*Rn000(4,p,arg) + 15*(3*X2 + Y2)*Rn000(5,p,arg) + 15*X2*((X2 + 3*Y2)*Rn000(6,p,arg) + X2*((X2 + 15*Y2)*Rn000(7,p,arg) + X2*Y2*Rn000(8,p,arg)));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=2,v=0)");
            }
            }
        }
        case 3:  {
            switch(t)
            {
            case 3:  {
                return X*Y*(9*Rn000(4,p,arg) + 3*(X*X + Y*Y)*Rn000(5,p,arg) + X*X*Y*Y*Rn000(6,p,arg));
            }
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return Y*(9*Rn000(4,p,arg) + 3*(6*X2 + Y2)*Rn000(5,p,arg) + X2*((3*X2 + 6*Y2)*Rn000(6,p,arg) + X2*Y2*Rn000(7,p,arg)));
            }
            case 5:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return X*Y*(45*Rn000(5,p,arg) + 15*(2*X2 + Y2)*Rn000(6,p,arg) + X2*((3*X2 + 10*Y2)*Rn000(7,p,arg) + X2*Y2*Rn000(8,p,arg)));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=3,v=0)");
            }
            }
        }
        case 4:  {
            switch(t)
            {
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return 9*Rn000(4,p,arg) + 18*(X2 + Y2)*Rn000(5,p,arg) + (3*(X2*X2 + Y2*Y2) + 36*X2*Y2)*Rn000(6,p,arg) + 6*X2*Y2*(X2 + Y2)*Rn000(7,p,arg) + X2*X2*Y2*Y2*Rn000(8,p,arg);
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=4,v=0)");
            }
            }
        }
        default:  {
            throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with v=0)");
        }
        }
    } 
    case 1:  {
        switch(u)
        {
        case 1:  {
            switch(t)
            {
            case 1:  {
                return X*Y*Z*Rn000(3,p,arg);
            }
            case 2:  {
                return Y*Z*(Rn000(3,p,arg) + X*X*Rn000(4,p,arg));
            }
            case 3:  {
                return X*Y*Z*(3*Rn000(4,p,arg) + X*X*Rn000(5,p,arg));
            }
            case 4:  {
                return Y*Z*(3*Rn000(4,p,arg) + X*X*(6*Rn000(5,p,arg) + X*X*Rn000(6,p,arg)));
            }
            case 5:  {
                return X*Y*Z*(15*Rn000(5,p,arg) + X*X*(10*Rn000(6,p,arg) + X*X*Rn000(7,p,arg)));
            }
            case 6:  {
                double X2 = X*X;
                return Y*Z*(15*Rn000(5,p,arg) + X2*(45*Rn000(6,p,arg) + X2*(15*Rn000(7,p,arg) + X2*Rn000(8,p,arg))));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=1,v=1)");
            }
            }
        }
        case 2:  {
            switch(t)
            {
            case 2:  {
                return Z*(Rn000(3,p,arg) + (X*X + Y*Y)*Rn000(4,p,arg) + X*X*Y*Y*Rn000(5,p,arg));
            }
            case 3:  {
                return X*Z*(3*Rn000(4,p,arg) + (X*X + 3*Y*Y)*Rn000(5,p,arg) + X*X*Y*Y*Rn000(6,p,arg));
            }
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return Z*(3*Rn000(4,p,arg) + 3*(2*X2 + Y2)*Rn000(5,p,arg) + X2*((X2 + 6*Y2)*Rn000(6,p,arg) + X2*Y2*Rn000(7,p,arg)));
            }
            case 5:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return X*Z*(15*Rn000(5,p,arg) + (10*X2 + 15*Y2)*Rn000(6,p,arg) + X2*((X2 + 10*Y2)*Rn000(7,p,arg) + X2*Y2*Rn000(8,p,arg)));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=2,v=1)");
            }
            }
        }
        case 3:  {
            switch(t)
            {
            case 3:  {
                double XY = X*Y;
                return XY*Z*(9*Rn000(5,p,arg) + 3*(X*X + Y*Y)*Rn000(6,p,arg) + XY*XY*Rn000(7,p,arg));
            }
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                return Y*Z*(9*Rn000(5,p,arg) + 3*(6*X2 + Y2)*Rn000(6,p,arg) + 3*X2*((X2 + Y2)*Rn000(7,p,arg) + X2*Y2*Rn000(8,p,arg)));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=3,v=1)");
            }
            }
        }
        default:  {
            throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with v=1)");
        }
        }
    } 
    case 2:  {
        switch(u)
        {
        case 2:  {
            switch(t)
            {
            case 2:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                double Z2 = Z*Z;
                return Rn000(3,p,arg) + (X2 + Y2 + Z2)*Rn000(4,p,arg) + (X2*(Y2 + Z2) + Y2*Z2)*Rn000(5,p,arg) + X2*Y2*Z2*Rn000(6,p,arg);
            }
            case 3:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                double Z2 = Z*Z;
                return X*(3*Rn000(4,p,arg) + (X2 + 3*(Y2 + Z2))*Rn000(5,p,arg) + (X2*(Y2 + Z2) + 3*Y2*Z2)*Rn000(6,p,arg) + X2*Y2*Z2*Rn000(7,p,arg));
            }
            case 4:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                double Z2 = Z*Z;
                double Y2Z2 = Y2*Z2;
                return 3*Rn000(4,p,arg) + 3*(2*X2 + Y2 + Z2)*Rn000(5,p,arg) + (X2*(X2 + 6*(Y2 + Z2)) + 3*Y2Z2)*Rn000(6,p,arg) + X2*((6*Y2Z2 + X2*(Y2 + Z2))*Rn000(7,p,arg) + X2*Y2Z2*Rn000(8,p,arg));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=2,v=2)");
            }
            }
        }
        case 3:  {
            switch(t)
            {
            case 3:  {
                double X2 = X*X;
                double Y2 = Y*Y;
                double Z2 = Z*Z;
                return X*Y*(9*Rn000(5,p,arg) + 3*(X2 + Y2 + 3*Z2)*Rn000(6,p,arg) + (X2*Y2 + 3*Z2*(X2 + Y2))*Rn000(7,p,arg) + X2*Y2*Z2*Rn000(8,p,arg));
            }
            default: {
                throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with u=3,v=2)");
            }
            }
        }
        default:  {
            throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, with v=2)");
        }
        }
    }
    default:  {
        throw std::invalid_argument("Coulomb2Centers::HermiteCoulomb error: the R^{0}_{t,u,v} coefficients are being evaluated for t+u+v > 8 (in particular, v>=3)");
    }
    }

}


}