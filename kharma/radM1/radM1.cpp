/* 
 *  File: radM1.cpp
 *  
 *  BSD 3-Clause License
 *  
 *  Copyright (c) 2020, AFD Group at UIUC
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1. Redistributions of source code must retain the above copyright notice, this
 *     list of conditions and the following disclaimer.
 *  
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from
 *     this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "radM1.hpp"
#include "kharma_driver.hpp"
#include "units.hpp"
#include "kharma.hpp"
#include "inverter.hpp"
#include <limits> 


std::shared_ptr<KHARMAPackage> RadM1::Initialize(ParameterInput *pin, std::shared_ptr<Packages_t>& packages)
{
    bool units_enabled = pin->GetOrAddBoolean("units", "on", false);
    bool correct_connections = pin->GetOrAddBoolean("coordinates", "correct_connections", false);

    // Check if the Units package is initialized, since we need it for the radiation four-force calculations.
    if (!units_enabled) {
        printf("\033[1;31mError: Units package not enabled! It must be enabled with/BEFORE RadM1.\033[0m\n");
        exit(1);
    }
    if (!correct_connections) {
        printf("\033[1;33mError: Connection coefficient corrections for GRMHD are disabled. RadM1 requires these connections to evolve the fields properly.\033[0m\n");
        exit(1);
    }

    auto pkg = std::make_shared<KHARMAPackage>("RadM1");
    Params &params = pkg->AllParams();

    auto& driver = packages->Get("Driver")->AllParams();
    auto driver_type = driver.Get<DriverType>("type");
    bool implicit_radm1 = (driver_type == DriverType::imex && pin->GetOrAddBoolean("RadM1", "implicit", false));

    MetadataFlag areWeImplicit = (implicit_radm1) ? Metadata::GetUserFlag("Implicit")
                                              : Metadata::GetUserFlag("Explicit");

    //Registers a new label called "RADM1" for the variables to be grouped together later?
    Metadata::AddUserFlag("RADM1");
    // Properties of these new flags
    // I believe Metadata::Cell means these variables are defined at cell centers, but I should ask Cora.
    // We also add the "areWeImplicit" flag, which is either "Implicit" or "Explicit" based on the user's choice in the input file.
    std::vector<MetadataFlag> flags_radm1 = {Metadata::Real, Metadata::Cell, areWeImplicit, Metadata::GetUserFlag("RADM1")};


    //Retrieves the existing flags for the primitive and conserved variables, and adds the new radM1 flags to them.
    //Then adds the new variables for the radiation primitives and conserved variables with these flags.
    // auto flags_prim = driver.Get<std::vector<MetadataFlag>>("prim_flags");
    // flags_prim.insert(flags_prim.end(), flags_radm1.begin(), flags_radm1.end());

    auto flags_prim = driver.Get<std::vector<MetadataFlag>>("prim_flags");
    flags_prim.insert(flags_prim.end(), flags_radm1.begin(), flags_radm1.end());

    // Explicitly tag these as Primitive so your PackVariables call finds them!
    flags_prim.push_back(Metadata::GetUserFlag("Primitive")); 

    // I think this will push this metadata to restart files
    flags_prim.push_back(Metadata::Restart);

    // sync variables across boundaries (ASK CORA)
    if (pin->GetOrAddBoolean("RadM1", "sync_utop_seed", true)) { 
        flags_prim.push_back(Metadata::FillGhost);
    }
    //
    auto flags_cons = driver.Get<std::vector<MetadataFlag>>("cons_flags");
    flags_cons.insert(flags_cons.end(), flags_radm1.begin(), flags_radm1.end());

    auto m_prim_scalar = Metadata(flags_prim);
    pkg->AddField("prims.u_rad", m_prim_scalar);

    auto m_cons_scalar = Metadata(flags_cons);
    pkg->AddField("cons.u_rad", m_cons_scalar);

    auto flags_prim_vec(flags_prim);
    flags_prim_vec.push_back(Metadata::Vector);

    auto flags_cons_vec(flags_cons);
    flags_cons_vec.push_back(Metadata::Vector);

    std::vector<int> s_vector({NVEC});
    auto m_prim_vector = Metadata(flags_prim_vec, s_vector);
    pkg->AddField("prims.uvec_rad", m_prim_vector);

    auto m_cons_vector = Metadata(flags_cons_vec, s_vector);
    pkg->AddField("cons.uvec_rad", m_cons_vector);



    Real u_rad_floor = pin->GetOrAddReal("radM1", "u_rad_floor", 1.e-8);
    pkg->AllParams().Add("u_rad_floor", u_rad_floor);

    Real src_rootfind_eps = pin->GetOrAddReal("RadM1", "src_rootfind_eps", 1e-8);
    Real src_rootfind_tol = pin->GetOrAddReal("RadM1", "src_rootfind_tol", 1e-8);
    int src_rootfind_maxiter = pin->GetOrAddInteger("RadM1", "src_rootfind_maxiter", 50);
    pkg->AllParams().Add("src_rootfind_eps", src_rootfind_eps);
    pkg->AllParams().Add("src_rootfind_tol", src_rootfind_tol);
    pkg->AllParams().Add("src_rootfind_maxiter", src_rootfind_maxiter);


    // Real kappa_ep = pin->GetOrAddReal("radM1", "kappa_ep", 0.0);
    // Real sigma = pin->GetOrAddReal("radM1", "sigma", 0.0);
    // pkg->AllParams().Add("kappa_ep", kappa_ep);
    // pkg->AllParams().Add("sigma", sigma);

    //Right now, to execute the torus problem with radM1, we need to initialize the radiation primitives in fm_torus.cpp (this is stupid) (ASK Cora)
    //I think this should be moved to RadM1 method, maybe call an initialization method like a task straight after initializing the torus?
    //Especially because we'll want to initialize the radiation field for other problems. 

    // //This method should allow you to add source terms to both plasma and radiation variables separately.
    // pkg->AddSource = RadM1::AddImplicitRadiationSourceTerms;

    //Add inversion to the tasks
    pkg->BlockUtoP = RadM1::BlockUtoP;

    pkg->BlockApplyFloors = RadM1::ApplyRadM1Floors;


    return pkg;
}

void RadM1::ApplyRadM1Floors(MeshBlockData<Real> *rc, IndexDomain domain)
{
    auto pmb = rc->GetBlockPointer();
    auto& params = pmb->packages.Get("RadM1")->AllParams();

    const Real erad_floor   = params.Get<Real>("u_rad_floor");
    PackIndexMap prims_map;
    auto P = rc->PackVariables({Metadata::GetUserFlag("Primitive")}, prims_map);
    const VarMap m_p(prims_map, false);

    // We need to check if we actually have B fields enabled to avoid segfaults
    const bool has_b_field = pmb->packages.AllPackages().count("B_FluxCT") || 
                             pmb->packages.AllPackages().count("B_CD");

    auto bounds = pmb->cellbounds;
    const IndexRange ib = bounds.GetBoundsI(domain);
    const IndexRange jb = bounds.GetBoundsJ(domain);
    const IndexRange kb = bounds.GetBoundsK(domain);

    pmb->par_for("ApplyRadM1Floors", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
            
            // Absolute Floor
            Real ehat = P(m_p.UU_RAD, k, j, i);
            if (ehat < erad_floor) {
                ehat = erad_floor;
                P(m_p.UU_RAD, k, j, i) = erad_floor;
            }
        }
    );
}


TaskStatus RadM1::BlockPtoU(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    // Pack Conserved Variables (Destination)
    auto& U = rc->PackVariables(std::vector<std::string>{"cons.u_rad", "cons.uvec_rad"});

    // Pack Primitive Variables (Source)
    PackIndexMap prim_map;
    auto P = rc->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    // Get Loop Bounds
    IndexRange3 b = KDomain::GetRange(rc, domain, coarse);

    // Parallel Loop
    pmb->par_for("RadM1_PtoU", b.ks, b.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {

        Real Erf = P(m_p.UU_RAD, k, j, i);
        Real uvec_radframe[4] = {0, P(m_p.U1_RAD, k, j, i), P(m_p.U2_RAD, k, j, i), P(m_p.U3_RAD, k, j, i)};
        const Real gamma = GRMHD::lorentz_calc(G, uvec_radframe, k, j, i, Loci::center);
        Real ucon_rad[GR_DIM];
        //GRMHD::calc_ucon(G, uvec_radframe, k, j, i, Loci::center, ucon_rad);
        calc_ucon_rad(G, P, m_p, k, j, i, Loci::center, ucon_rad);
        // R^t^mu
        Real R_t_con[GR_DIM];

        for (int mu=0; mu<GR_DIM; ++mu) {
            R_t_con[mu] = 4./3. * Erf * ucon_rad[0] * ucon_rad[mu] + 1./3. * Erf * G.gcon(Loci::center, j, i, 0, mu);
        }

        // Calculat R^t_mu
        Real R_t_cov[GR_DIM];
        G.lower(R_t_con, R_t_cov, k, j, i, Loci::center);

        U(0, k, j, i) = R_t_cov[0] * G.gdet(Loci::center, j, i); // cons.u_rad
        U(1, k, j, i) = R_t_cov[1] * G.gdet(Loci::center, j, i); // cons.uvec_rad 1
        U(2, k, j, i) = R_t_cov[2] * G.gdet(Loci::center, j, i); // cons.uvec_rad 2
        U(3, k, j, i) = R_t_cov[3] * G.gdet(Loci::center, j, i); // cons.uvec_rad 3
    });
    return TaskStatus::complete;
}


KOKKOS_INLINE_FUNCTION void ApplyColdClosureFix(const GRCoordinates& G, 
                                                const Real R_t_cov_orig[GR_DIM], 
                                                const double gammarel2_fixed,
                                                const int& j, const int& i, 
                                                Real& new_R_t_t, 
                                                Real& Erf) {
    
    Real gcon_tt = G.gcon(Loci::center, j, i, 0, 0);
    
    // Time-Space cross term: g^{ti} R_i
    Real dot_t_i = G.gcon(Loci::center, j, i, 0, 1) * R_t_cov_orig[1] +
                   G.gcon(Loci::center, j, i, 0, 2) * R_t_cov_orig[2] +
                   G.gcon(Loci::center, j, i, 0, 3) * R_t_cov_orig[3];

    // Purely Spatial term: g^{ij} R_i R_j
    Real dot_spatial = G.gcon(Loci::center, j, i, 1, 1) * (R_t_cov_orig[1] * R_t_cov_orig[1]) +
                 2.0 * G.gcon(Loci::center, j, i, 1, 2) * (R_t_cov_orig[1] * R_t_cov_orig[2]) +
                       G.gcon(Loci::center, j, i, 2, 2) * (R_t_cov_orig[2] * R_t_cov_orig[2]) +
                 2.0 * G.gcon(Loci::center, j, i, 1, 3) * (R_t_cov_orig[1] * R_t_cov_orig[3]) +
                 2.0 * G.gcon(Loci::center, j, i, 2, 3) * (R_t_cov_orig[2] * R_t_cov_orig[3]) +
                       G.gcon(Loci::center, j, i, 3, 3) * (R_t_cov_orig[3] * R_t_cov_orig[3]);

    // In HARM/Kharma, alpha = 1.0 / sqrt(-gcon_tt). 
    // Therefore utsq = gammarel2 / alpha^2 simplifies to:
    Real utsq = -gammarel2_fixed * gcon_tt; 

    // The massive C-string expands exactly to (dot_t_i^2 - gcon_tt * dot_spatial).
    Real radical_inside = (dot_t_i * dot_t_i - gcon_tt * dot_spatial) 
                        * utsq * (gcon_tt + utsq) 
                        * m::pow(gcon_tt + 4.0 * utsq, 2);
                        
    Real radical = m::sqrt(m::max(0.0, radical_inside));

    // Calculate new covariant R_t (Avcov[0] in C)
    new_R_t_t = 0.25 * (-4.0 * dot_t_i * utsq * (gcon_tt + utsq) + radical) 
              / (gcon_tt * utsq * (gcon_tt + utsq));

    // Calculate fluid-frame energy Erf
    Erf = 0.75 * radical / (utsq * (gcon_tt + utsq) * (gcon_tt + 4.0 * utsq));
}

KOKKOS_INLINE_FUNCTION double CalculateGammaRel2(const GRCoordinates& G, const Real R_t_cov[GR_DIM], const Real invariant_scalar, const int& j, const int& i) {
    Real gcon_tt = G.gcon(Loci::center, j, i, 0, 0);
    Real R_t_t = R_t_cov[0];
    
    // Time-Space cross term: g^{ti} R^t_i
    Real dot_t_i = G.gcon(Loci::center, j, i, 0, 1) * R_t_cov[1] +
                   G.gcon(Loci::center, j, i, 0, 2) * R_t_cov[2] +
                   G.gcon(Loci::center, j, i, 0, 3) * R_t_cov[3];

    // Purely Spatial term: g^{ij} R^t_i R^t_j
    Real dot_spatial = G.gcon(Loci::center, j, i, 1, 1) * (R_t_cov[1] * R_t_cov[1]) +
                   2.0 * G.gcon(Loci::center, j, i, 1, 2) * (R_t_cov[1] * R_t_cov[2]) +
                         G.gcon(Loci::center, j, i, 2, 2) * (R_t_cov[2] * R_t_cov[2]) +
                   2.0 * G.gcon(Loci::center, j, i, 1, 3) * (R_t_cov[1] * R_t_cov[3]) +
                   2.0 * G.gcon(Loci::center, j, i, 2, 3) * (R_t_cov[2] * R_t_cov[3]) +
                         G.gcon(Loci::center, j, i, 3, 3) * (R_t_cov[3] * R_t_cov[3]);

    // Calculate Roots for (u^t_R)^2
    Real radical_inside = 4.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                        + (dot_t_i * dot_t_i)
                        + gcon_tt * (8.0 * R_t_t * dot_t_i + 3.0 * dot_spatial);
    Real radical = m::sqrt(m::max(0.0, radical_inside)); 

    Real num_a = 2.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                + dot_t_i * (dot_t_i + radical)
                + gcon_tt * (4.0 * R_t_t * dot_t_i + dot_spatial + R_t_t * radical);

    Real num_b = -2.0 * (gcon_tt * gcon_tt) * (R_t_t * R_t_t)
                - gcon_tt * (4.0 * R_t_t * dot_t_i + dot_spatial)
                + gcon_tt * R_t_t * radical
                + dot_t_i * (-dot_t_i + radical);

    Real gamma2a = -0.25 * num_a / invariant_scalar;
    Real gamma2b =  0.25 * num_b / invariant_scalar;

    // Sadowski 2013: Select the largest real root for (u^t)^2
    Real gamma2 = gamma2a;
    if (m::isnan(gamma2a) || m::isinf(gamma2a)) {
        gamma2 = gamma2b;
    } else if (!m::isnan(gamma2b) && !m::isinf(gamma2b)) {
        gamma2 = m::max(gamma2a, gamma2b);
    }

    Real alpha_sq = -1.0 / gcon_tt;
    Real gammarel2 = gamma2 * alpha_sq; 

    // Hard floor for physical bounds (Lorentz factor squared MUST be >= 1.0)
    // TODO: FLOOR! CHANGE THIS
    Real GAMMA_SMALL_LIMIT = (1.0-1e-10);
    if (gammarel2 < 1.0 && gammarel2 > GAMMA_SMALL_LIMIT) {
        gammarel2 = 1.0;
    }
    
    return gammarel2;
}

KOKKOS_INLINE_FUNCTION TaskStatus CalculateRadPrimitive_M1(
    const GRCoordinates& G,
    const Real U_rad[4],
    Real P_rad[4],
    const int k, const int j, const int i) 
{
    Real gdet = G.gdet(Loci::center, j, i);

    //Recover R^t_mu from conserved state
    Real R_t_cov[4] = {
        U_rad[0] / gdet,
        U_rad[1] / gdet,
        U_rad[2] / gdet,
        U_rad[3] / gdet
    };

    //Raise index to get R^{t\mu}
    Real R_t_con[4];
    G.raise(R_t_cov, R_t_con, k, j, i, Loci::center);

    //Calculate Invariant Scalar S = R^t_\mu * R^{t\mu}
    Real invariant_scalar = 0.0;
    for(int mu = 0; mu < 4; ++mu) {
        invariant_scalar += R_t_cov[mu] * R_t_con[mu];
    }

    // Calculate gamma^2 for the radiation frame
    Real gammarel2 = CalculateGammaRel2(G, R_t_cov, invariant_scalar, j, i);

    //Pre-calculate alpha bounds and rest-frame energy E_rf
    Real alpha_sq = -1.0 / G.gcon(Loci::center, j, i, 0, 0); 
    Real alpha = m::sqrt(alpha_sq);
    Real E_rf = (3.0 * R_t_con[0] * alpha_sq) / (4.0 * gammarel2 - 1.0);
    
    // Limits
    const Real min_erad = 10. * 1.e-9;
    const Real GAMMAMAX = 20.0;
    int nonfailure = (gammarel2 >= 1.0) && (E_rf > min_erad) && 
                     (gammarel2 <= (GAMMAMAX * GAMMAMAX) / (min_erad * min_erad));
    
    Real uvec_radframe_con[4] = {0};

    // Evaluate valid primitives or apply cold closure fix
    if (nonfailure) {
        for(int mu = 0; mu < 4; ++mu) {
            uvec_radframe_con[mu] = alpha * (R_t_con[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2 - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2));
        }
    } else {
        // Attempt Cold Closure
        Real gammarel2_slow = m::pow(1.0 + 10.0 * std::numeric_limits<double>::epsilon(), 2.0);
        Real gammarel2_fast = GAMMAMAX * GAMMAMAX;

        Real R_t_t_slow, Erf_slow;
        ApplyColdClosureFix(G, R_t_cov, gammarel2_slow, j, i, R_t_t_slow, Erf_slow);

        Real R_t_t_fast, Erf_fast;
        ApplyColdClosureFix(G, R_t_cov, gammarel2_fast, j, i, R_t_t_fast, Erf_fast);

        Real R_t_t_new, gammarel2_new;
        
        if (m::abs(R_t_t_slow - R_t_cov[0]) > m::abs(R_t_t_fast - R_t_cov[0])) {
            R_t_t_new = R_t_t_fast;
            E_rf = Erf_fast;
            gammarel2_new = gammarel2_fast;
        } else {
            R_t_t_new = R_t_t_slow;
            E_rf = Erf_slow;
            gammarel2_new = gammarel2_slow;
        }

        // If even the closure fix yields a non-positive energy, this step is a failure. I'm not setting failure here
        // cause otherwise it won't run...
        if (E_rf <= 0.0) {
            P_rad[0] = min_erad;
            P_rad[1] = 0.0;
            P_rad[2] = 0.0;
            P_rad[3] = 0.0;
            return TaskStatus::complete;
        }

        Real R_t_cov_new[4] = {R_t_t_new, R_t_cov[1], R_t_cov[2], R_t_cov[3]};
        Real R_t_con_new[4];
        G.raise(R_t_cov_new, R_t_con_new, k, j, i, Loci::center);

        if (E_rf > 0.0) {
            for(int mu = 0; mu < 4; ++mu) {
                uvec_radframe_con[mu] = alpha * (R_t_con_new[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2_new - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2_new));
            }
        } else {
            for(int mu = 0; mu < 4; ++mu) uvec_radframe_con[mu] = 0.0;
        }
    }

    P_rad[0] = E_rf;
    P_rad[1] = uvec_radframe_con[1];
    P_rad[2] = uvec_radframe_con[2];
    P_rad[3] = uvec_radframe_con[3];

    return TaskStatus::complete;
}


TaskStatus RadM1::BlockUtoP(MeshBlockData<Real> *rc, IndexDomain domain, bool coarse)
{
    auto pmb = rc->GetBlockPointer();
    const auto& G = pmb->coords;

    //Pack Variables
    auto U = rc->PackVariables(std::vector<std::string>{"cons.u_rad", "cons.uvec_rad"});

    // Pack Primitive Variables
    PackIndexMap prim_map;
    auto P = rc->PackVariables(std::vector<MetadataFlag>{Metadata::GetUserFlag("Primitive")}, prim_map);
    const VarMap m_p(prim_map, false);

    //Get Loop Bounds
    IndexRange3 b = KDomain::GetRange(rc, domain, coarse);

    auto& params = pmb->packages.Get("RadM1")->AllParams();
    //TODO: FLOOR! CHANGE THIS
    const Real min_erad   = 10.*1.e-80;

    //Parallel Loop
    pmb->par_for("RadM1_UtoP", b.ks, b.ke, b.js, b.je, b.is, b.ie,
        KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {


            
            Real gdet = G.gdet(Loci::center, j, i);

            // Reconstruct Radiation Stress-Energy Tensor R^{mu, nu} 
            // For the M1 scheme, we assume the radiation is isotropic and satisfies the Eddington approximation (P^{ij} = (1/3) E delta^{ij} in the fluid frame)
            // but in the radiation frame. Therefore, \bar{R}^{tt} = E, \bar{R}^{ii} = \bar{E}/3, and every other component is zero. 
            // So we have Equation 27 in Sadowski et al. 2013 in the radiation rest frame, but since it's covariant, it's valid for every other frame (including lab frame).
            // The equation goes as follows:
            // R^{mu nu} = (4/3) \bar{E} u_R^mu u_R^nu + (1/3) \bar{E} g^{mu nu}, \bar{E} is always in the radiation rest frame.


            //Recover R^t_mu
            // U(0) is R^t_t, U(1..3) are R^t_i
            Real R_t_cov[GR_DIM] = {
                U(0, k, j, i) / gdet,
                U(1, k, j, i) / gdet,
                U(2, k, j, i) / gdet,
                U(3, k, j, i) / gdet
            };

            // Get R^{t\mu} from R^t_\mu by doing R^{t\mu} = g^{t\nu} R_{\nu\mu}
            // This is R^t_\mu
            Real R_t_con[GR_DIM];
            G.raise(R_t_cov, R_t_con, k, j, i, Loci::center);

            //Calculate Invariant Scalar S = R^t_mu * R^{t\mu}
            // Real invariant_scalar = 0.0;
            // for(int mu=0; mu<4; ++mu) {
            //      // Note: R_t_cov is R^t_mu (mixed), R_t_con is R^{t\mu} (upper). 
            //      // S = g_{mu nu} R^{t mu} R^{t nu} 
            //      for(int nu=0; nu<4; ++nu) {
            //         invariant_scalar += G.gcov(Loci::center, j, i, mu, nu) * R_t_con[mu] * R_t_con[nu];
            //      }
            // }
            Real invariant_scalar = 0.0;
            for(int mu=0; mu<4; ++mu) {
                invariant_scalar += R_t_cov[mu] * R_t_con[mu];
            }

            // From Sadowski et al. Eq 32 and 33
            // Solving by isolating Isolaring u^t_R^2 results in a much better equation, but I think (and I could be wrong)
            // that solving by isolating \bar{E} results in a much more numerical stable solution.

            // Calculating u^t_R^2
            Real gammarel2 = CalculateGammaRel2(G, R_t_cov, invariant_scalar, j, i);

            // Pre-calculate alpha_sq so E_rf can use it
            Real alpha_sq = -1.0 / G.gcon(Loci::center, j, i, 0, 0); 
            Real alpha = m::sqrt(alpha_sq);
            Real E_rf = (3.0 * R_t_con[0] * alpha_sq) / (4.0 * gammarel2 - 1.0);
            
            // TODO: FLOOR! CHANGE THIS
            Real GAMMAMAX = 20.;
            int nonfailure = gammarel2 >= 1.0 && E_rf > min_erad && gammarel2 <= (GAMMAMAX * GAMMAMAX)/(min_erad * min_erad);
            Real uvec_radframe_con[4] = {0};

            if(nonfailure){
                for(int mu=0; mu<4; ++mu) {
                    uvec_radframe_con[mu] = alpha * (R_t_con[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2 - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2));
                }
            } else {
                // Attempt Cold Closure
                Real gammarel2_slow = m::pow(1.0 +10.0 * std::numeric_limits<double>::epsilon(), 2.0);
                Real gammarel2_fast = GAMMAMAX * GAMMAMAX;

                Real R_t_t_slow, Erf_slow;
                ApplyColdClosureFix(G, R_t_cov, gammarel2_slow, j, i, R_t_t_slow, Erf_slow);

                Real R_t_t_fast, Erf_fast;
                ApplyColdClosureFix(G, R_t_cov, gammarel2_fast, j, i, R_t_t_fast, Erf_fast);

                Real R_t_t_new, gammarel2_new;
                
                if (m::abs(R_t_t_slow - R_t_cov[0]) > m::abs(R_t_t_fast - R_t_cov[0])) {
                    R_t_t_new = R_t_t_fast;
                    E_rf = Erf_fast;
                    gammarel2_new = gammarel2_fast;
                } else {
                    R_t_t_new = R_t_t_slow;
                    E_rf = Erf_slow;
                    gammarel2_new = gammarel2_slow;
                }

                Real R_t_cov_new[GR_DIM] = {R_t_t_new, R_t_cov[1], R_t_cov[2], R_t_cov[3]};

                Real R_t_con_new[GR_DIM];
                G.raise(R_t_cov_new, R_t_con_new, k, j, i, Loci::center);

                // Calculate primitive velocities with the new state
                if (E_rf > 0.0) {
                    for(int mu=0; mu<4; ++mu) {
                        uvec_radframe_con[mu] = alpha * (R_t_con_new[mu] + 1./3. * E_rf * G.gcon(Loci::center, j, i, 0, mu) * (4.0 * gammarel2_new - 1.0)) / (4./3. * E_rf * m::sqrt(gammarel2_new));
                    }
                } else {
                    for(int mu=0; mu<4; ++mu) uvec_radframe_con[mu] = 0.0;
                }
            }

            P(m_p.UU_RAD, k, j, i) = E_rf;
            P(m_p.U1_RAD, k, j, i) = uvec_radframe_con[1];
            P(m_p.U2_RAD, k, j, i) = uvec_radframe_con[2];
            P(m_p.U3_RAD, k, j, i) = uvec_radframe_con[3];

    });
    return TaskStatus::complete;
}


#define RAD_LARGE (0.1 * std::numeric_limits<Real>::max())
#define RAD_SMALL (10.0 * std::numeric_limits<Real>::min())
#define RAD_EPS   (10.0 * std::numeric_limits<Real>::epsilon())


KOKKOS_INLINE_FUNCTION
void adjoint(Real m[16], Real adjOut[16]) {
  adjOut[0] = MINOR(m, 1, 2, 3, 1, 2, 3);
  adjOut[1] = -MINOR(m, 0, 2, 3, 1, 2, 3);
  adjOut[2] = MINOR(m, 0, 1, 3, 1, 2, 3);
  adjOut[3] = -MINOR(m, 0, 1, 2, 1, 2, 3);

  adjOut[4] = -MINOR(m, 1, 2, 3, 0, 2, 3);
  adjOut[5] = MINOR(m, 0, 2, 3, 0, 2, 3);
  adjOut[6] = -MINOR(m, 0, 1, 3, 0, 2, 3);
  adjOut[7] = MINOR(m, 0, 1, 2, 0, 2, 3);

  adjOut[8] = MINOR(m, 1, 2, 3, 0, 1, 3);
  adjOut[9] = -MINOR(m, 0, 2, 3, 0, 1, 3);
  adjOut[10] = MINOR(m, 0, 1, 3, 0, 1, 3);
  adjOut[11] = -MINOR(m, 0, 1, 2, 0, 1, 3);

  adjOut[12] = -MINOR(m, 1, 2, 3, 0, 1, 2);
  adjOut[13] = MINOR(m, 0, 2, 3, 0, 1, 2);
  adjOut[14] = -MINOR(m, 0, 1, 3, 0, 1, 2);
  adjOut[15] = MINOR(m, 0, 1, 2, 0, 1, 2);
}

KOKKOS_INLINE_FUNCTION
Real MINOR(Real m[16], int r0, int r1, int r2, int c0, int c1, int c2) {
  return m[4 * r0 + c0] *
             (m[4 * r1 + c1] * m[4 * r2 + c2] - m[4 * r2 + c1] * m[4 * r1 + c2]) -
         m[4 * r0 + c1] *
             (m[4 * r1 + c0] * m[4 * r2 + c2] - m[4 * r2 + c0] * m[4 * r1 + c2]) +
         m[4 * r0 + c2] *
             (m[4 * r1 + c0] * m[4 * r2 + c1] - m[4 * r2 + c0] * m[4 * r1 + c1]);
}

KOKKOS_INLINE_FUNCTION
Real determinant(Real m[16]) {
  return m[0] * MINOR(m, 1, 2, 3, 1, 2, 3) - m[1] * MINOR(m, 1, 2, 3, 0, 2, 3) +
         m[2] * MINOR(m, 1, 2, 3, 0, 1, 3) - m[3] * MINOR(m, 1, 2, 3, 0, 1, 2);
}

KOKKOS_INLINE_FUNCTION
void matrixInverse4x4(Real mat[4][4], Real matinv[4][4]) {
  adjoint(&(mat[0][0]), &(matinv[0][0]));

  Real det = determinant(&(mat[0][0]));
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      matinv[i][j] /= det;
    }
  }
}



KOKKOS_INLINE_FUNCTION void ComputeCovariantFourForce(
    const GRCoordinates& G, 
    const Real P_mhd[4], 
    const Real P_rad[4], 
    const Real Gas_Rho,  
    const Real gam,      
    const int k, const int j, const int i, 
    Real dS[4]) 
{
    Real uvec[3] = {P_mhd[1], P_mhd[2], P_mhd[3]};
    Real ucon[4];
    
    GRMHD::calc_ucon(G, uvec, k, j, i, Loci::center, ucon);

    Real ucov[4];
    G.lower(ucon, ucov, k, j, i, Loci::center);


    Real Tg = (gam - 1.0) * (P_mhd[0] / Gas_Rho); 
    
    const Real sigma_rad = 3.085e9; // Test 1 value
    const Real kappa_rho = 0.4;     // Test 1 value
    
    Real JBB = sigma_rad * (Tg * Tg * Tg * Tg); // E_eq = sigma * T^4
    Real kappa_a = Gas_Rho * kappa_rho;
    Real kappa_tot = kappa_a; 
    
    // Safety cap for extremely optically thick regimes (Phoebus use the same thing)
    kappa_a = m::min(kappa_a, 1.e5);
    kappa_tot = m::min(kappa_tot, 1.e5);

    //u_mu * F^mu = 0 to find F_0
    Real F_hat_cov[4] = {0.0, P_rad[1], P_rad[2], P_rad[3]};
    Real F_hat_con[4];
    G.raise(F_hat_cov, F_hat_con, k, j, i, Loci::center);
    
    // F_hat_0 = - (u_i * F_hat^i) / u^0
    Real uF_space = ucov[1]*F_hat_con[1] + ucov[2]*F_hat_con[2] + ucov[3]*F_hat_con[3];
    F_hat_cov[0] = -uF_space / ucon[0];

    // G_mu = kappa_a * (E_eq - E_hat) * u_mu - kappa_tot * F_hat_mu
    Real coupling_term = kappa_a * (JBB - P_rad[0]);
    
    dS[0] = coupling_term * ucov[0] - kappa_tot * F_hat_cov[0]; // Energy exchange
    dS[1] = coupling_term * ucov[1] - kappa_tot * F_hat_cov[1]; // Momentum exchange X1
    dS[2] = coupling_term * ucov[2] - kappa_tot * F_hat_cov[2]; // Momentum exchange X2
    dS[3] = coupling_term * ucov[3] - kappa_tot * F_hat_cov[3]; // Momentum exchange X3
    // dS[0] = 0.0;
    // dS[1] = 0.0;
    // dS[2] = 0.0;
    // dS[3] = 0.0;
}


KOKKOS_INLINE_FUNCTION TaskStatus SolveImplicitRadiation4D(
    const GRCoordinates& G,
    const VariablePack<Real> U_init, 
    const VariablePack<Real> P_init, 
    VariablePack<Real> P_new,       
    VariablePack<Real> U_new,       
    const VarMap m_p,               
    const VarMap m_u,               
    const int k, const int j, const int i,
    const Real dt, const Real gam, const double src_rootfind_eps, const double src_rootfind_tol, const int src_rootfind_maxiter)
{
    //extract initial fluid and magnetic field primitives from the grid (P_init) to generate the initial guess for U_mhd and U_rad.
    const Real Gas_Rho = P_init(m_p.RHO, k, j, i);
    Real B_P[3]  = {P_init(m_p.B1, k, j, i), 
                          P_init(m_p.B2, k, j, i), 
                          P_init(m_p.B3, k, j, i)};

    
    Real P_mhd_guess[4] = {
        P_init(m_p.UU, k, j, i), // Index 0 is internal energy
        P_init(m_p.U1, k, j, i), 
        P_init(m_p.U2, k, j, i), 
        P_init(m_p.U3, k, j, i)
    };

    Real U_rad_0[4] = {
        U_init(m_u.UU_RAD, k, j, i),
        U_init(m_u.U1_RAD, k, j, i),
        U_init(m_u.U2_RAD, k, j, i),
        U_init(m_u.U3_RAD, k, j, i)
    };

    Real resid[4];
    bool success = false;

    Real ug_init = P_init(m_p.UU, k, j, i);

    Real U_mhd_0[4]; 
    Real rho_ut_dummy; 

    Real uvec[NVEC] = {P_init(m_p.U1, k, j, i), P_init(m_p.U2, k, j, i), P_init(m_p.U3, k, j, i)};
    
    GRMHD::p_to_u_mhd(G, Gas_Rho, ug_init, uvec, B_P, gam, k, j, i, 
               rho_ut_dummy, U_mhd_0, Loci::center);

    //Extract Radition conserved variables
    // Real R_init[4][4];
    // for (int dir = 0; dir < 4; ++dir) {
    //     calc_tensor_m1(G, P_init, m_p, dir, k, j, i, Loci::center, R_init[dir]);
    // }

    Real P_rad_m[4];
    Real P_rad_p[4];
    Real U_mhd_m[4];
    Real U_mhd_p[4];
    Real U_rad_m[4];
    Real U_rad_p[4];
    Real dS_m[4];
    Real dS_p[4];

    Real U_mhd_guess[4];
    Real U_rad_guess[4];
    Real P_rad_guess[4];
    Real dS_guess[4];
    Real dcov_rad[4];

    //Iteration 0
    
    // Convert current gas primitives guess to conserved variables
    uvec[0] = P_mhd_guess[1];
    uvec[1] = P_mhd_guess[2];
    uvec[2] = P_mhd_guess[3];
    GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_guess[0], uvec, B_P, gam, k, j, i, 
               rho_ut_dummy, U_mhd_guess, Loci::center);

    // Conservation law: Delta U_rad = - Delta U_mhd
    for (int n = 0; n < 4; n++) {
        U_rad_guess[n] = U_rad_0[n] - (U_mhd_guess[n] - U_mhd_0[n]);
    }
    Real gdet = G.gdet(Loci::center, j, i);

    // Convert the newly guessed U_rad to P_rad
    auto status = CalculateRadPrimitive_M1(G, U_rad_guess, P_rad_guess, k, j, i);
    ComputeCovariantFourForce(G, P_mhd_guess, P_rad_guess, Gas_Rho, gam, k, j, i, dS_guess);

    for (int n = 0; n < 4; n++) dS_guess[n] = gdet * dS_guess[n];

    for (int n = 0; n < 4; n++) {
        resid[n] = U_mhd_guess[n] - U_mhd_0[n] + dt * dS_guess[n];
    }

    Real err = 1.e100;
    int niter = 0;
    bool bad_guess = false;
    
    do {
        Real jac[4][4] = {0};

        // Find minimum non-zero magnitude from P_mhd_guess to scale FD step safely
        Real P_mhd_mag_min = RAD_LARGE; // Or std::numeric_limits<Real>::max()
        for (int m = 0; m < 4; m++) {
            if (m::abs(P_mhd_guess[m]) > 0.) {
                P_mhd_mag_min = m::min(P_mhd_mag_min, m::abs(P_mhd_guess[m]));
            }
        }

        bool bad_guess_m = false;
        bool bad_guess_p = false;

        // Loop over the 4 fluid variables to perturb each one
        for (int m = 0; m < 4; m++) {
            Real P_mhd_m[4] = {P_mhd_guess[0], P_mhd_guess[1], P_mhd_guess[2], P_mhd_guess[3]};
            Real P_mhd_p[4] = {P_mhd_guess[0], P_mhd_guess[1], P_mhd_guess[2], P_mhd_guess[3]};
            
            const Real fd_step = m::max(src_rootfind_eps * P_mhd_mag_min,
                                        src_rootfind_eps * m::abs(P_mhd_guess[m]));
            P_mhd_m[m] -= fd_step;
            P_mhd_p[m] += fd_step;

            // Evaluate minus perturbation
            Real rho_ut_dummy_m;
            uvec[0] = P_mhd_m[1];
            uvec[1] = P_mhd_m[2];
            uvec[2] = P_mhd_m[3];
            GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_m[0], uvec, B_P, gam, k, j, i, 
                        rho_ut_dummy_m, U_mhd_m, Loci::center);

            // Conservation law: Delta U_rad = - Delta U_mhd
            for (int n = 0; n < 4; n++) {
                U_rad_m[n] = U_rad_0[n] - (U_mhd_m[n] - U_mhd_0[n]);
            }

            // Recover rad primitives
            auto status_m = CalculateRadPrimitive_M1(G, U_rad_m, P_rad_m, k, j, i);
            if (status_m == TaskStatus::fail) {
                bad_guess_m = true;
            }

            
            ComputeCovariantFourForce(G, P_mhd_m, P_rad_m, Gas_Rho, gam, k, j, i, dS_m);
            for (int n = 0; n < 4; n++) dS_m[n] = gdet * dS_m[n];


            // Evaluate plus perturbation
            Real rho_ut_dummy_p;
            uvec[0] = P_mhd_p[1];
            uvec[1] = P_mhd_p[2];
            uvec[2] = P_mhd_p[3];
            GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_p[0], uvec, B_P, gam, k, j, i, 
                        rho_ut_dummy_p, U_mhd_p, Loci::center);

            for (int n = 0; n < 4; n++) {
                U_rad_p[n] = U_rad_0[n] - (U_mhd_p[n] - U_mhd_0[n]);
            }

            auto status_p = CalculateRadPrimitive_M1(G, U_rad_p, P_rad_p, k, j, i);
            if (status_p == TaskStatus::fail) {
                bad_guess_p = true;
            }            
            ComputeCovariantFourForce(G, P_mhd_p, P_rad_p, Gas_Rho, gam, k, j, i, dS_p);
            for (int n = 0; n < 4; n++) dS_p[n] = gdet * dS_p[n];


            //Populate Jacobian
            for (int n = 0; n < 4; n++) {
                Real fp = U_mhd_p[n] - U_mhd_0[n] + dt * dS_p[n];
                Real fm = U_mhd_m[n] - U_mhd_0[n] + dt * dS_m[n];
                jac[n][m] = (fp - fm) / (P_mhd_p[m] - P_mhd_m[m]);
            }
        }

        if(bad_guess_m == true && bad_guess_p == true) {
            bad_guess = true;
            break; // Exit the iteration loop if both perturbations yield bad guesses
        } else if(bad_guess_m == true){
            // If only - finite difference support point is bad, do one-sided
            // difference with + support point

            for (int m = 0; m < 4; m++) {
                Real P_mhd_p[4] = {P_mhd_guess[0], P_mhd_guess[1], P_mhd_guess[2], P_mhd_guess[3]};
                P_mhd_p[m] += std::max(src_rootfind_eps * P_mhd_mag_min, src_rootfind_eps * m::abs(P_mhd_p[m]));

                Real rho_ut_dummy_p;
                uvec[0] = P_mhd_p[1];
                uvec[1] = P_mhd_p[2];
                uvec[2] = P_mhd_p[3];
                GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_p[0], uvec, B_P, gam, k, j, i, 
                            rho_ut_dummy_p, U_mhd_p, Loci::center);
                for (int n = 0; n < 4; n++) {
                    U_rad_p[n] = U_rad_0[n] - (U_mhd_p[n] - U_mhd_0[n]);
                }
                auto status_p = CalculateRadPrimitive_M1(G, U_rad_p, P_rad_p, k, j, i);
                ComputeCovariantFourForce(G, P_mhd_p, P_rad_p, Gas_Rho, gam, k, j, i, dS_p);
                for (int n = 0; n < 4; n++) dS_p[n] = gdet * dS_p[n];
                
                PARTHENON_DEBUG_REQUIRE(status_p == TaskStatus::complete,
                                    "This inversion should have already worked!");
                for (int n = 0; n < 4; n++) {
                    Real fp = U_mhd_p[n] - U_mhd_0[n] + dt * dS_p[n];
                    Real fguess = U_mhd_guess[n] - U_mhd_0[n] + dt * dS_guess[n];
                    jac[n][m] = (fp - fguess) / (P_mhd_p[m] - P_mhd_guess[m]);
                }
            }
        } else if(bad_guess_p == true) {
            // If only + finite difference support point is bad, do one-sided
            // difference with - support point

            for (int m = 0; m < 4; m++) {
                Real P_mhd_m[4] = {P_mhd_guess[0], P_mhd_guess[1], P_mhd_guess[2], P_mhd_guess[3]};
                P_mhd_m[m] -= std::max(src_rootfind_eps * P_mhd_mag_min, src_rootfind_eps * m::abs(P_mhd_m[m]));

                Real rho_ut_dummy_m;
                uvec[0] = P_mhd_m[1];
                uvec[1] = P_mhd_m[2];
                uvec[2] = P_mhd_m[3];
                GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_m[0], uvec, B_P, gam, k, j, i, 
                            rho_ut_dummy_m, U_mhd_m, Loci::center);
                for (int n = 0; n < 4; n++) {
                    U_rad_m[n] = U_rad_0[n] - (U_mhd_m[n] - U_mhd_0[n]);
                }
                auto status_m = CalculateRadPrimitive_M1(G, U_rad_m, P_rad_m, k, j, i);
                ComputeCovariantFourForce(G, P_mhd_m, P_rad_m, Gas_Rho, gam, k, j, i, dS_m);
                for (int n = 0; n < 4; n++) dS_m[n] = gdet * dS_m[n];
                PARTHENON_DEBUG_REQUIRE(status_m == TaskStatus::complete,
                                    "This inversion should have already worked!");
                for (int n = 0; n < 4; n++) {
                    Real fm = U_mhd_m[n] - U_mhd_0[n] + dt * dS_m[n];
                    Real fguess = U_mhd_guess[n] - U_mhd_0[n] + dt * dS_guess[n];
                    jac[n][m] = (fguess - fm) / (P_mhd_guess[m] - P_mhd_m[m]);
                }
            }
        }
    
        Real jacinv[4][4];
        matrixInverse4x4(jac, jacinv);

        if(bad_guess == false){
            Real ug0 = P_mhd_guess[0];
            Real ur0 = P_rad_guess[0];

            //update guess
            for(int m = 0; m < 4; m++) {
                for(int n = 0; n < 4; n++) {
                    P_mhd_guess[m] -= jacinv[m][n] * resid[n];
                }
            }

            // Re-evaluate residual with updated guess
            Real rho_ut_dummy_new;
            uvec[0] = P_mhd_guess[1];
            uvec[1] = P_mhd_guess[2];
            uvec[2] = P_mhd_guess[3];
            GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_guess[0], uvec, B_P, gam, k, j, i, 
                        rho_ut_dummy_new, U_mhd_guess, Loci::center);
            for (int n = 0; n < 4; n++) {
                U_rad_guess[n] = U_rad_0[n] - (U_mhd_guess[n] - U_mhd_0[n]);
            }
            auto status = CalculateRadPrimitive_M1(G, U_rad_guess, P_rad_guess, k, j, i);
            ComputeCovariantFourForce(G, P_mhd_guess, P_rad_guess, Gas_Rho, gam, k, j, i, dS_guess);

            for (int n = 0; n < 4; n++) dS_guess[n] = gdet * dS_guess[n];

            if(status == TaskStatus::fail){
                bad_guess = true;
            }

            if(bad_guess){
                // Try reducing the step size.
                constexpr Real umin = 1.e-12;
                constexpr Real Emin = 1.e-60;
                const Real gamma_max_sq = 400.0; // Corresponds to Gamma_max = 20

                Real scaling_factor = 0.0;

                // Check Gas Internal Energy Violation
                if (P_mhd_guess[0] < umin) {
                    // If it went negative or too small, calculate relative overstep
                    scaling_factor = m::max(scaling_factor, (ug0 - umin) / (ug0 - P_mhd_guess[0] + 1e-20));
                }

                //Check Radiation Rest-Frame Energy Violation
                if (P_rad_guess[0] < Emin) {
                    scaling_factor = m::max(scaling_factor, (ur0 - Emin) / (ur0 - P_rad_guess[0] + 1e-20));
                }

                // Check Velocity / Lorentz Factor Violation
                // Re-calculate the Lorentz factor squared of your guess to see if it went superluminal
                Real invariant_scalar_guess = 0.0;
                Real R_t_cov_guess[4] = {U_rad_guess[0]/gdet, U_rad_guess[1]/gdet, U_rad_guess[2]/gdet, U_rad_guess[3]/gdet};
                Real R_t_con_guess[4];
                G.raise(R_t_cov_guess, R_t_con_guess, k, j, i, Loci::center);
                for(int mu = 0; mu < 4; ++mu) invariant_scalar_guess += R_t_cov_guess[mu] * R_t_con_guess[mu];
                
                Real gamma_sq_guess = CalculateGammaRel2(G, R_t_cov_guess, invariant_scalar_guess, j, i);

                if (gamma_sq_guess > gamma_max_sq || gamma_sq_guess < 1.0) {
                    // If velocity exploded, aggressively damp the step (e.g., cut it in half)
                    scaling_factor = m::max(scaling_factor, 0.5);
                }

                // Verify the scaling factor is sane
                if (!(scaling_factor > 0.0 && scaling_factor <= 1.0)) {
                    bad_guess = true;
                    break; // Step is completely unrecoverable, abort to avoid NaN cascading
                }

                //Retain a 50% safety buffer away from the boundary edge
                scaling_factor *= 0.5;

                // Roll back to old guess, then take the scaled/damped step
                for (int m = 0; m < 4; m++) {
                    for (int n = 0; n < 4; n++) {
                        // Notice the += here!
                        P_mhd_guess[m] += (1.0 - scaling_factor) * jacinv[m][n] * resid[n];
                    }
                }

                // Re-evaluate unperturbed state with the newly dampened primitives
                uvec[0] = P_mhd_guess[1];
                uvec[1] = P_mhd_guess[2];
                uvec[2] = P_mhd_guess[3];
                GRMHD::p_to_u_mhd(G, Gas_Rho, P_mhd_guess[0], uvec, B_P, gam, k, j, i, rho_ut_dummy_new, U_mhd_guess, Loci::center);
                for (int n = 0; n < 4; n++) U_rad_guess[n] = U_rad_0[n] - (U_mhd_guess[n] - U_mhd_0[n]);
                
                status = CalculateRadPrimitive_M1(G, U_rad_guess, P_rad_guess, k, j, i);
                ComputeCovariantFourForce(G, P_mhd_guess, P_rad_guess, Gas_Rho, gam, k, j, i, dS_guess);
                for (int n = 0; n < 4; n++) dS_guess[n] = gdet * dS_guess[n];

                // If the scaled step lands in a physically valid regime, we cleared the error flag!
                if (status != TaskStatus::fail && P_mhd_guess[0] > umin) {
                    bad_guess = false; 
                }
            }
        }

        // if guess was invalid, break and mark this zone as a failure
        if (bad_guess){
            break;
        }

        //Update residuals
        for (int n = 0; n < 4; n++) {
            resid[n] = U_mhd_guess[n] - U_mhd_0[n] + dt * dS_guess[n];
        
            if (std::isnan(resid[n])){
                bad_guess = true;
                break;
            }
        }

        if (bad_guess == true){
            break;
        }

        //Calculate error now

        err = RAD_SMALL;
        Real max_divisor = RAD_SMALL;

        for (int n = 0; n < 4; n++) {
            max_divisor = std::max<Real>(max_divisor, std::fabs(U_mhd_guess[n]) +
                                        std::fabs(U_mhd_0[n]) +
                                        std::fabs(dt * dS_guess[n]));

        }

        for (int n = 0; n < 4; n++){
            Real suberr = std::fabs(resid[n]) / max_divisor;
            if(suberr > err) {
                err = suberr;
            }
        }

        niter++;
    } while (err > src_rootfind_tol && niter < src_rootfind_maxiter);

    if (niter == src_rootfind_maxiter || err > src_rootfind_tol ||
    std::isnan(U_rad_guess[0]) || std::isnan(U_rad_guess[1]) || 
    std::isnan(U_rad_guess[2]) || std::isnan(U_rad_guess[3]) || 
    bad_guess) {
        success = false;
    }else{
        success = true;

        dcov_rad[0] = U_rad_guess[0] - U_rad_0[0]; // Energy delta
        dcov_rad[1] = U_rad_guess[1] - U_rad_0[1]; // Momentum 1 delta
        dcov_rad[2] = U_rad_guess[2] - U_rad_0[2]; // Momentum 2 delta
        dcov_rad[3] = U_rad_guess[3] - U_rad_0[3]; // Momentum 3 delta
    }

    if (success == false){
        //TODO: Fall back to 1d solver, if 4d encounters an issue
        // For now, I'll only implement the 4D solver, so if it fails just return fail
        return TaskStatus::fail;
    }

    bool successful_prim_recovery = false;

    if(success == true){
        // Copy ALL variables from _init to _new. 
        const int nvar_u = U_init.GetDim(4);
        for (int v = 0; v < nvar_u; ++v) {
            U_new(v, k, j, i) = U_init(v, k, j, i);
        }

        const int nvar_p = P_init.GetDim(4);
        for (int v = 0; v < nvar_p; ++v) {
            P_new(v, k, j, i) = P_init(v, k, j, i);
        }

        U_new(m_u.UU_RAD, k, j, i) = U_init(m_u.UU_RAD, k, j, i) + dcov_rad[0];
        U_new(m_u.U1_RAD, k, j, i) = U_init(m_u.U1_RAD, k, j, i) + dcov_rad[1];
        U_new(m_u.U2_RAD, k, j, i) = U_init(m_u.U2_RAD, k, j, i) + dcov_rad[2];
        U_new(m_u.U3_RAD, k, j, i) = U_init(m_u.U3_RAD, k, j, i) + dcov_rad[3];
        // Here we are gonna overwrite only the active fluid conserved variables modified by the solver
        U_new(m_u.UU, k, j, i) = U_init(m_u.UU, k, j, i) - dcov_rad[0];
        U_new(m_u.U1, k, j, i) = U_init(m_u.U1, k, j, i) - dcov_rad[1];
        U_new(m_u.U2, k, j, i) = U_init(m_u.U2, k, j, i) - dcov_rad[2];
        U_new(m_u.U3, k, j, i) = U_init(m_u.U3, k, j, i) - dcov_rad[3];

        // Do GRMHD U2P with the new guess for the fluid primitives.
        // Real U_mhd_final[4] = {U_new(m_u.UU, k, j, i), U_new(m_u.U1, k, j, i), 
        //                             U_new(m_u.U2, k, j, i), U_new(m_u.U3, k, j, i)};
        //Real P_mhd_final[4];

        auto mhd_inverter_status = Inverter::u_to_p<Inverter::Type::kastaun>(
            G, U_new, m_u, gam, k, j, i, P_new, m_p, Loci::center, 
            25, 1e-12, false);
            
        if (mhd_inverter_status != static_cast<int>(Inverter::Status::success)) {
            successful_prim_recovery = false;
            // It failed, go back to what it was!
            for (int v = 0; v < nvar_p; ++v) {
                P_new(v, k, j, i) = P_init(v, k, j, i);
            }

        } else {
            successful_prim_recovery = true;

            // We don't need this, u_to_p will update P_new directly, I think.
            // P_new(m_p.UU, k, j, i) = P_mhd_final[0];
            // P_new(m_p.U1, k, j, i) = P_mhd_final[1];
            // P_new(m_p.U2, k, j, i) = P_mhd_final[2];
            // P_new(m_p.U3, k, j, i) = P_mhd_final[3];

            // Now since the P2U for MHD was successful, do it for radiation:
            Real U_rad_final[4] = {U_new(m_u.UU_RAD, k, j, i), U_new(m_u.U1_RAD, k, j, i), 
                                    U_new(m_u.U2_RAD, k, j, i), U_new(m_u.U3_RAD, k, j, i)};
            Real P_rad_final[4];

            auto rad_status = CalculateRadPrimitive_M1(G, U_rad_final, P_rad_final, k, j, i);

            if (rad_status == TaskStatus::fail) {
                successful_prim_recovery = false;
            } else {
                P_new(m_p.UU_RAD,  k, j, i) = P_rad_final[0];
                P_new(m_p.U1_RAD, k, j, i) = P_rad_final[1];
                P_new(m_p.U2_RAD, k, j, i) = P_rad_final[2];
                P_new(m_p.U3_RAD, k, j, i) = P_rad_final[3];
            }
        }

        if(!successful_prim_recovery) {
            // Define a safe numerical floor
            const Real floor_val = 10 * SMALL;

            // Reset fluid to floors
            // Density and energy to minimum, velocity to zero.
            P_new(m_p.RHO, k, j, i) = floor_val; 
            P_new(m_p.UU,  k, j, i) = floor_val;
            P_new(m_p.U1,  k, j, i) = 0.0;
            P_new(m_p.U2,  k, j, i) = 0.0;
            P_new(m_p.U3,  k, j, i) = 0.0;

            // Fetch B-fields to recalculate fluid conserved state
            
            B_P[0] = P_init(m_p.B1, k, j, i);
            B_P[1] = P_init(m_p.B2, k, j, i);
            B_P[2] = P_init(m_p.B3, k, j, i);
            Real U_mhd_floor[4];
            Real rho_ut_dummy;
            
            uvec[0] = P_new(m_p.U1, k, j, i);
            uvec[1] = P_new(m_p.U2, k, j, i);
            uvec[2] = P_new(m_p.U3, k, j, i);
            GRMHD::p_to_u_mhd(G, floor_val, floor_val, uvec, B_P, gam, k, j, i, 
                    rho_ut_dummy, U_mhd_floor, Loci::center);

            U_new(m_u.UU, k, j, i) = U_mhd_floor[0];
            U_new(m_u.U1, k, j, i) = U_mhd_floor[1];
            U_new(m_u.U2, k, j, i) = U_mhd_floor[2];
            U_new(m_u.U3, k, j, i) = U_mhd_floor[3];

            // Reset radiation floors
            // Radiation energy to minimum, flux to zero.
            P_new(m_p.UU_RAD, k, j, i) = floor_val;
            P_new(m_p.U1_RAD, k, j, i) = 0.0;
            P_new(m_p.U2_RAD, k, j, i) = 0.0;
            P_new(m_p.U3_RAD, k, j, i) = 0.0;

            // Because fluid velocity is exactly zero, fluid rest frame == lab frame.
            // All we have to do is multiply by gdet to make them Conserved Variables.
            U_new(m_u.UU_RAD, k, j, i) = floor_val * gdet;
            U_new(m_u.U1_RAD, k, j, i) = 0.0;
            U_new(m_u.U2_RAD, k, j, i) = 0.0;
            U_new(m_u.U3_RAD, k, j, i) = 0.0;
        }
    }else{
        return TaskStatus::fail;
    }
    return TaskStatus::complete;
}



 TaskStatus RadM1::Step(MeshData<Real> *md_full_init, MeshData<Real> *md_sub_init, MeshData<Real> *md_sub_final, const Real dt)
{
    for (int b = 0; b < md_sub_final->NumBlocks(); ++b) {
        auto pmb_data = md_sub_final->GetBlockData(b);
        auto pmb      = pmb_data->GetBlockPointer();
        auto& params  = pmb->packages.Get("RadM1")->AllParams();

        // // Fetch parameters
        const Real src_rootfind_eps   = params.Get<Real>("src_rootfind_eps");
        const Real src_rootfind_tol   = params.Get<Real>("src_rootfind_tol");
        const int src_rootfind_maxiter = params.Get<int>("src_rootfind_maxiter");
        const Real gam                 = pmb->packages.Get("GRMHD")->AllParams().Get<Real>("gamma");
        const auto& G                  = pmb->coords;

        PackIndexMap prims_map, cons_map;
        auto P_new = pmb_data->PackVariables({Metadata::GetUserFlag("Primitive")}, prims_map);
        auto U_new = pmb_data->PackVariables({Metadata::Conserved, Metadata::Cell}, cons_map);
        const VarMap m_p(prims_map, false);
        const VarMap m_u(cons_map, true);

        auto pmb_init_data = md_sub_init->GetBlockData(b);
        
        PackIndexMap dummy_map; // We can reuse or ignore this map for the initial state
        auto P_init = pmb_init_data->PackVariables({Metadata::GetUserFlag("Primitive")}, prims_map);
        auto U_init = pmb_init_data->PackVariables({Metadata::Conserved, Metadata::Cell}, cons_map);

        auto bounds = pmb->cellbounds;
        const IndexRange ib = bounds.GetBoundsI(IndexDomain::interior);
        const IndexRange jb = bounds.GetBoundsJ(IndexDomain::interior);
        const IndexRange kb = bounds.GetBoundsK(IndexDomain::interior);

        pmb->par_for("RadM1_Implicit_Solver4D", kb.s, kb.e, jb.s, jb.e, ib.s, ib.e,
            KOKKOS_LAMBDA (const int &k, const int &j, const int &i) {
                
                // Directly runs your 4D point-wise equation solver updating P_new and U_new in-place
                SolveImplicitRadiation4D(
                    G, U_init, P_init, P_new, U_new, m_p, m_u, k, j, i,
                    dt, gam, src_rootfind_eps, src_rootfind_tol, src_rootfind_maxiter
                );
            }
        );
    }

    return TaskStatus::complete;
}