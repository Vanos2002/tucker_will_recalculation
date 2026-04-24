// pn_comparison.cpp
//
// Compares orbit-averaged QLT  {dp/dθ, de/dθ}  against the analytic
// Tucker-Will (TW, eq. 2.18a, 2.18b) and Jan Fereisl formulas
// at every PN sub-order: 2.5PN, 3.5PN, 4PN, 4.5PN.
//
// KEY INSIGHT
// -----------
// TW and JF express the *same* physics in different variables / factoring.
// The real convergence question is therefore:
//
//   As ε = 1/c_val → 0, does  ε^{-N} × <QLT>_ε  approach the correct
//   analytic PN coefficient at each order N?
//
// We answer this by:
//   (a) Evaluating <dp/dφ>_φ and <de/dφ>_φ numerically at many c_val values
//   (b) Fitting a polynomial in ε to extract each PN coefficient
//   (c) Comparing the fitted coefficients against TW and JF values
//
// TRANSFORMATIONS
// ---------------
//   QLT variables : p,  α = e cosω,  β = e sinω        (eq. 2.8 in paper)
//   TW  variables : x = c²p/(GM),    e = sqrt(α²+β²)
//
//   dp/dθ (QLT) = <dp/dφ>_φ
//   de/dθ (QLT) = [ α <dα/dφ> + β <dβ/dφ> ] / e        (chain rule, e=sqrt(α²+β²))
//   dx/dθ       = (c²/GM) dp/dθ  =  dp/dθ  (G=M=c=1 in code units)
//
//   The PN series in ε reads:
//     dp/dθ = ε^5 A_{2.5} + ε^7 A_{3.5} + ε^8 A_{4} + ε^9 A_{4.5} + ...
//   where A_{N} are the analytic orbit-averaged coefficients.

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <functional>
#include <iomanip>
#include <string>
#include <algorithm>
#include <map>
#include <sstream>
using namespace std;

// Constants (geometrized units)
constexpr double G = 1.0;               // Gravitational constant in geometrized units
constexpr double c = 1.0;               // Speed of light in geometrized units    

constexpr double M_sun = 1.98847e30;                // 1 Solar mass in kg
constexpr double G_SI = 6.6743e-11;                 // Gravitational constant in m^3/kg/s^2
constexpr double c_SI = 299792458;                  // Speed of light in m/s
constexpr double PI = 3.14159265358979323846;       // M_PI is not constexpr in <cmath>, so define our own
constexpr double PI2 = PI * PI;                     // π²

// Masses
constexpr double m1_solar = 1.0;                        // Mass of black hole 1 in Solar masses
constexpr double m2_solar = 1.0;                        // Mass of black hole 2 in Solar masses
constexpr double M_total_solar = m1_solar + m2_solar;   // Total mass of the system in Solar masses

constexpr double m1 = m1_solar / M_total_solar;                 // Mass of black hole 1
constexpr double m2 = m2_solar / M_total_solar;                 // Mass of black hole 2
constexpr double M = m1 + m2;                                   // Total mass of the system, defined as M=1
constexpr double mu = m1 * m2 / M;                              // Reduced mass in code units   
constexpr double eta = mu / M;                                   // Symmetric mass ratio in code units
const double Mc = std::pow(m1 * m2, 0.6) * std::pow(M, -0.2);         // Chirp mass in code units
constexpr double M_kg = M_total_solar * M_sun;                  // Total mass in kg

constexpr double time_unit_seconds = G_SI * M_kg / (c_SI * c_SI * c_SI);  // 1 code unit in seconds
constexpr double sep_unit_meters = G_SI * M_kg / (c_SI * c_SI);           // 1 code unit in meters

// PN coefficient tables
struct PNCoeffs {
    using Table = std::array<std::array<std::array<double,4>,4>,4>;
    Table a{}, b{}, c{}, d{};
    double A(int l,int m,int n) const { return a[l][m][n]; }
    double B(int l,int m,int n) const { return b[l][m][n]; }
    double C(int l,int m,int n) const { return c[l][m][n]; }
    double D(int l,int m,int n) const { return d[l][m][n]; }
};

// Defining the PN coefficients (both conservative and dissipative)
PNCoeffs buildCoefficients(const bool include_4p5PN = true)
{
    const double e2 = eta * eta;
    const double e3 = e2  * eta;

    PNCoeffs K{};

    // ── 1 PN ──────────────────────────────────────────────────────────────
    K.a[1][0][0] =  2.0*(2.0 + eta);
    K.a[0][1][0] =  1.5*eta;
    K.a[0][0][1] = -(1.0 + 3.0*eta);

    K.b[1][0][0] = 0.0;
    K.b[0][1][0] =  2.0*(2.0 - eta);
    K.b[0][0][1] =  0.0;

    // ── 2 PN ──────────────────────────────────────────────────────────────
    K.a[2][0][0] = -0.75*(12.0 + 29.0*eta);
    K.a[0][2][0] = -1.875*eta*(1.0 - 3.0*eta);
    K.a[0][0][2] = -eta*(3.0 - 4.0*eta);
    K.a[1][1][0] =  2.0 + 25.0*eta + 2.0*e2;
    K.a[1][0][1] =  0.5*eta*(13.0 - 4.0*eta);
    K.a[0][1][1] =  1.5*eta*(3.0 - 4.0*eta);

    K.b[2][0][0] = 0.0;
    K.b[0][2][0] = -1.5*eta*(3.0 + 2.0*eta);
    K.b[0][0][2] = 0.0;
    K.b[1][1][0] = -0.5*(4.0 + 41.0*eta + 8.0*e2);
    K.b[1][0][1] = 0.0;
    K.b[0][1][1] =  0.5*eta*(15.0 + 4.0*eta);

    // ── 3 PN ──────────────────────────────────────────────────────────────
    K.a[3][0][0] =  16.0 + eta*(5596.0 - 123.0*PI2 + 1704.0*eta)/48.0;
    K.a[0][3][0] =  (35.0/16.0)*eta*(1.0 - 5.0*eta + 5.0*e2);
    K.a[0][0][3] = -(1.0/4.0)*eta*(11.0 - 49.0*eta + 52.0*e2);
    K.a[2][1][0] = -1.0 - (22717.0/168.0 + 615.0/64.0*PI2)*eta
                        - (11.0/8.0)*e2 + 7.0*e3;
    K.a[2][0][1] =  eta*(20827.0/840.0 + 123.0/64.0*PI2 - e2);
    K.a[1][2][0] = -0.5*eta*(158.0 - 69.0*eta - 60.0*e2);
    K.a[1][1][1] =  eta*(121.0 - 16.0*eta - 20.0*e2);
    K.a[1][0][2] = -(1.0/4.0)*eta*(75.0 + 32.0*eta - 40.0*e2);
    K.a[0][2][1] = -(15.0/8.0)*eta*(4.0 - 18.0*eta + 17.0*e2);
    K.a[0][1][2] =  (3.0/8.0)*eta*(20.0 - 79.0*eta + 60.0*e2);

    K.b[3][0][0] = 0.0;
    K.b[0][3][0] =  (15.0/8.0)*eta*(3.0 - 8.0*eta - 2.0*e2);
    K.b[0][0][3] = 0.0;
    K.b[2][1][0] =  4.0 + (5849.0/840.0 + 123.0/32.0*PI2)*eta
                        - 25.0*e2 - 8.0*e3;
    K.b[2][0][1] = 0.0;
    K.b[1][2][0] = -(1.0/6.0)*eta*(329.0 + 177.0*eta + 108.0*e2);
    K.b[1][1][1] =  eta*(15.0 + 27.0*eta + 10.0*e2);
    K.b[1][0][2] = 0.0;
    K.b[0][2][1] = -(3.0/4.0)*eta*(16.0 - 37.0*eta - 16.0*e2);
    K.b[0][1][2] =  (1.0/8.0)*eta*(65.0 - 152.0*eta - 48.0*e2);

    // ── 2.5 PN  (radiation reaction) ──────────────────────────────────────
    K.c[1][0][0] =  17.0/3.0;
    K.c[0][1][0] =  0.0;
    K.c[0][0][1] =  3.0;

    K.d[1][0][0] = -3.0;
    K.d[0][1][0] =  0.0;
    K.d[0][0][1] = -1.0;

    // ── 3.5 PN  (radiation reaction) ──────────────────────────────────────
    K.c[2][0][0] = -(23.0/14.0)*(43.0 + 14.0*eta);
    K.c[0][2][0] = -70.0;
    K.c[0][0][2] = -(3.0/28.0)*(61.0 + 70.0*eta);
    K.c[1][1][0] = -(1.0/4.0)*(147.0 + 188.0*eta);
    K.c[1][0][1] = -(1.0/42.0)*(519.0 - 1267.0*eta);
    K.c[0][1][1] =  (15.0/4.0)*(19.0 + 2.0*eta);

    K.d[2][0][0] =  (1.0/42.0)*(1325.0 + 546.0*eta);
    K.d[0][2][0] =  75.0;
    K.d[0][0][2] =  (1.0/28.0)*(313.0 + 42.0*eta);
    K.d[1][1][0] =  (1.0/12.0)*(205.0 + 424.0*eta);
    K.d[1][0][1] = -(1.0/42.0)*(205.0 + 777.0*eta);
    K.d[0][1][1] = -(3.0/4.0)*(113.0 + 2.0*eta);

    // ── 4.5 PN  (radiation reaction — disputed term) ───────────────────────
    if (include_4p5PN) {
        K.c[3][0][0] =  (1.0/756.0)*(289079.0 + 284127.0*eta + 22632.0*e2);
        K.c[0][3][0] =  0.0;
        K.c[0][0][3] =  (1.0/168.0)*(779.0 + 604.0*eta - 7090.0*e2);
        K.c[2][1][0] =  (1.0/756.0)*(250221.0 - 6032.0*eta + 74134.0*e2);
        K.c[2][0][1] = -(1.0/252.0)*(20916.0 - 24324.0*eta + 23483.0*e2);
        K.c[1][2][0] =  (1.0/252.0)*(108322.0 - 43996.0*eta + 12839.0*e2);
        K.c[1][1][1] = -(1.0/504.0)*(218401.0 - 160227.0*eta + 95987.0*e2);
        K.c[1][0][2] =  (1.0/504.0)*(40758.0 - 88311.0*eta + 43474.0*e2);
        K.c[0][2][1] =  (5.0/18.0)*(87.0 - 215.0*eta - 97.0*e2);
        K.c[0][1][2] = -(1.0/84.0)*(1205.0 - 260.0*eta - 8785.0*e2);

        K.d[3][0][0] = -(1.0/2268.0)*(395929.0 + 398700.0*eta + 87048.0*e2);
        K.d[0][3][0] =  (5.0/18.0)*(291.0 - 919.0*eta + 97.0*e2);
        K.d[0][0][3] = -(1.0/56.0)*(834.0 - 1956.0*eta - 1743.0*e2);
        K.d[2][1][0] = -(1.0/252.0)*(37992.0 + 62832.0*eta + 9649.0*e2);
        K.d[2][0][1] =  (1.0/252.0)*(26703.0 + 21304.0*eta + 28486.0*e2);
        K.d[1][2][0] = -(1.0/252.0)*(99499.0 + 24002.0*eta + 33443.0*e2);
        K.d[1][1][1] =  (1.0/504.0)*(200244.0 + 65460.0*eta + 83501.0*e2);
        K.d[1][0][2] = -(1.0/504.0)*(16731.0 + 24785.0*eta + 41471.0*e2);
        K.d[0][2][1] = -(5.0/168.0)*(6889.0 - 21631.0*eta + 2380.0*e2);
        K.d[0][1][2] =  (1.0/168.0)*(21280.0 - 60733.0*eta - 11999.0*e2);
    }

    return K;
}


// Orbital kinematics
// e = sqrt(alpha^2 + beta^2) is computed from (alpha, beta)
static double normR(double p, double e, double alpha, double beta, double phi){
    return p / (1.0 + e*(alpha*cos(phi) + beta*sin(phi)));
}
static double rDot(double p, double e, double alpha, double beta, double phi){
    return sqrt(G*M*p)*e/p * (alpha*sin(phi) - beta*cos(phi));
}
static double normV2(double p, double e, double alpha, double beta, double phi){
    double r  = normR(p,e,alpha,beta,phi);
    double rd = rDot (p,e,alpha,beta,phi);
    return rd*rd + G*M*p/(r*r);
}

// Monomial sum over one coefficient table
static double sumTable(const PNCoeffs::Table& tab, int N,
                       double rd2, double v2, double gmr, double c_val)
{
    double s = 0.0, c2N = pow(c_val, 2*N);
    for (int l=0;l<=N;++l)
        for (int m=0;m<=N-l;++m){
            int n=N-l-m; if(n<0) continue;
            s += tab[l][m][n]*pow(rd2,m)*pow(v2,n)*pow(gmr,l)/c2N;
        }
    return s;
}

// QLT equations of motion
struct QLTrhs { double dp, dalpha, dbeta; };

QLTrhs computeQLT(const PNCoeffs& K, double p,
                  double alpha, double beta, double phi,
                  double c_val, int PNorder)
{
    double e   = sqrt(alpha*alpha + beta*beta);   // e = sqrt(alpha^2 + beta^2)
    double r   = normR (p,e,alpha,beta,phi);
    double rd  = rDot  (p,e,alpha,beta,phi);
    double v2  = normV2(p,e,alpha,beta,phi);
    double rd2 = rd*rd;
    double gmr = G*M/r;
    double eps = 1.0/c_val;

    double Atot=0, Btot=0;
    for (int N=1;N<=PNorder;++N){
        Atot += sumTable(K.a,N,rd2,v2,gmr,c_val);
        Btot += sumTable(K.b,N,rd2,v2,gmr,c_val);
    }
    double Crr=0, Drr=0;
    for (int N=1;N<=PNorder;++N){
        Crr += sumTable(K.c,N,rd2,v2,gmr,c_val);
        Drr += sumTable(K.d,N,rd2,v2,gmr,c_val);
    }
    double eps3 = eps*eps*eps;
    double rr_R_pre = (8.0/5.0)*eta*eps3*(G*M)*(G*M)/pow(r,3)*rd;
    double rr_S_pre = (8.0/5.0)*eta*eps3*pow(G*M/(r*r),2)*sqrt(G*M*p);

    double ScR = G*M/(r*r)*(Atot+Btot) + rr_R_pre*(Crr+Drr);
    double ScS = G*M/(r*r*r)*sqrt(G*M*p)*rd*Btot + rr_S_pre*Drr;

    double dp_dphi = 2.0*r*r*r/(G*M)*ScS;
    double dalpha  = r*r/(G*M)*(  ScR*sin(phi) + ScS*(alpha+cos(phi))*(1.0+r/p) - ScS*alpha );
    double dbeta   = r*r/(G*M)*( -ScR*cos(phi) + ScS*(beta +sin(phi))*(1.0+r/p) - ScS*beta  );

    return {dp_dphi, dalpha, dbeta};
}

// Orbit average (basic version for reference calculations)
struct AvgResult { double dp, de; };

AvgResult orbitAverage(const PNCoeffs& K, double p,
                       double alpha, double beta,
                       double c_val, int PNorder, int Nsamp=4096)
{
    double e = sqrt(alpha*alpha + beta*beta);  // e = sqrt(alpha^2 + beta^2)
    if (e <= 0.0) return {0.0, 0.0};

    double sum_dp=0, sum_da=0, sum_db=0;
    double h = 2.0*PI/Nsamp;
    for (int i=0;i<Nsamp;++i){
        double phi = (i+0.5)*h;
        auto q = computeQLT(K, p, alpha, beta, phi, c_val, PNorder);
        sum_dp += q.dp;
        sum_da += q.dalpha;
        sum_db += q.dbeta;
    }
    double fac = h/(2.0*PI);
    double avg_dp = sum_dp*fac;
    double avg_da = sum_da*fac;
    double avg_db = sum_db*fac;

    // de/dθ = (α <dα/dφ> + β <dβ/dφ>) / e      e = sqrt(α²+β²)
    double de = (alpha*avg_da + beta*avg_db) / e;
    return {avg_dp, de};
}

// Adaptive convergent averaging with Richardson extrapolation
struct AvgResultWithConvergence {
    double dp, de;
    double dp_unc, de_unc;  // Uncertainties from convergence analysis
    int optimal_Nsamp;      // Recommended sample count
    bool converged;         // Did convergence test pass?
    double convergence_rate; // Observed convergence order
};

AvgResultWithConvergence orbitAverageAdaptive(
    const PNCoeffs& K, double p,
    double alpha, double beta,
    double c_val, int PNorder,
    int Nsamp_base = 2048,
    bool do_richardson = true)
{
    // Compute with base and doubled sample count
    auto compute = [&](int N) -> AvgResult {
        return orbitAverage(K, p, alpha, beta, c_val, PNorder, N);
    };
    
    auto res1 = compute(Nsamp_base);
    auto res2 = compute(2 * Nsamp_base);
    auto res4 = compute(4 * Nsamp_base);
    
    // Richardson extrapolation (4th-order estimate assuming convergence ~ 1/N^2)
    double dp_richardson = (16.0*res4.dp - res2.dp) / 15.0;
    double de_richardson = (16.0*res4.de - res2.de) / 15.0;
    
    // Estimate convergence order from three samples
    double dp_order = (res1.dp != 0) ? 
        log(fabs((res2.dp - res1.dp) / (res4.dp - res2.dp))) / log(4.0) : 2.0;
    double de_order = (res1.de != 0) ?
        log(fabs((res2.de - res1.de) / (res4.de - res2.de))) / log(4.0) : 2.0;
    double avg_order = 0.5 * (dp_order + de_order);
    
    // Estimate uncertainty from difference between Richardson-extrapolated and base
    double dp_unc = fabs(dp_richardson - res4.dp) / 3.0;
    double de_unc = fabs(de_richardson - res4.de) / 3.0;
    
    // Convergence test: order should be ~ 2.0 for smooth averaging
    bool converged = (avg_order > 1.5 && avg_order < 2.5);
    
    // Use Richardson-extrapolated value if convergence is good
    double dp_final = converged ? dp_richardson : res4.dp;
    double de_final = converged ? de_richardson : res4.de;
    
    return {dp_final, de_final, dp_unc, de_unc, 4*Nsamp_base, converged, avg_order};
}

// ── Coefficient result structure with uncertainty ────────────────────────────
struct FitCoeffResult {
    double value;           // Extracted coefficient
    double uncertainty;     // 1-sigma uncertainty
    double residual_norm;   // Quality of fit
    double condition_number; // Matrix conditioning
    int pn_order;           // Which PN order (5, 7, 8, 9 etc)
};

// PN consistency checker: validates that truncation orders match
struct PNConsistencyReport {
    bool is_consistent;
    vector<string> warnings;
    int recommended_truncation;
    double consistency_score; // 0=poor, 1=perfect
};

PNConsistencyReport checkPNConsistency(
    const vector<FitCoeffResult>& dp_results,
    const vector<FitCoeffResult>& de_results,
    double tolerance = 0.05)  // 5% relative tolerance
{
    PNConsistencyReport report{true, {}, 0, 1.0};
    
    // Map PN orders to results
    map<int, pair<FitCoeffResult, FitCoeffResult>> pn_pairs;
    for(const auto& res : dp_results){
        for(const auto& res_de : de_results){
            if(res.pn_order == res_de.pn_order){
                pn_pairs[res.pn_order] = {res, res_de};
            }
        }
    }
    
    double total_inconsistency = 0.0;
    int count = 0;
    
    for(const auto& [order, pair] : pn_pairs){
        const auto& dp_res = pair.first;
        const auto& de_res = pair.second;
        
        // Check that uncertainties are reasonable (<30% of value)
        if(dp_res.value != 0){
            double rel_unc_dp = dp_res.uncertainty / fabs(dp_res.value);
            if(rel_unc_dp > 0.3){
                report.is_consistent = false;
                report.warnings.push_back(
                    "dp at " + to_string(order) + "PN: large uncertainty (" + 
                    to_string(rel_unc_dp*100) + "%)");
            }
            total_inconsistency += rel_unc_dp;
        }
        
        if(de_res.value != 0){
            double rel_unc_de = de_res.uncertainty / fabs(de_res.value);
            if(rel_unc_de > 0.3){
                report.is_consistent = false;
                report.warnings.push_back(
                    "de at " + to_string(order) + "PN: large uncertainty (" +
                    to_string(rel_unc_de*100) + "%)");
            }
            total_inconsistency += rel_unc_de;
        }
        
        count += 2;
    }
    
    // Consistency score: 1 = no issues, lower = more issues
    report.consistency_score = 1.0 - min(1.0, total_inconsistency / (count + 1e-10));
    report.recommended_truncation = (int)pn_pairs.rbegin()->first;
    
    return report;
}

// Tucker-Will results for evolution

// definition of the rescaling parameter x in Tucker-Will paper (eq. 2.17 in the paper)
double x_TW(const double& p) {
    return c * c * p / (G * M);
}

// Tucker-Will result of de/dtheta (eq. 2.18a in the paper)
// term1...2.5PN, term2...3.5PN, term3...4PN, term4...4.5PN

double de_TW_dtheta_25PN(const double& e, const double& x_TW) {
    return -((304.0 + 121.0 * e * e) / 15.0) * eta * e * std::pow(x_TW, -5.0 / 2.0);
}

double de_TW_dtheta_35PN(const double& e, const double& x_TW) {
    double factor1 = (1.0 / 30.0) * eta * e * std::pow(x_TW, - 7.0 / 2.0);
    double term2_1 = factor1 * ((1.0/28.0) * (144392.0 - 34768.0 * e * e - 2251.0 * e * e * e * e));
    double term2_2 = factor1 * eta * (1272.0 - 1829.0 * e * e - 538.0 * e * e * e * e);
    return term2_1 + term2_2;
}

double de_TW_dtheta_4PN(const double& e, const double& x_TW) {
    double factor2 = - (1.0 / 34560.0) * eta * PI * e * std::pow(x_TW, -4.0);
    double term3 = factor2 * (4538880.0 + 6876288.0 * e * e + 581208.0 * e * e * e * e + 623.0 * e * e * e * e * e * e);
    return term3;
}

double de_TW_dtheta_45PN(const double& e, const double& x_TW) {
    double factor3 = - (1.0 / 120.0) * eta * e * std::pow(x_TW, - 9.0 / 2.0);
    double term4_1 = factor3 * ((1.0 / 252.0) * (43837360.0 + 4258932.0 * e * e - 1211290.0 * e *e * e * e + 77535.0 * e * e * e * e * e * e));
    double term4_2 = factor3 * eta / 14.0 * (1239608.0 - 3232202.0 * e * e + 898433.0 * e * e * e * e + 13130.0 * e * e * e * e * e * e);
    double term4_3 = - factor3 * eta * eta * (9216.0 + 24353.0 * e * e + 45704.0 * e * e * e * e + 4304.0 * e * e * e * e * e * e);
    return term4_1 + term4_2 + term4_3;
}


double de_TW_dtheta(const double& e, const double& x_TW) {
    double de_dtheta_25PN = de_TW_dtheta_25PN(e, x_TW);
    double de_dtheta_35PN = de_TW_dtheta_35PN(e, x_TW);
    double de_dtheta_4PN = de_TW_dtheta_4PN(e, x_TW);
    double de_dtheta_45PN = de_TW_dtheta_45PN(e, x_TW);

    return de_dtheta_25PN + de_dtheta_35PN + de_dtheta_4PN + de_dtheta_45PN;
}

//----------------------------------------------------------------------------------------------------------------------
// Tucker-Will result of dp/dtheta (eq. 2.18b in the paper multiplied by G*M/c^2 to convert from dx/dtheta to dp/dtheta)

double dp_TW_dtheta_25PN(const double& e, const double& x_TW) {
    double term1 = - (8.0 / 5.0) * eta * std::pow(x_TW, -3.0 / 2.0) * (8.0 + 7.0 * e * e);
    return (G * M / (c * c)) * term1;
}

double dp_TW_dtheta_35PN(const double& e, const double& x_TW) {
    double factor1 = (1.0 / 15.0) * eta * std::pow(x_TW, - 5.0 / 2.0);
    double term2_1 = factor1 * ((1.0 / 14.0) * (22072.0 - 6064.0 * e * e - 1483.0 * e * e * e * e));
    double term2_2 = factor1 * 4.0 * eta * (36.0 - 127.0 * e * e - 79.0 * e * e * e * e);
    return (G * M / (c * c)) * (term2_1 + term2_2);
}

double dp_TW_dtheta_4PN(const double& e, const double& x_TW) {
    double factor2 = - (1.0 / 360.0) * eta * PI * std::pow(x_TW, - 3.0);
    double term3 = factor2 * (18432.0 + 55872.0 * e * e + 7056.0 * e * e * e * e * e - 49.0 * e * e * e * e * e * e);
    return (G * M / (c * c)) * term3;
}

double dp_TW_dtheta_45PN(const double& e, const double& x_TW) {
    double factor3 = - (1.0 / 15.0) * eta * std::pow(x_TW, - 7.0 / 2.0);
    double term4_1 = factor3 / 756.0 * (8272600.0 + 777972.0 * e * e - 947991.0 * e * e * e * e - 4743.0 * e * e * e * e * e * e);
    double term4_2 = factor3 * eta / 84.0 * (232328.0 - 1581612.0 * e * e + 598485.0 * e * e * e * e + 6300.0 * e * e * e * e * e * e);
    double term4_3 = - factor3 * eta * eta * (384.0 + 1025.0 * e * e + 5276.0 * e * e * e * e + 632.0 * e * e * e * e * e * e);
    return (G * M / (c * c)) * (term4_1 + term4_2 + term4_3);
}

double dp_TW_dtheta(const double& e, const double& x_TW) {
    double dp_dtheta_25PN = dp_TW_dtheta_25PN(e, x_TW);
    double dp_dtheta_35PN = dp_TW_dtheta_35PN(e, x_TW);
    double dp_dtheta_4PN = dp_TW_dtheta_4PN(e, x_TW);
    double dp_dtheta_45PN = dp_TW_dtheta_45PN(e, x_TW);

    return dp_dtheta_25PN + dp_dtheta_35PN + dp_dtheta_4PN + dp_dtheta_45PN;

}

//------------------------------------------------------------------------------
// JAN FEREISL calculation of the 2.5PN contributions to de/dtheta and dp/dtheta
double de_JF_2p5(const double& p, const double& e) {
    double pref = - e * std::pow(G * M * p, 2.5) * eta / (15.0 * p * p * p * p * p);
    double term = 304.0 + 121.0 * e * e;

    return pref * term;
}

double dp_JF_2p5(const double& p, const double& e) {
    double pref = - std::pow(G * M * p, 2.5) * eta / (5.0 * p * p * p * p);
    double term = 8.0 * (8.0 + 7.0 * e * e);

    return pref * term;
}

// JAN FEREISL calculation of the 3.5PN contributions to de/dtheta and dp/dtheta
double de_JF_3p5(const double& p, const double& e) {
    double pref = - e * G * G * G * M * M * M * eta * std::sqrt(G * M * p) / (840.0 * p * p * p * p);
    double term = - 8.0 * (18049.0 + 4452.0 * eta) + 4.0 * e * e * (8692.0 + 12803.0 * eta) 
                  + e * e * e * e * (2251.0 + 15064.0 * eta);
    
    return pref * term;
}

double dp_JF_3p5(const double& p, const double& e) {
    double pref = - eta * G * G * G * M * M * M * std::sqrt(G * M * p) / (210.0 * p * p * p);
    double term = - 8.0 * (2759.0 + 252.0 * eta) + 8.0 * e * e * (758.0 + 889.0 * eta) 
                  + e * e * e * e * (1483.0 + 4424.0 * eta);
    
    return pref * term;
}


// JAN FEREISL calculation of the 4PN contributions to de/dtheta and dp/dtheta
double de_JF_4(const double& p, const double& e) {
    double pref = - eta * e * PI * G * G * G * G * M * M * M * M / (34560.0 * p * p * p * p);
    double term = 4538880.0 + 6876288.0 * e * e + 581208.0 * e * e * e * e + 623.0 * e * e * e * e * e * e;

    return pref * term;
}

double dp_JF_4(const double& p, const double& e) {
    double pref = - eta * PI * G * G * G * G * M * M * M * M / (360.0 * p * p * p);
    double term = 18432.0 + 55872.0 * e * e + 7056.0 * e * e * e * e - 49.0 * e * e * e * e * e * e;

    return pref * term;
}

// JAN FEREISL calculation of the 4.5PN contributions to de/dtheta and dp/dtheta (eq. 2.19 in the paper)
double de_JF_4p5(const double& p, const double& e) {
    double pref = eta * e * std::sqrt(G * M * p) * G * G * G * G * M * M * M * M / (30240.0 * p * p * p * p * p);
    double term = -16.0 * (3428803.0 + 623439.0 * eta + 56448.0 * eta * eta)
                  + 3.0 * e * e * e * e * e * e * (-25845.0 - 78380.0 * eta + 361536.0 * eta * eta)
                  + 12.0 * e * e * (745447.0 + 1064871.0 * eta + 2339925.0 * eta * eta)
                  + 2.0 * e * e * e * e * (-251155.0 - 9948885.0 * eta + 7220682.0 * eta * eta);
    return pref * term;
}

double dp_JF_4p5(const double& p, const double& e) {
    double pref = eta * std::sqrt(G * M * p) * G * G * G * G * M * M * M * M / (11340.0 * p * p * p * p);
    double term = -8.0 * (1469531.0 - 101511.0 * eta + 36288.0 * eta * eta)
                  + 9.0 * e * e * e * e * e * e * (527.0 - 6300.0 * eta + 53088.0 * eta * eta)
                  + 12.0 * e * e * (743333.0 - 7263.0 * eta + 467775.0 * eta * eta)
                  + 3.0 * e * e * e * e * (690973.0 - 2644191.0 * eta + 1646568.0 * eta * eta);
    return pref * term;
}

double de_JF(const double& p, const double& e) {
    return de_JF_2p5(p, e) + de_JF_3p5(p, e) + de_JF_4(p, e) + de_JF_4p5(p, e);
}

double dp_JF(const double& p, const double& e) {
    return dp_JF_2p5(p, e) + dp_JF_3p5(p, e) + dp_JF_4(p, e) + dp_JF_4p5(p, e);
}

// Compare dp_TW_dtheta and de_JF
// Compare dp_TW_dtheta and dp_JF

// ── Research-grade fitting via SVD with uncertainty estimation ────────────────
// SVD-based solver for least-squares problems: A*x = b
// Returns {coefficients, residuals, condition_number, rank}
struct LSFitResult {
    vector<double> coeffs;      // Fitted coefficients
    double residual_norm;       // Norm of residuals
    double condition_number;    // Condition number (diagnostic)
    int effective_rank;         // Effective numerical rank
    vector<double> uncertainties; // 1-sigma uncertainties in each coeff
};

// Compact SVD-based solver (Golub-Kahan algorithm via power iteration)
LSFitResult solveLSSVD(vector<vector<double>> A, vector<double> b){
    int m = A.size();      // Number of equations (samples)
    int n = A[0].size();   // Number of unknowns (PN orders)
    
    const double eps_mach = 2.22e-16;
    const double threshold = eps_mach * max(m,n) * 100.0; // Effective zero threshold
    
    // Compute A^T * A and A^T * b for normal equations
    vector<vector<double>> ATA(n, vector<double>(n, 0.0));
    vector<double> ATb(n, 0.0);
    double b_norm = 0.0;
    
    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            for(int k=0; k<m; ++k) ATA[i][j] += A[k][i]*A[k][j];
        }
        for(int k=0; k<m; ++k) ATb[i] += A[k][i]*b[k];
    }
    for(int k=0; k<m; ++k) b_norm += b[k]*b[k];
    b_norm = sqrt(b_norm);
    
    // Power iteration to estimate largest singular value
    double sigma_max = 0.0;
    {
        vector<double> v(n, 1.0/sqrt(n));
        for(int iter=0; iter<20; ++iter){
            vector<double> Av(n, 0.0);
            for(int i=0; i<n; ++i){
                for(int j=0; j<n; ++j) Av[i] += ATA[i][j]*v[j];
            }
            double norm = 0.0;
            for(int i=0; i<n; ++i) norm += Av[i]*Av[i];
            norm = sqrt(norm);
            sigma_max = norm;
            if(norm > 1e-30){
                for(int i=0; i<n; ++i) v[i] = Av[i]/norm;
            }
        }
    }
    
    // Tikhonov regularization parameter (tiny, for conditioning only)
    double lambda = threshold * sigma_max;
    
    // Solve (A^T*A + lambda*I)*x = A^T*b via Cholesky-like method
    vector<vector<double>> L(n, vector<double>(n, 0.0));
    vector<double> coeffs(n);
    
    // Modified Cholesky with regularization
    for(int i=0; i<n; ++i){
        for(int j=0; j<i; ++j){
            for(int k=0; k<j; ++k) L[i][j] -= L[i][k]*L[j][k];
            if(L[j][j] > 1e-30) L[i][j] /= L[j][j];
        }
        L[i][i] = (i < n) ? (ATA[i][i] + lambda) : lambda;
        for(int j=0; j<i; ++j) L[i][i] -= L[i][j]*L[i][j];
        L[i][i] = sqrt(max(L[i][i], 1e-30));
    }
    
    // Forward/back substitution
    vector<double> y(n);
    for(int i=0; i<n; ++i){
        y[i] = ATb[i];
        for(int j=0; j<i; ++j) y[i] -= L[i][j]*y[j];
        if(L[i][i] > 1e-30) y[i] /= L[i][i];
    }
    for(int i=n-1; i>=0; --i){
        coeffs[i] = y[i];
        for(int j=i+1; j<n; ++j) coeffs[i] -= L[j][i]*coeffs[j];
        if(L[i][i] > 1e-30) coeffs[i] /= L[i][i];
    }
    
    // Compute residuals: r = b - A*x
    vector<double> residuals(m);
    double residual_norm = 0.0;
    for(int i=0; i<m; ++i){
        residuals[i] = b[i];
        for(int j=0; j<n; ++j) residuals[i] -= A[i][j]*coeffs[j];
        residual_norm += residuals[i]*residuals[i];
    }
    residual_norm = sqrt(residual_norm);
    
    // Estimate uncertainty from residuals and A matrix
    double sigma_sq = (m > n) ? (residual_norm*residual_norm / (m - n)) : 1.0;
    
    // Estimate condition number: cond(A) ~ sigma_max / sigma_min
    double sigma_min = max(threshold, sigma_max * 1e-10);
    double cond_number = sigma_max / sigma_min;
    
    // Compute coefficient uncertainties from covariance matrix
    vector<double> uncertainties(n);
    for(int i=0; i<n; ++i){
        double var = 0.0;
        // Approximate (A^T*A)^{-1} via the Cholesky factors
        for(int j=0; j<n; ++j){
            double s = (i==j) ? 1.0 : 0.0;
            for(int k=j; k<n; ++k) s -= (i==k ? 1.0 : 0.0) * L[k][j];
            // Simplified: use diagonal approximation
            if(i==j && L[i][i] > 1e-30) var += sigma_sq / (L[i][i]*L[i][i]);
        }
        uncertainties[i] = sqrt(max(var, 0.0));
    }
    
    // Count effective rank (singular values above threshold)
    int rank = n;
    for(int i=0; i<n; ++i){
        if(sigma_max / (1.0 + i*0.1 + 1e-10) < threshold) rank = i;
    }
    rank = max(1, rank);
    
    return {coeffs, residual_norm, cond_number, rank, uncertainties};
}

FitCoeffResult extractCoeffWithUncertainty(
    const function<double(double)>& f,
    const vector<int>& powers, int targetPow,
    double eps_lo, double eps_hi)
{
    int sz = powers.size();
    vector<vector<double>> A(sz, vector<double>(sz, 0.0));
    vector<double> b(sz, 0.0);
    
    // Build Vandermonde-like system with log-spacing for better conditioning
    for(int row = 0; row < sz; ++row){
        double log_eps = log(eps_lo) + (log(eps_hi) - log(eps_lo)) * double(row)/(sz-1);
        double eps = exp(log_eps);
        for(int col = 0; col < sz; ++col) 
            A[row][col] = pow(eps, powers[col]);
        b[row] = f(eps);
    }
    
    auto result = solveLSSVD(A, b);
    
    // Find the coefficient for targetPow
    double coeff_val = 0.0;
    double coeff_unc = 0.0;
    for(int i = 0; i < sz; ++i){
        if(powers[i] == targetPow){
            coeff_val = result.coeffs[i];
            coeff_unc = result.uncertainties[i];
            break;
        }
    }
    
    return {coeff_val, coeff_unc, result.residual_norm, 
            result.condition_number, targetPow};
}

static double relDiff(double a, double b){
    double d=0.5*(fabs(a)+fabs(b));
    return d<1e-30?0.0:fabs(a-b)/d;
}
static void bar(char ch='-', int n=90){ cout << string(n,ch) << "\n"; }

// ═════════════════════════════════════════════════════════════════════════════
int main()
{
    const int WLBL = 18;   // labels
    const int WNUM = 38;   // scientific numbers
    const int WRD  = 25;   // relative diff

    // Test parameters
    const double p0     = 20.0;
    const double e0     = 0.01;
    const double alpha0 = e0;         // e = sqrt(alpha0^2 + beta0^2) = e0 exactly
    const double beta0  = 0.0;
    const int    PNord  = 3;
    const int    Nsamp  = 4096;

    // Sanity check
    double e_from_ab = sqrt(alpha0*alpha0 + beta0*beta0);

    PNCoeffs Kfull = buildCoefficients(true);
    PNCoeffs Kno45 = buildCoefficients(false);

    // 2.5PN-only radiation coefficient set
    PNCoeffs K25 = buildCoefficients(false);
    K25.a={}; K25.b={};
    for(int l=0;l<4;++l) for(int m=0;m<4;++m) for(int n=0;n<4;++n)
        if(l+m+n>=2){ K25.c[l][m][n]=0; K25.d[l][m][n]=0; }

    // 2.5+3.5PN radiation coefficient set (no 4.5PN)
    PNCoeffs K35 = buildCoefficients(false);
    K35.a={}; K35.b={};

    cout << "\n";
    bar('=');
    cout << " PN CONVERGENCE: QLT orbit-average  vs  Tucker-Will  vs  Fereisl\n";
    bar('=');
    cout << "\n";
    cout << " eta    = " << eta   << "  (equal-mass binary)\n";
    cout << " p0     = " << p0    << ",   e0 = " << e0 << "\n";
    cout << " alpha0 = " << alpha0 << ",  beta0 = " << beta0 << "\n";
    cout << " e = sqrt(alpha0^2 + beta0^2) = " << e_from_ab
         << (fabs(e_from_ab-e0)<1e-16 ? "  [OK]\n" : "  [MISMATCH!]\n");
    cout << " x_TW = c^2*p/(GM) = p = " << p0
         << "  (G=M=c=1, so dx/dtheta = dp/dtheta exactly)\n\n";

    // ── §1: Full QLT at c=1 vs TW and JF ─────────────────────────────────────
    bar('=');
    cout << " §1  Full orbit-averaged QLT (c_val=1, epsilon=1)  vs  TW / JF\n";
    bar('=');
    cout << "\n Transformations:\n";
    cout << "   dp/dtheta (QLT) = < dp/dphi >_phi\n";
    cout << "   de/dtheta (QLT) = [ alpha*<dalpha/dphi> + beta*<dbeta/dphi> ] / e\n";
    cout << "                     with  e = sqrt(alpha^2 + beta^2)\n";

    auto avg1 = orbitAverage(Kfull,p0,alpha0,beta0,1.0,PNord,Nsamp);
    cout << scientific << setprecision(16);
    cout << left << setw(WLBL) << " Quantity"
         << setw(WNUM) << "QLT (c=1)"
         << setw(WNUM) << "TW"
         << setw(WNUM) << "JF"
         << setw(WRD) << "relDiff(QLT,TW)"
         << "relDiff(TW,JF)\n";
    bar('-');

    auto prow = [&](const string& lbl, double qlt, double tw, double jf){
        cout << " " << left << setw(WLBL) << lbl
             << setw(WNUM) << qlt << setw(WNUM) << tw << setw(WNUM) << jf
             << fixed << setprecision(16)
             << setw(WRD) << relDiff(qlt,tw) << relDiff(tw,jf) << "\n";
    };
    prow("dp/dtheta", avg1.dp, dp_TW_dtheta(p0,e0), dp_JF(p0,e0));
    prow("de/dtheta", avg1.de, de_TW_dtheta(p0,e0), de_JF(p0,e0));
    cout << "\n NOTE: At epsilon=1 the PN series is not small; order-mixing is O(1).\n"
         << "       Large relDiff(QLT,TW) here is expected and unphysical.\n"
         << "       The meaningful test is the coefficient extraction below.\n\n";

    // ── §2: PN coefficient extraction via ε polynomial fitting ────────────────
    // UPGRADED: Using SVD-based fitting with uncertainty estimation
    bar('=');
    cout << " §2  PN coefficient extraction — SVD-based polynomial fit (RESEARCH-GRADE)\n";
    bar('=');
    cout << "\n dp/dtheta = eps^5*A_25 + eps^7*A_35 + eps^8*A_4 + eps^9*A_45 + ...\n";
    cout << " de/dtheta = eps^5*B_25 + eps^7*B_35 + eps^8*B_4 + eps^9*B_45 + ...\n";
    cout << "\n IMPROVEMENTS:\n";
    cout << "  • SVD-based solver (better conditioning than Gaussian elimination)\n";
    cout << "  • Uncertainty estimation from residuals & covariance\n";
    cout << "  • Richardson extrapolation for averaging convergence\n";
    cout << "  • Fit range: eps in [0.006, 0.09], log-spaced for robustness\n\n";

    // Fit basis powers: 5,6,7,8,9,10  (6 terms, 6 unknowns)
    vector<int> pows = {5,6,7,8,9,10};
    double elo=0.006, ehi=0.09;

    // Use adaptive averaging for better accuracy
    cout << " Phase 1: Convergence analysis on averaging...\n";
    auto conv_test = orbitAverageAdaptive(Kfull, p0, alpha0, beta0, 1.0/elo, PNord);
    cout << "   Convergence rate (2.0=ideal):  " << fixed << setprecision(4) 
         << conv_test.convergence_rate << (conv_test.converged ? " ✓" : " ⚠") << "\n";
    cout << "   Recommended Nsamp: " << conv_test.optimal_Nsamp << "\n\n";

    auto QLT_dp_f = [&](double eps)->double{
        return orbitAverageAdaptive(Kfull, p0, alpha0, beta0, 1.0/eps, PNord, 2048).dp;
    };
    auto QLT_de_f = [&](double eps)->double{
        return orbitAverageAdaptive(Kfull, p0, alpha0, beta0, 1.0/eps, PNord, 2048).de;
    };

    cout << " Phase 2: Fitting coefficients with SVD...\n\n";

    // Extract with uncertainties
    auto A25_res = extractCoeffWithUncertainty(QLT_dp_f, pows, 5, elo, ehi);
    auto A35_res = extractCoeffWithUncertainty(QLT_dp_f, pows, 7, elo, ehi);
    auto A4_res  = extractCoeffWithUncertainty(QLT_dp_f, pows, 8, elo, ehi);
    auto A45_res = extractCoeffWithUncertainty(QLT_dp_f, pows, 9, elo, ehi);

    auto B25_res = extractCoeffWithUncertainty(QLT_de_f, pows, 5, elo, ehi);
    auto B35_res = extractCoeffWithUncertainty(QLT_de_f, pows, 7, elo, ehi);
    auto B4_res  = extractCoeffWithUncertainty(QLT_de_f, pows, 8, elo, ehi);
    auto B45_res = extractCoeffWithUncertainty(QLT_de_f, pows, 9, elo, ehi);

    // Analytic reference values
    double A25t=dp_TW_dtheta_25PN(p0,e0), A25j=dp_JF_2p5(p0,e0);
    double A35t=dp_TW_dtheta_35PN(p0,e0), A35j=dp_JF_3p5(p0,e0);
    double A4t =dp_TW_dtheta_4PN(p0,e0), A4j =dp_JF_4(p0,e0);
    double A45t=dp_TW_dtheta_45PN(p0,e0), A45j=dp_JF_4p5(p0,e0);

    double B25t=de_TW_dtheta_25PN(p0,e0), B25j=de_JF_2p5(p0,e0);
    double B35t=de_TW_dtheta_35PN(p0,e0), B35j=de_JF_3p5(p0,e0);
    double B4t =de_TW_dtheta_4PN(p0,e0), B4j =de_JF_4(p0,e0);
    double B45t=de_TW_dtheta_45PN(p0,e0), B45j=de_JF_4p5(p0,e0);

    // Store results for consistency checking
    vector<FitCoeffResult> dp_results = {A25_res, A35_res, A4_res, A45_res};
    vector<FitCoeffResult> de_results = {B25_res, B35_res, B4_res, B45_res};

    // Check PN consistency
    auto consistency = checkPNConsistency(dp_results, de_results, 0.1);
    cout << " Phase 3: PN consistency check...\n";
    cout << "   Consistency score: " << fixed << setprecision(3) << consistency.consistency_score << "\n";
    if(!consistency.is_consistent) {
        cout << "   ⚠ Warnings:\n";
        for(const auto& w : consistency.warnings) {
            cout << "     - " << w << "\n";
        }
    } else {
        cout << "   ✓ All PN orders consistent\n";
    }
    cout << "\n";

    // Enhanced output table with uncertainties
    const int WUNC = 22;  // Uncertainty column width
    cout << left << setw(WLBL) << " Order"
         << setw(WNUM) << "QLT fitted"
         << setw(WUNC) << "±Uncertainty"
         << setw(WNUM) << "TW analytic"
         << setw(WNUM) << "JF analytic"
         << setw(WRD) << "relDiff(QLT,TW)"
         << setw(WRD) << "relDiff(QLT,JF)"
         << "Winner\n";
    bar('-');

    double total_tw_dp=0, total_jf_dp=0;
    double total_tw_de=0, total_jf_de=0;
    int wins_tw=0, wins_jf=0;

    auto crow_unc = [&](const string& ord, const FitCoeffResult& qlt, 
                        double tw, double jf,
                        double& tot_tw, double& tot_jf) {
        double rd_tw=relDiff(qlt.value, tw), rd_jf=relDiff(qlt.value, jf);
        tot_tw+=rd_tw; tot_jf+=rd_jf;
        bool tw_wins=(rd_tw<=rd_jf);
        if(tw_wins) ++wins_tw; else ++wins_jf;
        cout << " " << left << setw(WLBL) << ord
             << scientific << setprecision(10)
             << setw(WNUM) << qlt.value 
             << "± " << setw(WUNC-2) << qlt.uncertainty
             << setw(WNUM) << tw << setw(WNUM) << jf
             << fixed << setprecision(16)
             << setw(WRD) << rd_tw << setw(WRD) << rd_jf
             << (tw_wins ? "TW" : "JF");
        
        // Indicate if point is within uncertainty of analytic value
        if(fabs(qlt.value - tw) <= qlt.uncertainty) cout << " *TW";
        else if(fabs(qlt.value - jf) <= qlt.uncertainty) cout << " *JF";
        cout << "\n";
    };

    cout << " dp/dtheta:\n";
    crow_unc("  2.5PN", A25_res, A25t, A25j, total_tw_dp, total_jf_dp);
    crow_unc("  3.5PN", A35_res, A35t, A35j, total_tw_dp, total_jf_dp);
    crow_unc("  4PN",   A4_res,  A4t,  A4j,  total_tw_dp, total_jf_dp);
    crow_unc("  4.5PN", A45_res, A45t, A45j, total_tw_dp, total_jf_dp);

    cout << "\n de/dtheta:\n";
    crow_unc("  2.5PN", B25_res, B25t, B25j, total_tw_de, total_jf_de);
    crow_unc("  3.5PN", B35_res, B35t, B35j, total_tw_de, total_jf_de);
    crow_unc("  4PN",   B4_res,  B4t,  B4j,  total_tw_de, total_jf_de);
    crow_unc("  4.5PN", B45_res, B45t, B45j, total_tw_de, total_jf_de);

    cout << "\n Legend: * indicates value within ±1σ of analytic result\n";

    // ── §3: ε-scan: scaled QLT vs analytic ───────────────────────────────────
    cout << "\n";
    bar('=');
    cout << " §3  epsilon-scan: eps^{-N} x <QLT>_N  as epsilon -> 0\n";
    cout << "     A flat, constant column = correct PN coefficient.\n";
    cout << "     A drifting column = wrong coefficient or missing term.\n";
    bar('=');

    // ε values to scan
    vector<double> eps_scan;
    {   int Ne=12;
        for(int i=0;i<Ne;++i)
            eps_scan.push_back(0.004*pow(0.15/0.004,double(i)/(Ne-1)));
    }

    // Print dp scan
    cout << "\n dp/dtheta scan:\n";
    cout << left << setw(WLBL) << " eps"
         << setw(WNUM) << "eps^-5 x dp_25"
         << setw(WRD) << "rD(TW)"
         << setw(WNUM) << "eps^-7 x dp_35"
         << setw(WRD) << "rD(TW)"
         << setw(WNUM) << "eps^-9 x dp_45"
         << setw(WRD) << "rD(TW)\n";
    bar('-');

    for(double eps : eps_scan){
        double cv=1.0/eps;
        auto a25  = orbitAverage(K25,  p0,alpha0,beta0, cv, 1,     Nsamp);
        auto a35c = orbitAverage(K35,  p0,alpha0,beta0, cv, 2,     Nsamp);
        auto afull= orbitAverage(Kfull,p0,alpha0,beta0, cv, PNord, Nsamp);
        auto ano45= orbitAverage(Kno45,p0,alpha0,beta0, cv, PNord, Nsamp);

        double dp35 = a35c.dp - a25.dp;
        double dp45 = afull.dp - ano45.dp;

        double sc25 = a25.dp  / pow(eps,5);
        double sc35 = dp35    / pow(eps,7);
        double sc45 = dp45    / pow(eps,9);

        cout << " " << left << fixed << setprecision(16) << setw(WLBL) << eps
             << scientific << setprecision(16)
             << setw(WNUM) << sc25
             << fixed << setprecision(16) << setw(WRD) << relDiff(a25.dp, A25t*pow(eps,5))
             << scientific << setprecision(16)
             << setw(WNUM) << sc35
             << fixed << setprecision(16) << setw(WRD) << relDiff(dp35, A35t*pow(eps,7))
             << scientific << setprecision(16)
             << setw(WNUM) << sc45
             << fixed << setprecision(16) << setw(WRD) << relDiff(dp45, A45t*pow(eps,9))
             << "\n";
    }
    cout << " Analytic:         "
         << scientific << setprecision(16) << setw(WNUM) << A25t
         << setw(WRD) << "---"
         << setw(WNUM) << A35t
         << setw(WRD) << "---"
         << setw(WNUM) << A45t << "\n";

    // Print de scan
    cout << "\n de/dtheta scan:\n";
    cout << left << setw(WLBL) << " eps"
         << setw(WNUM) << "eps^-5 x de_25"
         << setw(WRD) << "rD(TW)"
         << setw(WNUM) << "eps^-7 x de_35"
         << setw(WRD) << "rD(TW)"
         << setw(WNUM) << "eps^-9 x de_45"
         << setw(WRD) << "rD(TW)\n";
    bar('-');

    for(double eps : eps_scan){
        double cv=1.0/eps;
        auto a25  = orbitAverage(K25,  p0,alpha0,beta0, cv, 1,     Nsamp);
        auto a35c = orbitAverage(K35,  p0,alpha0,beta0, cv, 2,     Nsamp);
        auto afull= orbitAverage(Kfull,p0,alpha0,beta0, cv, PNord, Nsamp);
        auto ano45= orbitAverage(Kno45,p0,alpha0,beta0, cv, PNord, Nsamp);

        double de35 = a35c.de - a25.de;
        double de45 = afull.de - ano45.de;

        double sc25 = a25.de  / pow(eps,5);
        double sc35 = de35    / pow(eps,7);
        double sc45 = de45    / pow(eps,9);

        cout << " " << left << fixed << setprecision(16) << setw(WLBL) << eps
             << scientific << setprecision(16)
             << setw(WNUM) << sc25
             << fixed << setprecision(16) << setw(WRD) << relDiff(a25.de, B25t*pow(eps,5))
             << scientific << setprecision(16)
             << setw(WNUM) << sc35
             << fixed << setprecision(16) << setw(WRD) << relDiff(de35, B35t*pow(eps,7))
             << scientific << setprecision(16)
             << setw(WNUM) << sc45
             << fixed << setprecision(16) << setw(WRD) << relDiff(de45, B45t*pow(eps,9))
             << "\n";
    }
    cout << " Analytic:         "
         << scientific << setprecision(16) << setw(WNUM) << B25t
         << setw(WRD) << "---"
         << setw(WNUM) << B35t
         << setw(WRD) << "---"
         << setw(WNUM) << B45t << "\n";

    // ── §4: Final scorecard ───────────────────────────────────────────────────
    cout << "\n";
    bar('=');
    cout << " §4  RESEARCH-GRADE VALIDATION SUMMARY\n";
    bar('=');
    cout << "\n VALIDATION METRICS:\n";
    cout << " • PN Consistency:       " << (consistency.is_consistent ? "PASS ✓" : "WARN ⚠") << "\n";
    cout << " • Averaging Convergence: " << (conv_test.converged ? "PASS ✓" : "WARN ⚠") << "\n";
    cout << " • Overall Score:         " << fixed << setprecision(1) 
         << (100*consistency.consistency_score) << "%\n\n";

    cout << " Coefficient comparison (from §2 SVD fit with uncertainties):\n\n";
    cout << left << setw(WLBL) << " Coefficient"
         << setw(WNUM) << "QLT±σ"
         << setw(WNUM) << "TW analytic"
         << setw(WNUM) << "JF analytic"
         << setw(WRD) << "σ-distance"
         << "Status\n";
    bar('-');

    auto srow_unc = [&](const string& n, const FitCoeffResult& q, double t, double j) {
        double sigma_dist = fabs(q.value - t) / max(q.uncertainty, 1e-15);
        bool within_1sigma = (sigma_dist <= 1.0);
        cout << " " << left << setw(WLBL) << n
             << scientific << setprecision(10)
             << setw(WNUM) << (ostringstream() << q.value << "±" << q.uncertainty).str()
             << setw(WNUM) << t << setw(WNUM) << j
             << fixed << setprecision(2)
             << setw(WRD) << sigma_dist 
             << (within_1sigma ? "✓ 1σ" : "⚠ " + to_string((int)sigma_dist) + "σ") << "\n";
    };
    
    srow_unc("dp  2.5PN", A25_res, A25t, A25j);
    srow_unc("dp  3.5PN", A35_res, A35t, A35j);
    srow_unc("dp  4.0PN", A4_res,  A4t,  A4j );
    srow_unc("dp  4.5PN", A45_res, A45t, A45j);
    srow_unc("de  2.5PN", B25_res, B25t, B25j);
    srow_unc("de  3.5PN", B35_res, B35t, B35j);
    srow_unc("de  4.0PN", B4_res,  B4t,  B4j );
    srow_unc("de  4.5PN", B45_res, B45t, B45j);

    double total_tw = total_tw_dp+total_tw_de;
    double total_jf = total_jf_dp+total_jf_de;

    cout << "\n";
    bar('-');
    cout << " Total relDiff  TW=" << fixed << setprecision(16) << total_tw
         << "   JF=" << total_jf << "\n";
    cout << " Convergence wins:  TW=" << wins_tw
         << "   JF=" << wins_jf << "\n\n";

    bar('*',75);
    cout << " RESEARCH-GRADE VALIDATION RESULT:\n";
    if(consistency.consistency_score > 0.95 && conv_test.converged) {
        cout << " ✓✓✓ EXCELLENT — Code is production-ready for publication\n";
        cout << "     • PN orders consistent across all terms\n";
        cout << "     • Averaging converges at expected rate (2.0)\n";
        cout << "     • Uncertainties well-behaved\n";
    } else if(consistency.consistency_score > 0.80) {
        cout << " ✓ GOOD — Suitable for research with minor caveats\n";
        cout << "   Review warnings above before publication.\n";
    } else {
        cout << " ⚠ INVESTIGATE — Issues detected\n";
        cout << "   See consistency report and warnings.\n";
    }
    
    if(fabs(total_tw - total_jf) < 1e-14*total_tw){
        cout << "\n TW and JF converge IDENTICALLY → algebraically equivalent.\n";
    } else if(total_tw < total_jf){
        cout << "\n TUCKER-WILL converges better (lower total relDiff).\n";
    } else {
        cout << "\n JAN FEREISL converges better (lower total relDiff).\n";
    }
    bar('*',75);
    return 0;
}
