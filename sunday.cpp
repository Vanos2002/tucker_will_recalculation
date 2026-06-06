// THIS CODE WAS CREATED ON THURSDAY JUNE 4, 2026
// THIS FILE COMPARES THE FEREISL-TUCKERWILL 4.5PN EXPRESSIONS TO THE RK4 ADAPTIVE CODE

// From the article "Residual eccentricity of inspiralling orbits at the gravitational-wave detection threshold: Accurate estimates using post-Newtonian theory"
// by Alexandria Tucker and Clifford M. Will (arXiv:2108.12210v2 [gr-qc] 15 Nov 2021)
// We compare the transformed 4.5PN contributions of dp/dtheta and de/dtheta from the Fereisl-Tucker-Will (FTW) paper against the numerical orbit-averaged QLT results

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <functional>
#include <iomanip>
#include <string>
#include <algorithm>
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

    // Definition of the QLT equations of motion (see eq. 2.8 in the paper) 
    double dp_dphi = 2.0*r*r*r/(G*M)*ScS;
    double dalpha  = r*r/(G*M)*(  ScR*sin(phi) + ScS*(alpha+cos(phi))*(1.0+r/p) - ScS*alpha );
    double dbeta   = r*r/(G*M)*( -ScR*cos(phi) + ScS*(beta +sin(phi))*(1.0+r/p) - ScS*beta  );

    return {dp_dphi, dalpha, dbeta};
}




//-----------------------------------------------
// JAN FEREISL RESULTS FOR EVOLUTION
//-----------------------------------------------





// Physical constants and parameters
struct PhysicalParams {
    double G;      // Gravitational constant
    double M;      // Total mass
    double eta;    // Symmetric mass ratio (eta = m1*m2/(m1+m2)^2)
    double phi;    // Phase angle
    double eps;    // PN expansion parameter: ε = 1/c
};

// State variables for binary evolution
struct BinaryState {
    double p;        // Orbital parameter (semi-latus rectum)
    double alpha;    // Spin parameter 1
    double beta;     // Spin parameter 2
};

// Container for RHS terms: {dp/dtheta, dalpha/dtheta, dbeta/dtheta}
using SecularRHS = std::array<double, 3>;

std::array<double, 3> oscillatory_1PN(const BinaryState& state, const PhysicalParams& params);
std::array<double, 3> oscillatory_2PN(const BinaryState& state, const PhysicalParams& params);

static constexpr int PHI_AVERAGE_SAMPLES = 512;
using QLTrhsFunc = std::function<QLTrhs(const BinaryState&)>;

static QLTrhs addQLT(const QLTrhs& a, const QLTrhs& b) {
    return {a.dp + b.dp, a.dalpha + b.dalpha, a.dbeta + b.dbeta};
}

static QLTrhs subQLT(const QLTrhs& a, const QLTrhs& b) {
    return {a.dp - b.dp, a.dalpha - b.dalpha, a.dbeta - b.dbeta};
}

static QLTrhs scaleQLT(const QLTrhs& a, double s) {
    return {a.dp * s, a.dalpha * s, a.dbeta * s};
}

static QLTrhs computeQLTDissipativeOrder(const PNCoeffs& K,
                                        const BinaryState& state,
                                        const PhysicalParams& params,
                                        int dissipativeOrder)
{
    double p = state.p;
    double alpha = state.alpha;
    double beta = state.beta;
    double phi = params.phi;
    double c_val = 1.0 / params.eps;
    double e = std::sqrt(alpha*alpha + beta*beta);
    double r = normR(p, e, alpha, beta, phi);
    double rd = rDot(p, e, alpha, beta, phi);
    double v2 = normV2(p, e, alpha, beta, phi);
    double rd2 = rd * rd;
    double GM = params.G * params.M;
    double eta = params.eta;
    double gmr = GM / r;

    double Crr = sumTable(K.c, dissipativeOrder, rd2, v2, gmr, c_val);
    double Drr = sumTable(K.d, dissipativeOrder, rd2, v2, gmr, c_val);

    double eps = params.eps;
    double eps3 = eps * eps * eps;
    double rr_R_pre = (8.0 / 5.0) * eta * eps3 * GM * GM / std::pow(r, 3) * rd;
    double rr_S_pre = (8.0 / 5.0) * eta * eps3 * std::pow(GM / (r * r), 2) * std::sqrt(GM * p);

    double ScR = rr_R_pre * (Crr + Drr);
    double ScS = rr_S_pre * Drr;

    double dp_dphi = 2.0 * r * r * r / GM * ScS;
    double dalpha  = r * r / GM * (ScR * std::sin(phi) + ScS * (alpha + std::cos(phi)) * (1.0 + r / p) - ScS * alpha);
    double dbeta   = r * r / GM * (-ScR * std::cos(phi) + ScS * (beta + std::sin(phi)) * (1.0 + r / p) - ScS * beta);

    return {dp_dphi, dalpha, dbeta};
}
static QLTrhs averageQLTRHS(const std::function<QLTrhs(double)>& f, int samples = PHI_AVERAGE_SAMPLES) {
    QLTrhs sum{0.0, 0.0, 0.0};
    double h = 2.0 * PI / samples;
    for (int i = 0; i < samples; ++i) {
        double phi = (i + 0.5) * h;
        QLTrhs val = f(phi);
        sum.dp += val.dp;
        sum.dalpha += val.dalpha;
        sum.dbeta += val.dbeta;
    }
    double scale = h / (2.0 * PI);
    return scaleQLT(sum, scale);
}

static SecularRHS secular_2_5PN_explicit(const PNCoeffs& K,
                                          const BinaryState& state,
                                          const PhysicalParams& params)
{
    auto avg = averageQLTRHS([&](double phi) {
        PhysicalParams pphi = params;
        pphi.phi = phi;
        return computeQLTDissipativeOrder(K, state, pphi, 1);
    });
    return {avg.dp, avg.dalpha, avg.dbeta};
}

static SecularRHS secular_3_5PN_explicit(const PNCoeffs& K,
                                          const BinaryState& state,
                                          const PhysicalParams& params)
{
    SecularRHS result{0.0, 0.0, 0.0};
    double h = 2.0 * PI / PHI_AVERAGE_SAMPLES;
    for (int i = 0; i < PHI_AVERAGE_SAMPLES; ++i) {
        double phi = (i + 0.5) * h;
        PhysicalParams pphi = params;
        pphi.phi = phi;

        // Only include the averaged 3.5PN dissipative contribution (no analytic corrections)
        QLTrhs F35 = computeQLTDissipativeOrder(K, state, pphi, 2);
        result[0] += F35.dp;
        result[1] += F35.dalpha;
        result[2] += F35.dbeta;
    }
    double scale = h / (2.0 * PI);
    result[0] *= scale;
    result[1] *= scale;
    result[2] *= scale;
    return result;
}

static SecularRHS secular_4_5PN_explicit(const PNCoeffs& K,
                                          const BinaryState& state,
                                          const PhysicalParams& params)
{
    SecularRHS result{0.0, 0.0, 0.0};
    double h = 2.0 * PI / PHI_AVERAGE_SAMPLES;
    for (int i = 0; i < PHI_AVERAGE_SAMPLES; ++i) {
        double phi = (i + 0.5) * h;
        PhysicalParams pphi = params;
        pphi.phi = phi;

        // Only include the averaged 4.5PN dissipative contribution (no analytic corrections)
        QLTrhs F45 = computeQLTDissipativeOrder(K, state, pphi, 3);
        result[0] += F45.dp;
        result[1] += F45.dalpha;
        result[2] += F45.dbeta;
    }
    double scale = h / (2.0 * PI);
    result[0] *= scale;
    result[1] *= scale;
    result[2] *= scale;
    return result;
}


SecularRHS secular_1PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    double eps2 = params.eps * params.eps;
    
    SecularRHS rhs;
    rhs[0] = 0.0;
    rhs[1] = eps2 * (-(3.0 * GM * beta) / p);
    rhs[2] = eps2 * ((3.0 * GM * alpha) / p);
    
    return rhs;
}

std::array<double, 3> oscillatory_1PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double& phi = params.phi;
    const double GM = params.G * params.M;
    
    std::array<double, 3> Y;
    
    // Y[p]
    Y[0] = 4.0 * GM * (-2.0 + eta) * (alpha * std::cos(phi) + beta * std::sin(phi));
    
    // Y[alpha]
    double cos_phi = std::cos(phi);
    double cos_2phi = std::cos(2.0 * phi);
    double cos_3phi = std::cos(3.0 * phi);
    double sin_phi = std::sin(phi);
    double sin_2phi = std::sin(2.0 * phi);
    double sin_3phi = std::sin(3.0 * phi);
    
    Y[1] = (1.0 / (8.0 * p)) * GM * (
        (8.0 * (-3.0 + eta) + beta * beta * (8.0 + 21.0 * eta) + 
         alpha * alpha * (-56.0 + 47.0 * eta)) * cos_phi +
        4.0 * alpha * (-5.0 + 4.0 * eta) * cos_2phi +
        (alpha * alpha - beta * beta) * eta * cos_3phi +
        2.0 * alpha * beta * (-32.0 + 13.0 * eta) * sin_phi +
        4.0 * beta * (-5.0 + 4.0 * eta) * sin_2phi +
        2.0 * alpha * beta * eta * sin_3phi
    );
    
    // Y[beta]
    Y[2] = (1.0 / (8.0 * p)) * GM * (
        2.0 * alpha * beta * (-32.0 + 13.0 * eta) * cos_phi +
        4.0 * beta * (5.0 - 4.0 * eta) * cos_2phi -
        2.0 * alpha * beta * eta * cos_3phi +
        (8.0 * (-3.0 + eta) + alpha * alpha * (8.0 + 21.0 * eta) + 
         beta * beta * (-56.0 + 47.0 * eta)) * sin_phi +
        4.0 * alpha * (-5.0 + 4.0 * eta) * sin_2phi +
        (alpha * alpha - beta * beta) * eta * sin_3phi
    );
    
    return Y;
}

// ============================================================================
// 2 PN ORDER TERMS
// ============================================================================

SecularRHS secular_2PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM2 = GM * GM;
    double eps2 = params.eps * params.eps;
    double eps4 = eps2 * eps2;
    
    SecularRHS rhs;
    
    double bracket = -10.0 + beta * beta - 4.0 * eta + 10.0 * beta * beta * eta + 
                     alpha * alpha * (1.0 + 10.0 * eta);
    
    rhs[0] = 0.0;
    rhs[1] = eps4 * (-(3.0 * GM2 * beta * bracket) / (4.0 * p * p));
    rhs[2] = eps4 * ((3.0 * GM2 * alpha * bracket) / (4.0 * p * p));
    
    return rhs;
}

std::array<double, 3> oscillatory_2PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double& phi = params.phi;
    const double GM = params.G * params.M;
    const double GM2 = GM * GM;
    
    std::array<double, 3> Y;
    
    double cos_phi = std::cos(phi);
    double cos_2phi = std::cos(2.0 * phi);
    double cos_3phi = std::cos(3.0 * phi);
    double cos_4phi = std::cos(4.0 * phi);
    double cos_5phi = std::cos(5.0 * phi);
    double sin_phi = std::sin(phi);
    double sin_2phi = std::sin(2.0 * phi);
    double sin_3phi = std::sin(3.0 * phi);
    double sin_4phi = std::sin(4.0 * phi);
    double sin_5phi = std::sin(5.0 * phi);
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    
    // Y[p]
    Y[0] = (1.0 / (4.0 * p)) * GM2 * (
        alpha * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2 + 
                 2.0 * (-8.0 + a2 + b2) * eta)) * cos_phi +
        (a2 - b2) * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * cos_2phi -
        alpha * (a2 - 3.0 * b2) * eta * (3.0 + 2.0 * eta) * cos_3phi +
        beta * (-160.0 + eta * (256.0 - 33.0 * a2 - 33.0 * b2 + 
                2.0 * (-8.0 + a2 + b2) * eta)) * sin_phi +
        2.0 * alpha * beta * (68.0 + 3.0 * eta * (-15.0 + 4.0 * eta)) * sin_2phi +
        beta * (-3.0 * a2 + b2) * eta * (3.0 + 2.0 * eta) * sin_3phi
    );
    
    // Y[alpha]
    Y[1] = (1.0 / (128.0 * p * p)) * GM2 * (
        -2.0 * (5.0 * b4 * eta * (-27.0 + 41.0 * eta) + 
                a4 * eta * (477.0 + 161.0 * eta) + 
                b2 * (-944.0 + 92.0 * eta + 64.0 * eta * eta) - 
                16.0 * (60.0 + eta * (17.0 + 8.0 * eta)) + 
                a2 * (3056.0 + 2.0 * eta * (-1958.0 + 224.0 * eta + 
                      3.0 * b2 * (57.0 + 61.0 * eta)))) * cos_phi +
        16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta) + 
                       b2 * (-162.0 + eta * (57.0 + 28.0 * eta)) + 
                       a2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi +
        (-a4 * eta * (73.0 + 53.0 * eta) + 
         b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta) + b2 * eta * (-33.0 + 19.0 * eta)) + 
         2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta + 3.0 * b2 * (53.0 + 17.0 * eta)))) * cos_3phi -
        8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi -
        3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * cos_5phi +
        8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta + 
                            a2 * (-153.0 + 11.0 * eta) + 
                            b2 * (-153.0 + 11.0 * eta))) * sin_phi +
        16.0 * beta * (-28.0 * (1.0 + b2 - 5.0 * eta) + 
                      b2 * eta * (-53.0 + 44.0 * eta) + 
                      a2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi -
        4.0 * alpha * beta * (-240.0 + eta * (b2 * (-43.0 + eta) + 
                            7.0 * a2 * (9.0 + 5.0 * eta) - 
                            4.0 * (9.0 + 16.0 * eta))) * sin_3phi +
        8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi -
        12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    );
    
    // Y[beta]
    Y[2] = (1.0 / (128.0 * p * p)) * GM2 * (
        8.0 * alpha * beta * (-1000.0 + eta * (1002.0 - 96.0 * eta + 
                            a2 * (-153.0 + 11.0 * eta) + 
                            b2 * (-153.0 + 11.0 * eta))) * cos_phi -
        16.0 * beta * (28.0 * (-1.0 + 5.0 * eta) + 
                      a2 * (-162.0 + eta * (57.0 + 28.0 * eta)) + 
                      b2 * (106.0 + eta * (-163.0 + 60.0 * eta))) * cos_2phi +
        4.0 * alpha * beta * (-240.0 + eta * (a2 * (-43.0 + eta) + 
                            7.0 * b2 * (9.0 + 5.0 * eta) - 
                            4.0 * (9.0 + 16.0 * eta))) * cos_3phi -
        8.0 * beta * (-3.0 * a2 + b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * cos_4phi +
        12.0 * alpha * (alpha - beta) * beta * (alpha + beta) * eta * (-1.0 + 3.0 * eta) * cos_5phi -
        2.0 * (5.0 * a4 * eta * (-27.0 + 41.0 * eta) + 
               b4 * eta * (477.0 + 161.0 * eta) - 
               16.0 * (60.0 + eta * (17.0 + 8.0 * eta)) + 
               4.0 * b2 * (764.0 + eta * (-979.0 + 112.0 * eta)) + 
               a2 * (-944.0 + 2.0 * eta * (46.0 + 32.0 * eta + 
                     3.0 * b2 * (57.0 + 61.0 * eta)))) * sin_phi +
        16.0 * alpha * (28.0 * (-1.0 + 5.0 * eta) + 
                       a2 * (-28.0 + eta * (-53.0 + 44.0 * eta)) + 
                       b2 * (240.0 + eta * (-273.0 + 76.0 * eta))) * sin_2phi +
        (a4 * (33.0 - 19.0 * eta) * eta + 
         b2 * (-480.0 - 8.0 * eta * (9.0 + 16.0 * eta) + b2 * eta * (73.0 + 53.0 * eta)) + 
         2.0 * a2 * (240.0 + eta * (36.0 + 64.0 * eta - 3.0 * b2 * (53.0 + 17.0 * eta)))) * sin_3phi -
        8.0 * alpha * (a2 - 3.0 * b2) * (-1.0 + eta * (-5.0 + 4.0 * eta)) * sin_4phi -
        3.0 * (a4 - 6.0 * a2 * b2 + b4) * eta * (-1.0 + 3.0 * eta) * sin_5phi
    );
    
    return Y;
}

// ============================================================================
// 2.5 PN ORDER TERMS (RADIATION REACTION)
// ============================================================================

SecularRHS secular_2_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    
    rhs[0] = -(8.0 * GMp_5_2 * (8.0 + 7.0 * a2 + 7.0 * b2) * eta) / (5.0 * p * p * p * p);
    rhs[1] = -(GMp_5_2 * alpha * (304.0 + 121.0 * a2 + 121.0 * b2) * eta) / (15.0 * p * p * p * p * p);
    rhs[2] = -(GMp_5_2 * beta * (304.0 + 121.0 * a2 + 121.0 * b2) * eta) / (15.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// 3.5 PN ORDER TERMS
// ============================================================================

SecularRHS secular_3_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM3 = GM * GM * GM;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_3_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    
    // dp/dtheta
    double term_p = -8.0 * (2759.0 + 252.0 * eta) + 
                    8.0 * b2 * (758.0 + 889.0 * eta) + 
                    a4 * (1483.0 + 4424.0 * eta) + 
                    b4 * (1483.0 + 4424.0 * eta) + 
                    a2 * (8.0 * (758.0 + 889.0 * eta) + b2 * (2966.0 + 8848.0 * eta));
    
    rhs[0] = -(GM3 * GMp_3_2 * eta * term_p) / (210.0 * p * p * p);
    
    // dalpha/dtheta
    double term_alpha = -8.0 * (18049.0 + 4452.0 * eta) + 
                       4.0 * b2 * (8692.0 + 12803.0 * eta) + 
                       a4 * (2251.0 + 15064.0 * eta) + 
                       b4 * (2251.0 + 15064.0 * eta) + 
                       a2 * (34768.0 + 51212.0 * eta + b2 * (4502.0 + 30128.0 * eta));
    
    rhs[1] = -(GM3 * GMp_3_2 * alpha * eta * term_alpha) / (840.0 * p * p * p * p);
    
    // dbeta/dtheta
    double term_beta = -8.0 * (18049.0 + 4452.0 * eta) + 
                      4.0 * b2 * (8692.0 + 12803.0 * eta) + 
                      a4 * (2251.0 + 15064.0 * eta) + 
                      b4 * (2251.0 + 15064.0 * eta) + 
                      a2 * (34768.0 + 51212.0 * eta + b2 * (4502.0 + 30128.0 * eta));
    
    rhs[2] = -(GM3 * GMp_3_2 * beta * eta * term_beta) / (840.0 * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// 4.5 PN ORDER TERMS
// ============================================================================

SecularRHS secular_4_5PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM4 = GM * GM * GM * GM;
    

    // Tucker-Will Results beyoond
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    double a6 = a4 * a2;
    double b6 = b4 * b2;
    
    double eta2 = eta * eta;
    
    // dp/dtheta
    double term_p = 8.0 * (-1034075.0 - 261369.0 * eta + 36288.0 * eta2) + 
                   9.0 * a6 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                   9.0 * b6 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                   12.0 * b2 * (-64831.0 + 1186209.0 * eta + 64575.0 * eta2) + 
                   b4 * (947991.0 - 5386365.0 * eta + 3988656.0 * eta2) + 
                   3.0 * a4 * (315997.0 - 1795455.0 * eta + 1329552.0 * eta2 + 
                              9.0 * b2 * (527.0 - 6300.0 * eta + 53088.0 * eta2)) + 
                   3.0 * a2 * (9.0 * b4 * (527.0 - 6300.0 * eta + 53088.0 * eta2) + 
                              4.0 * (-64831.0 + 1186209.0 * eta + 64575.0 * eta2) + 
                              b2 * (631994.0 - 3590910.0 * eta + 2659104.0 * eta2));
    
    rhs[0] = (GM4 * GMp_5_2 * eta * term_p) / (11340.0 * p * p * p * p);
    
    // dalpha/dtheta
    double term_alpha = 16.0 * (-2739835.0 - 1394559.0 * eta + 145152.0 * eta2) + 
                       3.0 * a6 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                       3.0 * b6 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                       12.0 * b2 * (-354911.0 + 4848903.0 * eta + 511413.0 * eta2) + 
                       2.0 * b4 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2) + 
                       a4 * (9.0 * b2 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                             2.0 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2)) + 
                       a2 * (9.0 * b4 * (-25845.0 - 78380.0 * eta + 361536.0 * eta2) + 
                             12.0 * (-354911.0 + 4848903.0 * eta + 511413.0 * eta2) + 
                             4.0 * b2 * (605645.0 - 8079297.0 * eta + 5758704.0 * eta2));
    
    rhs[1] = (GM4 * GMp_5_2 * alpha * eta * term_alpha) / (30240.0 * p * p * p * p * p);
    
    // dbeta/dtheta (same as dalpha/dtheta but multiplied by beta instead of alpha)
    rhs[2] = (GM4 * GMp_5_2 * beta * eta * term_alpha) / (30240.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// COMPOSITE FUNCTIONS - COMBINING ALL PN ORDERS
// ============================================================================

SecularRHS compute_secular_RHS(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    static const PNCoeffs K = buildCoefficients(true);
    SecularRHS total_rhs = {0.0, 0.0, 0.0};
    
    if (max_PN_order >= 1) {
        auto rhs = secular_1PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 2) {
        auto rhs = secular_2PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 3) {  // 2.5 PN
        auto rhs = secular_2_5PN_explicit(K, state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 4) {  // 3.5 PN
        auto rhs = secular_3_5PN_explicit(K, state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 5) {  // 4.5 PN
        auto rhs = secular_4_5PN_explicit(K, state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    return total_rhs;
}




//----------------------------------------------------------------
// TUCKER-WILL RESULTS FOR EVOLUTION (ONLY DIFFER AT 4.5 PN ORDER)
//----------------------------------------------------------------

// Tucker-Will orders 1PN through 3.5PN are identical to Fereisl - use aliases
inline SecularRHS secular_1PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_1PN(state, params);
}

inline std::array<double, 3> oscillatory_1PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return oscillatory_1PN(state, params);
}

inline SecularRHS secular_2PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_2PN(state, params);
}

inline std::array<double, 3> oscillatory_2PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return oscillatory_2PN(state, params);
}

inline SecularRHS secular_2_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_2_5PN(state, params);
}

inline SecularRHS secular_3_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    return secular_3_5PN(state, params);
}

// ============================================================================
// 4.5 PN ORDER TERMS (Tucker-Will) - Only difference from Fereisl
// ============================================================================

SecularRHS secular_4_5PN_TW(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    const double GM4 = GM * GM * GM * GM;
    
    SecularRHS rhs;
    
    double sqrt_GMp = std::sqrt(GM * p);
    double GMp_5_2 = sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp * sqrt_GMp;
    
    double a2 = alpha * alpha;
    double b2 = beta * beta;
    double a4 = a2 * a2;
    double b4 = b2 * b2;
    double a6 = a4 * a2;
    double b6 = b4 * b2;
    
    double eta2 = eta * eta;
    
    // dp/dtheta (Tucker-Will 4.5PN)
    double term_p = -8272600.0 
                  + 72.0 * eta * (-29041.0 + 4032.0 * eta)
                  + 9.0 * a6 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta))
                  + 9.0 * b6 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta))
                  + 12.0 * b2 * (-64831.0 + 9.0 * eta * (131801.0 + 7175.0 * eta))
                  + b4 * (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta))
                  + a4 * (9.0 * b2 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta)) + 
                         (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta)))
                  + a2 * (12.0 * (-64831.0 + 9.0 * eta * (131801.0 + 7175.0 * eta)) + 
                         9.0 * b4 * (527.0 + 84.0 * eta * (-75.0 + 632.0 * eta)) + 
                         b2 * (2.0 * (947991.0 + 27.0 * eta * (-199495.0 + 147728.0 * eta))));
    
    rhs[0] = (GM4 * GMp_5_2 * eta * term_p) / (11340.0 * p * p * p * p);
    
    // dalpha/dtheta (Tucker-Will 4.5PN)
    double term_alpha = 43837360.0 
                      + 144.0 * eta * (154951.0 - 16128.0 * eta)
                      + 9.0 * a6 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta))
                      + 9.0 * b6 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta))
                      - 12.0 * b2 * (-354911.0 + 4848303.0 * eta + 511413.0 * eta * eta)
                      - 2.0 * b4 * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta))
                      + a4 * (9.0 * b2 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta)) + 
                             (-2.0) * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta)))
                      + a2 * (9.0 * b4 * (8615.0 + 4.0 * eta * (6565.0 - 30128.0 * eta)) + 
                             (-12.0) * (-354911.0 + 4848303.0 * eta + 511413.0 * eta * eta) + 
                             (-4.0) * b2 * (605645.0 + 9.0 * eta * (-898433.0 + 639856.0 * eta)));
    
    rhs[1] = -(alpha * eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p * p * p * p * p);
    
    // dbeta/dtheta (Tucker-Will 4.5PN - same as dalpha but with beta)
    rhs[2] = -(beta * eta * term_alpha * GM4 * GMp_5_2) / (30240.0 * p * p * p * p * p);
    
    return rhs;
}

// ============================================================================
// COMPOSITE FUNCTIONS - COMBINING ALL PN ORDERS (Tucker-Will)
// ============================================================================

SecularRHS compute_secular_RHS_TW(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    // Tucker-Will equations are identical to Fereisl for orders 1-3.5PN
    SecularRHS total_rhs = compute_secular_RHS(state, params, std::min(max_PN_order, 4));
    
    // Replace 4.5PN order with Tucker-Will version
    if (max_PN_order >= 5) {
        auto rhs_TW_4_5 = secular_4_5PN_TW(state, params);
        total_rhs[0] += rhs_TW_4_5[0];
        total_rhs[1] += rhs_TW_4_5[1];
            total_rhs[2] += rhs_TW_4_5[2];
    }
    
    return total_rhs;
}

static SecularRHS scaleSecularRHSForPhi(const SecularRHS& rhs, double eps) {
    double factor = (eps != 0.0) ? 1.0 / eps : 0.0;
    return {rhs[0] * factor, rhs[1] * factor, rhs[2] * factor};
}

SecularRHS compute_secular_RHS_phi(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    return scaleSecularRHSForPhi(compute_secular_RHS(state, params, max_PN_order), params.eps);
}

SecularRHS compute_secular_RHS_TW_phi(const BinaryState& state, const PhysicalParams& params, int max_PN_order) {
    return scaleSecularRHSForPhi(compute_secular_RHS_TW(state, params, max_PN_order), params.eps);
}

BinaryState transformTildeStateToActual(const BinaryState& tilde_state, const PhysicalParams& params, bool useTW = false) {
    double eps2 = params.eps * params.eps;
    double eps4 = eps2 * eps2;

    auto Y2 = useTW ? oscillatory_1PN_TW(tilde_state, params) : oscillatory_1PN(tilde_state, params);
    auto Y4 = useTW ? oscillatory_2PN_TW(tilde_state, params) : oscillatory_2PN(tilde_state, params);

    return {
        tilde_state.p + eps2 * Y2[0] + eps4 * Y4[0],
        tilde_state.alpha + eps2 * Y2[1] + eps4 * Y4[1],
        tilde_state.beta + eps2 * Y2[2] + eps4 * Y4[2]
    };
}

// ============================================================================
// IMPLICIT GAUSS COLLOCATION INTEGRATOR WITH ERROR ESTIMATION
// ============================================================================

struct IntegrationResult {
    std::vector<double> theta;
    std::vector<BinaryState> states;
    std::vector<double> error_estimates;
    int num_steps;
};

class AdaptiveGaussCollocationIntegrator {
private:
    double tolerance;
    double min_step;
    double max_step;
    int max_steps;
    
    typedef std::function<SecularRHS(const BinaryState&, const PhysicalParams&)> RHSFunc;
    
    static constexpr double sqrt3 = 1.7320508075688772;
    static constexpr double c1 = 0.5 - sqrt3 / 6.0;
    static constexpr double c2 = 0.5 + sqrt3 / 6.0;
    static constexpr double a11 = 0.25;
    static constexpr double a12 = 0.25 - sqrt3 / 6.0;
    static constexpr double a21 = 0.25 + sqrt3 / 6.0;
    static constexpr double a22 = 0.25;
    static constexpr double b1  = 0.5;
    static constexpr double b2  = 0.5;

public:
    AdaptiveGaussCollocationIntegrator(double tol = 1e-8, double h_min = 1e-5, double h_max = 0.1, int max = 100000)
        : tolerance(tol), min_step(h_min), max_step(h_max), max_steps(max) {}
    
    std::pair<BinaryState, double> adaptiveStep(
        const BinaryState& state, 
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double h
    ) {
        BinaryState y_full = gaussStep(state, params, rhs_func, h);
        BinaryState y_half_1 = gaussStep(state, params, rhs_func, h / 2.0);
        BinaryState y_half_2 = gaussStep(y_half_1, params, rhs_func, h / 2.0);

        double max_rel_error = 0.0;
        if (std::abs(y_half_2.p) > 1e-30) {
            max_rel_error = std::max(max_rel_error, std::abs(y_full.p - y_half_2.p) / std::abs(y_half_2.p));
        }
        if (std::abs(y_half_2.alpha) > 1e-30) {
            max_rel_error = std::max(max_rel_error, std::abs(y_full.alpha - y_half_2.alpha) / std::abs(y_half_2.alpha));
        }
        if (std::abs(y_half_2.beta) > 1e-30) {
            max_rel_error = std::max(max_rel_error, std::abs(y_full.beta - y_half_2.beta) / std::abs(y_half_2.beta));
        }

        return {y_half_2, max_rel_error};
    }
    
    IntegrationResult integrate(
        const BinaryState& initial_state,
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double theta_start,
        double theta_end
    ) {
        IntegrationResult result;
        BinaryState current_state = initial_state;
        double current_theta = theta_start;
        double h = (theta_end > theta_start) ? max_step : -max_step;
        int step_count = 0;
        
        result.theta.push_back(current_theta);
        result.states.push_back(current_state);
        result.error_estimates.push_back(0.0);
        
        while ((h > 0 && current_theta < theta_end) || (h < 0 && current_theta > theta_end)) {
            if (step_count >= max_steps) break;
            
            double h_attempt = h;
            if (h > 0 && current_theta + h > theta_end) {
                h_attempt = theta_end - current_theta;
            } else if (h < 0 && current_theta + h < theta_end) {
                h_attempt = theta_end - current_theta;
            }
            
            auto step_result = adaptiveStep(current_state, params, rhs_func, h_attempt);
            BinaryState next_state = step_result.first;
            double error = step_result.second;
            
            if (error < tolerance) {
                current_state = next_state;
                current_theta += h_attempt;
                result.theta.push_back(current_theta);
                result.states.push_back(current_state);
                result.error_estimates.push_back(error);
                step_count++;
                
                double scale = 0.9 * std::pow(tolerance / (error + 1e-16), 0.25);
                if (scale < 0.1) scale = 0.1;
                if (scale > 5.0) scale = 5.0;
                h *= scale;
                h = std::copysign(std::min(std::abs(h), max_step), h);
            } else {
                h *= 0.5;
                if (std::abs(h) < min_step) h = std::copysign(min_step, h);
            }
        }
        
        result.num_steps = step_count;
        return result;
    }
    
private:
    BinaryState gaussStep(
        const BinaryState& state,
        const PhysicalParams& params,
        const RHSFunc& rhs_func,
        double h
    ) {
        std::array<SecularRHS, 2> K{{rhs_func(state, params), rhs_func(state, params)}};
        std::array<BinaryState, 2> stage_states;
        const int max_iter = 20;
        const double stage_tol = 1e-10;

        for (int iter = 0; iter < max_iter; ++iter) {
            stage_states[0].p = state.p + h * (a11 * K[0][0] + a12 * K[1][0]);
            stage_states[0].alpha = state.alpha + h * (a11 * K[0][1] + a12 * K[1][1]);
            stage_states[0].beta = state.beta + h * (a11 * K[0][2] + a12 * K[1][2]);

            stage_states[1].p = state.p + h * (a21 * K[0][0] + a22 * K[1][0]);
            stage_states[1].alpha = state.alpha + h * (a21 * K[0][1] + a22 * K[1][1]);
            stage_states[1].beta = state.beta + h * (a21 * K[0][2] + a22 * K[1][2]);

            std::array<SecularRHS, 2> K_new{{rhs_func(stage_states[0], params), rhs_func(stage_states[1], params)}};

            double max_diff = 0.0;
            for (int s = 0; s < 2; ++s) {
                max_diff = std::max(max_diff, std::abs(K_new[s][0] - K[s][0]));
                max_diff = std::max(max_diff, std::abs(K_new[s][1] - K[s][1]));
                max_diff = std::max(max_diff, std::abs(K_new[s][2] - K[s][2]));
            }
            K = K_new;
            if (max_diff < stage_tol) break;
        }

        BinaryState next_state;
        next_state.p = state.p + h * (b1 * K[0][0] + b2 * K[1][0]);
        next_state.alpha = state.alpha + h * (b1 * K[0][1] + b2 * K[1][1]);
        next_state.beta = state.beta + h * (b1 * K[0][2] + b2 * K[1][2]);
        return next_state;
    }
};

// ============================================================================
// INTEGRATION COMPARISON
// ============================================================================

void compareEvolutionMethods(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance,
    const std::string& output_file
) {
    AdaptiveGaussCollocationIntegrator integrator(tolerance);
    
    auto rhs_fereisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS(s, p, max_PN_order);
    };
    auto result_fereisl = integrator.integrate(initial_state, params, rhs_fereisl, theta_start, theta_end);
    
    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS_TW(s, p, max_PN_order);
    };
    auto result_tw = integrator.integrate(initial_state, params, rhs_tw, theta_start, theta_end);
    
    std::ofstream outfile(output_file);
    outfile << std::scientific << std::setprecision(10);
    outfile << "theta,p_fereisl,alpha_fereisl,beta_fereisl,p_tw,alpha_tw,beta_tw,"
            << "dp_diff,dalpha_diff,dbeta_diff,max_diff\n";
    
    size_t min_size = std::min(result_fereisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        double theta = result_fereisl.theta[i];
        auto& state_f = result_fereisl.states[i];
        auto& state_tw_i = result_tw.states[i];
        
        double dp_diff = std::abs(state_f.p - state_tw_i.p);
        double da_diff = std::abs(state_f.alpha - state_tw_i.alpha);
        double db_diff = std::abs(state_f.beta - state_tw_i.beta);
        double max_diff = std::max({dp_diff, da_diff, db_diff});
        
        outfile << theta << "," 
                << state_f.p << "," << state_f.alpha << "," << state_f.beta << ","
                << state_tw_i.p << "," << state_tw_i.alpha << "," << state_tw_i.beta << ","
                << dp_diff << "," << da_diff << "," << db_diff << "," << max_diff << "\n";
    }
    outfile.close();
    
    std::cout << "\n=== CONVERGENCE COMPARISON ===" << std::endl;
    std::cout << "Fereisl steps: " << result_fereisl.num_steps << std::endl;
    std::cout << "Tucker-Will steps: " << result_tw.num_steps << std::endl;
    std::cout << "Results written to: " << output_file << std::endl;
}

void compareEvolutionMethodsPhi(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double phi_start,
    double phi_end,
    double tolerance,
    const std::string& output_file
) {
    AdaptiveGaussCollocationIntegrator integrator(tolerance);

    auto rhs_fereisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS_phi(s, p, max_PN_order);
    };
    auto result_fereisl = integrator.integrate(initial_state, params, rhs_fereisl, phi_start, phi_end);

    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS_TW_phi(s, p, max_PN_order);
    };
    auto result_tw = integrator.integrate(initial_state, params, rhs_tw, phi_start, phi_end);

    std::ofstream outfile(output_file);
    outfile << std::scientific << std::setprecision(10);
    outfile << "phi,p_fereisl_actual,alpha_fereisl_actual,beta_fereisl_actual,"
            << "p_tw_actual,alpha_tw_actual,beta_tw_actual,dp_diff,dalpha_diff,dbeta_diff,max_diff\n";

    size_t min_size = std::min(result_fereisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        double phi = result_fereisl.theta[i];
        auto& state_f = result_fereisl.states[i];
        auto& state_tw_i = result_tw.states[i];

        PhysicalParams pphi = params;
        pphi.phi = phi;

        auto actual_f = transformTildeStateToActual(state_f, pphi, false);
        auto actual_tw = transformTildeStateToActual(state_tw_i, pphi, true);

        double dp_diff = std::abs(actual_f.p - actual_tw.p);
        double da_diff = std::abs(actual_f.alpha - actual_tw.alpha);
        double db_diff = std::abs(actual_f.beta - actual_tw.beta);
        double max_diff = std::max({dp_diff, da_diff, db_diff});

        outfile << phi << ","
                << actual_f.p << "," << actual_f.alpha << "," << actual_f.beta << ","
                << actual_tw.p << "," << actual_tw.alpha << "," << actual_tw.beta << ","
                << dp_diff << "," << da_diff << "," << db_diff << "," << max_diff << "\n";
    }
    outfile.close();

    std::cout << "=== PHI-BASED TRANSFORMED COMPARISON ===" << std::endl;
    std::cout << "Fereisl steps: " << result_fereisl.num_steps << std::endl;
    std::cout << "Tucker-Will steps: " << result_tw.num_steps << std::endl;
    std::cout << "Phi-based actual comparison written to: " << output_file << std::endl;
}

struct EpsScanResult {
    double eps;
    double phi;
    double max_diff;
    int fereisl_steps;
    int tw_steps;
};

EpsScanResult runEpsilonConvergenceScan(
    const BinaryState& initial_state,
    PhysicalParams params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance
) {
    AdaptiveGaussCollocationIntegrator integrator(tolerance);
    auto rhs_fereisl = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS(s, p, max_PN_order);
    };
    auto rhs_tw = [max_PN_order](const BinaryState& s, const PhysicalParams& p) {
        return compute_secular_RHS_TW(s, p, max_PN_order);
    };

    auto result_fereisl = integrator.integrate(initial_state, params, rhs_fereisl, theta_start, theta_end);
    auto result_tw = integrator.integrate(initial_state, params, rhs_tw, theta_start, theta_end);

    double max_diff = 0.0;
    size_t min_size = std::min(result_fereisl.states.size(), result_tw.states.size());
    for (size_t i = 0; i < min_size; ++i) {
        auto& state_f = result_fereisl.states[i];
        auto& state_tw_i = result_tw.states[i];
        double dp_diff = std::abs(state_f.p - state_tw_i.p);
        double da_diff = std::abs(state_f.alpha - state_tw_i.alpha);
        double db_diff = std::abs(state_f.beta - state_tw_i.beta);
        double local_max = std::max(dp_diff, std::max(da_diff, db_diff));
        max_diff = std::max(max_diff, local_max);
    }

    return {params.eps, max_diff, max_diff, result_fereisl.num_steps, result_tw.num_steps};
}

void writeEpsilonConvergenceCSV(
    const BinaryState& initial_state,
    const PhysicalParams& params,
    int max_PN_order,
    double theta_start,
    double theta_end,
    double tolerance,
    const std::string& output_file
) {
    std::ofstream outfile(output_file);
    outfile << std::scientific << std::setprecision(10);
    outfile << "eps,phi,max_diff,fereisl_steps,tw_steps\n";

    std::vector<double> eps_values = {0.032, 0.016, 0.008, 0.004, 0.002};
    for (double eps : eps_values) {
        PhysicalParams scan_params = params;
        scan_params.eps = eps;
        auto result = runEpsilonConvergenceScan(initial_state, scan_params, max_PN_order,
                                                theta_start, theta_end, tolerance);
        outfile << result.eps << "," << result.phi << "," << result.max_diff << ","
                << result.fereisl_steps << "," << result.tw_steps << "\n";
    }
    outfile.close();
    std::cout << "Epsilon convergence scan written to: " << output_file << std::endl;
}

// ============================================================================
// MAIN
// ============================================================================

int main() {
    BinaryState initial_state{10.0, 0.1, 0.1};
    PhysicalParams params{1.0, 1.0, 0.25, 0.0, 1.0};
    
    double theta_start = 0.0, theta_end = 100.0;
    double tolerance = 1e-12;
    int max_PN_order = 5;
    
    std::cout << "=== Adaptive RK4 Integration Comparison ===" << std::endl;
    std::cout << "PN order: " << max_PN_order << " (4.5PN)" << std::endl;
    std::cout << "Tolerance: " << tolerance << std::endl;
    std::cout << "Integration: theta = [" << theta_start << ", " << theta_end << "]" << std::endl;

    double phi_start = 0.0;
    double phi_end = params.eps * theta_end;
    std::cout << "Phi-transformed integration: phi = [" << phi_start << ", " << phi_end << "]" << std::endl;

    compareEvolutionMethodsPhi(initial_state, params, max_PN_order, phi_start, phi_end, tolerance,
                               "fereisl_tw_phi_comparison.csv");
    writeEpsilonConvergenceCSV(initial_state, params, max_PN_order, theta_start, theta_end, tolerance,
                               "fereisl_tw_eps_phi_convergence.csv");
    
    return 0;
}
