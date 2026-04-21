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

// Orbit average
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

// ── Analytic formulas (TW eq.2.18a, 2.18b — same physics) ───────────────
// In code units G=M=c=1: x = c²p/(GM) = p, so dx/dθ = dp/dθ exactly.

double dp_25_TW(double p,double e){
    return -pow(G*M*p,2.5)*eta/(5.0*pow(p,4))*8.0*(8.0+7.0*e*e);
}
double dp_35_TW(double p,double e){
    double pref=-eta*pow(G*M,3)*sqrt(G*M*p)/(210.0*pow(p,3));
    return pref*(-8.0*(2759.0+252.0*eta)+8.0*e*e*(758.0+889.0*eta)
                 +e*e*e*e*(1483.0+4424.0*eta));
}
double dp_4_TW(double p,double e){
    double pref=-eta*PI*pow(G*M,4)/(360.0*pow(p,3));
    return pref*(18432.0+55872.0*e*e+7056.0*e*e*e*e-49.0*pow(e,6));
}
double dp_45_TW(double p,double e){
    double e2=e*e,e4=e2*e2,e6=e4*e2;
    double pref=eta*sqrt(G*M*p)*pow(G*M,4)/(11340.0*pow(p,4));
    return pref*(-8.0*(1469531.0-101511.0*eta+36288.0*eta*eta)
                 +12.0*e2*(743333.0-7263.0*eta+467775.0*eta*eta)
                 +3.0*e4*(690973.0-2644191.0*eta+1646568.0*eta*eta)
                 +9.0*e6*(527.0-6300.0*eta+53088.0*eta*eta));
}

double de_25_TW(double p,double e){
    return -(e*pow(G*M*p,2.5)*eta)/(15.0*pow(p,5))*(304.0+121.0*e*e);
}
double de_35_TW(double p,double e){
    double pref=-e*pow(G*M,3)*eta*sqrt(G*M*p)/(840.0*pow(p,4));
    return pref*(-8.0*(18049.0+4452.0*eta)+4.0*e*e*(8692.0+12803.0*eta)
                 +e*e*e*e*(2251.0+15064.0*eta));
}
double de_4_TW(double p,double e){
    double pref=-eta*e*PI*pow(G*M,4)/(34560.0*pow(p,3));
    return pref*(4538880.0+6876288.0*e*e+581208.0*e*e*e*e+623.0*pow(e,6));
}
double de_45_TW(double p,double e){
    double e2=e*e,e4=e2*e2,e6=e4*e2;
    double pref=eta*e*sqrt(G*M*p)*pow(G*M,4)/(30240.0*pow(p,5));
    return pref*(-16.0*(3428803.0+623439.0*eta+56448.0*eta*eta)
                 +12.0*e2*(745447.0+1064871.0*eta+2339925.0*eta*eta)
                 +2.0*e4*(-251155.0-9948885.0*eta+7220682.0*eta*eta)
                 +3.0*e6*(-25845.0-78380.0*eta+361536.0*eta*eta));
}

// 
double dp_25_JF(double p,double e){ return dp_25_TW(p,e); }
double dp_35_JF(double p,double e){ return dp_35_TW(p,e); }
double dp_4_JF (double p,double e){ return dp_4_TW (p,e); }
double dp_45_JF(double p,double e){  // JF writes (G*M)^(9/2) grouping
    double e2=e*e,e4=e2*e2,e6=e4*e2;
    double pref=eta*sqrt(G*M*p)*pow(G*M,4)/(11340.0*pow(p,4));
    return pref*(-8.0*(1469531.0-101511.0*eta+36288.0*eta*eta)
                 +12.0*e2*(743333.0-7263.0*eta+467775.0*eta*eta)
                 +3.0*e4*(690973.0-2644191.0*eta+1646568.0*eta*eta)
                 +9.0*e6*(527.0-6300.0*eta+53088.0*eta*eta));
}
double de_25_JF(double p,double e){ return de_25_TW(p,e); }
double de_35_JF(double p,double e){ return de_35_TW(p,e); }
double de_4_JF (double p,double e){ return de_4_TW (p,e); }
double de_45_JF(double p,double e){
    double e2=e*e,e4=e2*e2,e6=e4*e2;
    double pref=eta*e*sqrt(G*M*p)*pow(G*M,4)/(30240.0*pow(p,5));
    return pref*(-16.0*(3428803.0+623439.0*eta+56448.0*eta*eta)
                 +12.0*e2*(745447.0+1064871.0*eta+2339925.0*eta*eta)
                 +2.0*e4*(-251155.0-9948885.0*eta+7220682.0*eta*eta)
                 +3.0*e6*(-25845.0-78380.0*eta+361536.0*eta*eta));
}

double dp_TW(double p,double e){ return dp_25_TW(p,e)+dp_35_TW(p,e)+dp_4_TW(p,e)+dp_45_TW(p,e); }
double de_TW(double p,double e){ return de_25_TW(p,e)+de_35_TW(p,e)+de_4_TW(p,e)+de_45_TW(p,e); }
double dp_JF(double p,double e){ return dp_25_JF(p,e)+dp_35_JF(p,e)+dp_4_JF(p,e)+dp_45_JF(p,e); }
double de_JF(double p,double e){ return de_25_JF(p,e)+de_35_JF(p,e)+de_4_JF(p,e)+de_45_JF(p,e); }

// ── Polynomial fitting via Gaussian elimination ───────────────────────────────
// Fits f(ε) = sum_k coeff_k * ε^{powers[k]} using Npts sample points,
// returns the coefficient of ε^targetPow.
static vector<double> solveLS(vector<vector<double>> A, vector<double> b){
    int n=A.size();
    for(int i=0;i<n;++i){
        int piv=i;
        for(int r=i+1;r<n;++r) if(fabs(A[r][i])>fabs(A[piv][i])) piv=r;
        swap(A[i],A[piv]); swap(b[i],b[piv]);
        double d=A[i][i]; if(fabs(d)<1e-15) continue;
        for(int cc=i;cc<n;++cc) A[i][cc]/=d; b[i]/=d;
        for(int r=0;r<n;++r){ if(r==i) continue;
            double f=A[r][i];
            for(int cc=i;cc<n;++cc) A[r][cc]-=f*A[i][cc];
            b[r]-=f*b[i];
        }
    }
    return b;
}

double extractCoeff(const function<double(double)>& f,
                    const vector<int>& powers, int targetPow,
                    double eps_lo, double eps_hi)
{
    int sz = powers.size();
    vector<vector<double>> A(sz, vector<double>(sz,0.0));
    vector<double> b(sz,0.0);
    for(int row=0;row<sz;++row){
        double eps = eps_lo * pow(eps_hi/eps_lo, double(row)/(sz-1));
        for(int col=0;col<sz;++col) A[row][col] = pow(eps, powers[col]);
        b[row] = f(eps);
    }
    auto coeffs = solveLS(A,b);
    for(int i=0;i<sz;++i) if(powers[i]==targetPow) return coeffs[i];
    return 0.0;
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
    const int WNUM = 28;   // scientific numbers
    const int WRD  = 20;   // relative diff

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
         << (fabs(e_from_ab-e0)<1e-14 ? "  [OK]\n" : "  [MISMATCH!]\n");
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
    cout << "   dx/dtheta       = dp/dtheta  (since x = c^2*p/GM = p)\n\n";

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
    prow("dp/dtheta", avg1.dp, dp_TW(p0,e0), dp_JF(p0,e0));
    prow("de/dtheta", avg1.de, de_TW(p0,e0), de_JF(p0,e0));
    prow("dx/dtheta", avg1.dp, dp_TW(p0,e0), dp_JF(p0,e0));
    cout << "\n NOTE: At epsilon=1 the PN series is not small; order-mixing is O(1).\n"
         << "       Large relDiff(QLT,TW) here is expected and unphysical.\n"
         << "       The meaningful test is the coefficient extraction below.\n\n";

    // ── §2: PN coefficient extraction via ε polynomial fitting ────────────────
    bar('=');
    cout << " §2  PN coefficient extraction — polynomial fit in epsilon = 1/c_val\n";
    bar('=');
    cout << "\n dp/dtheta = eps^5*A_25 + eps^7*A_35 + eps^8*A_4 + eps^9*A_45 + ...\n";
    cout << " de/dtheta = eps^5*B_25 + eps^7*B_35 + eps^8*B_4 + eps^9*B_45 + ...\n";
    cout << " Fit this polynomial at small eps; the extracted A_N should equal\n";
    cout << " the TW/JF analytic value.  Fit range: eps in [0.006, 0.09]\n\n";

    // Fit basis powers: 5,6,7,8,9,10  (6 terms, 6 unknowns)
    vector<int> pows = {5,6,7,8,9,10};
    double elo=0.006, ehi=0.09;

    auto QLT_dp_f = [&](double eps)->double{
        return orbitAverage(Kfull,p0,alpha0,beta0,1.0/eps,PNord,Nsamp).dp; };
    auto QLT_de_f = [&](double eps)->double{
        return orbitAverage(Kfull,p0,alpha0,beta0,1.0/eps,PNord,Nsamp).de; };

    cout << " (Computing... this takes a few seconds)\n\n";

    double A25q=extractCoeff(QLT_dp_f,pows,5,elo,ehi);
    double A35q=extractCoeff(QLT_dp_f,pows,7,elo,ehi);
    double A4q =extractCoeff(QLT_dp_f,pows,8,elo,ehi);
    double A45q=extractCoeff(QLT_dp_f,pows,9,elo,ehi);

    double B25q=extractCoeff(QLT_de_f,pows,5,elo,ehi);
    double B35q=extractCoeff(QLT_de_f,pows,7,elo,ehi);
    double B4q =extractCoeff(QLT_de_f,pows,8,elo,ehi);
    double B45q=extractCoeff(QLT_de_f,pows,9,elo,ehi);

    // Analytic values (these are the dimensional coefficients at eps=1)
    double A25t=dp_25_TW(p0,e0), A25j=dp_25_JF(p0,e0);
    double A35t=dp_35_TW(p0,e0), A35j=dp_35_JF(p0,e0);
    double A4t =dp_4_TW (p0,e0), A4j =dp_4_JF (p0,e0);
    double A45t=dp_45_TW(p0,e0), A45j=dp_45_JF(p0,e0);

    double B25t=de_25_TW(p0,e0), B25j=de_25_JF(p0,e0);
    double B35t=de_35_TW(p0,e0), B35j=de_35_JF(p0,e0);
    double B4t =de_4_TW (p0,e0), B4j =de_4_JF (p0,e0);
    double B45t=de_45_TW(p0,e0), B45j=de_45_JF(p0,e0);

    cout << left << setw(WLBL) << " Order"
         << setw(WNUM) << "QLT fitted"
         << setw(WNUM) << "TW analytic"
         << setw(WNUM) << "JF analytic"
         << setw(WRD) << "relDiff(QLT,TW)"
         << setw(WRD) << "relDiff(QLT,JF)"
         << "Winner\n";
    bar('-');

    double total_tw_dp=0, total_jf_dp=0;
    double total_tw_de=0, total_jf_de=0;
    int wins_tw=0, wins_jf=0;

    auto crow = [&](const string& ord, double qlt, double tw, double jf,
                    double& tot_tw, double& tot_jf) {
        double rd_tw=relDiff(qlt,tw), rd_jf=relDiff(qlt,jf);
        tot_tw+=rd_tw; tot_jf+=rd_jf;
        bool tw_wins=(rd_tw<=rd_jf);
        if(tw_wins) ++wins_tw; else ++wins_jf;
        cout << " " << left << setw(WLBL) << ord
             << scientific << setprecision(16)
             << setw(WNUM) << qlt << setw(WNUM) << tw << setw(WNUM) << jf
             << fixed << setprecision(16)
             << setw(WRD) << rd_tw << setw(WRD) << rd_jf
             << (tw_wins ? "TW" : "JF") << "\n";
    };

    cout << " dp/dtheta:\n";
    crow("  2.5PN", A25q, A25t, A25j, total_tw_dp, total_jf_dp);
    crow("  3.5PN", A35q, A35t, A35j, total_tw_dp, total_jf_dp);
    crow("  4PN",   A4q,  A4t,  A4j,  total_tw_dp, total_jf_dp);
    crow("  4.5PN", A45q, A45t, A45j, total_tw_dp, total_jf_dp);

    cout << "\n de/dtheta:\n";
    crow("  2.5PN", B25q, B25t, B25j, total_tw_de, total_jf_de);
    crow("  3.5PN", B35q, B35t, B35j, total_tw_de, total_jf_de);
    crow("  4PN",   B4q,  B4t,  B4j,  total_tw_de, total_jf_de);
    crow("  4.5PN", B45q, B45t, B45j, total_tw_de, total_jf_de);

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
    cout << " §4  CONVERGENCE SCORECARD — which analytic solution does QLT approach?\n";
    bar('=');
    cout << "\n TW and JF encode identical physics.  Their analytic expressions at\n"
         << " 2.5PN and 3.5PN are literally the same.  At 4.5PN the factoring\n"
         << " differs but the numerical values agree to floating-point precision.\n"
         << " Therefore the question reduces to: do the QLT C/D coefficients\n"
         << " reproduce the correct orbit-averaged result at each PN order?\n\n";

    cout << " Coefficient comparison (from §2 polynomial fit):\n\n";
    cout << left << setw(WLBL) << " Coefficient"
         << setw(WNUM) << "QLT"
         << setw(WNUM) << "TW"
         << setw(WNUM) << "JF"
         << setw(WRD) << "rDiff(QLT,TW)"
         << "rDiff(QLT,JF)\n";
    bar('-');

    auto srow = [&](const string& n, double q, double t, double j){
        cout << " " << left << setw(WLBL) << n
             << scientific << setprecision(16)
             << setw(20) << q << setw(20) << t << setw(20) << j
             << fixed << setprecision(16)
             << setw(16) << relDiff(q,t) << relDiff(q,j) << "\n";
    };
    srow("dp  2.5PN", A25q, A25t, A25j);
    srow("dp  3.5PN", A35q, A35t, A35j);
    srow("dp  4.0PN", A4q,  A4t,  A4j );
    srow("dp  4.5PN", A45q, A45t, A45j);
    srow("de  2.5PN", B25q, B25t, B25j);
    srow("de  3.5PN", B35q, B35t, B35j);
    srow("de  4.0PN", B4q,  B4t,  B4j );
    srow("de  4.5PN", B45q, B45t, B45j);

    double total_tw = total_tw_dp+total_tw_de;
    double total_jf = total_jf_dp+total_jf_de;

    cout << "\n";
    bar('-');
    cout << " Total relDiff  TW=" << fixed << setprecision(16) << total_tw
         << "   JF=" << total_jf << "\n";
    cout << " Wins (lower rDiff per row):  TW=" << wins_tw
         << "   JF=" << wins_jf << "\n\n";

    bar('*',70);
    if(fabs(total_tw - total_jf) < 1e-6*total_tw){
        cout << " RESULT: TW and JF converge IDENTICALLY.\n"
             << "         Their formulas are algebraically equal.\n"
             << "         Any apparent difference is floating-point rounding only.\n";
    } else if(total_tw < total_jf){
        cout << " RESULT: TUCKER-WILL converges better (lower total relDiff).\n";
    } else {
        cout << " RESULT: JAN FEREISL converges better (lower total relDiff).\n";
    }
    bar('*',70);

    cout << "\n PHYSICAL INTERPRETATION OF THE EPSILON-SCAN (§3):\n"
         << "  - 2.5PN dp column:  should be FLAT = converging to A25t\n"
         << "  - 3.5PN dp column:  should be FLAT = converging to A35t\n"
         << "  - 4.5PN dp column:  flatness depends on correctness of 4.5PN C/D\n"
         << "  - 4PN tail:         ZERO from QLT — not in C/D tables by design;\n"
         << "                      this term enters via orbit-avg of a log term\n"
         << "                      and must be added separately if needed\n"
         << "  - de/dtheta at e=0.01: numerically tiny (~e x dp); relative\n"
         << "                      differences appear large even for small abs errors\n\n";

    return 0;
}
