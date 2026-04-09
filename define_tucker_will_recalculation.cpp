// This code compares the 4.5PN coefficients from Tucker & Will (2021) with the ones from the Mathematica notebook made by Jan Fereisl

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>
#include <functional>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;


// Constants (geometrized units)
constexpr double G = 1.0;               // Gravitational constant in geometrized units
constexpr double c = 1.0;               // Speed of light in geometrized units    

constexpr double M_sun = 1.98847e30;                // 1 Solar mass in kg
constexpr double G_SI = 6.6743e-11;                 // Gravitational constant in m^3/kg/s^2
constexpr double c_SI = 299792458;       

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



const int PNorderSetting = 9;
const int PNorder = PNorderSetting - 3;


// M_PI is not constexpr in <cmath>, so define our own
constexpr double PI  = 3.14159265358979323846;
constexpr double PI2 = PI * PI;

struct PNCoeffs {
    using Table = std::array<std::array<std::array<double,4>,4>,4>;

    Table a{}, b{}, c{}, d{};

    constexpr double A(int l, int m, int n) const { return a[l][m][n]; }
    constexpr double B(int l, int m, int n) const { return b[l][m][n]; }
    constexpr double C(int l, int m, int n) const { return c[l][m][n]; }
    constexpr double D(int l, int m, int n) const { return d[l][m][n]; }

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

// We omit including the 4.5PN gauge coefficients as they dissappear orbit-averaged equations for the orbit elements (as stated in the paper and simultaneously confirmed by Jan Fereisl)

// Defining the PN terms from the previous coefficients (eq. 2.2 in paper for direct use in the equations of motion)
// Conservative term (1PN, 2PN, 3PN)
double A_c(int N, const PNCoeffs& coeffs,const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val) {
    double sum = 0.0;
    for (int l = 0; l <= N; ++l) {
        for (int m = 0; m <= N - l; ++m) {
            int n = N - l - m;
            if (n >= 0) {
                sum += coeffs.A(l, m, n) * std::pow(r_dot_sq, m) * std::pow(v_dot_v, n) * std::pow(GM_over_r, l) / std::pow(c_val, 2 * N);
            }
        }
    }
    return sum;
}

// Conservative term (1PN, 2PN, 3PN)
double B_c(int N, const PNCoeffs& coeffs,const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val) {
    double sum = 0.0;
    for (int l = 0; l <= N; ++l) {
        for (int m = 0; m <= N - l; ++m) {
            int n = N - l - m;
            if (n >= 0) {
                sum += coeffs.B(l, m, n) * std::pow(r_dot_sq, m) * std::pow(v_dot_v, n) * std::pow(GM_over_r, l) / std::pow(c_val, 2 * N);
            }
        }
    }
    return sum;
}

// Radiation reaction term (2.5PN, 3.5PN, 4.5PN)
double A_rr(int N, const PNCoeffs& coeffs,const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val) {
    double sum = 0.0;
    for (int l = 0; l <= N; ++l) {
        for (int m = 0; m <= N - l; ++m) {
            int n = N - l - m;
            if (n >= 0) {
                sum += coeffs.C(l, m, n) * std::pow(r_dot_sq, m) * std::pow(v_dot_v, n) * std::pow(GM_over_r, l) / std::pow(c_val, 2 * N);
            }
        }
    }
    return sum;
}

// Radiation reaction term (2.5PN, 3.5PN, 4.5PN)
double B_rr(int N, const PNCoeffs& coeffs,const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val) {
    double sum = 0.0;
    for (int l = 0; l <= N; ++l) {
        for (int m = 0; m <= N - l; ++m) {
            int n = N - l - m;
            if (n >= 0) {
                sum += coeffs.D(l, m, n) * std::pow(r_dot_sq, m) * std::pow(v_dot_v, n) * std::pow(GM_over_r, l) / std::pow(c_val, 2 * N);
            }
        }
    }
    return sum;
}

// The conservative term at 3PN is dropped out as it is not used in the calculation.
// For it to show effect in the secular evolution, we would need to calculate up to 3 + 2.5 = 5.5 PN order

// Sum of conservative terms (1PN and 2PN) - the 3PN term is not included as it is not used in the calculation
double A_tot(const PNCoeffs& coeffs,const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val, int PNorder) {
    double sum = 0.0;
    for (int N = 1; N <= PNorder; ++N) {
        sum += A_c(N, coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val);
    }
    return sum;
}

// Sum of conservative terms (1PN and 2PN) - the 3PN term is not included as it is not used in the calculation
double B_tot(const PNCoeffs& coeffs, const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val, int PNorder) {
    double sum = 0.0;
    for (int N = 1; N <= PNorder; ++N) {
        sum += B_c(N, coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val);
    }
    return sum;
}

// Sum of radiation reaction terms (2.5PN, 3.5PN, 4.5PN)
double A_rrtot(const PNCoeffs& coeffs, const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val, int PNorder) {
    double sum = 0.0;
    for (int N = 1; N <= PNorder; ++N) {
        sum += A_rr(N, coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val);
    }
    return sum;
}

// Sum of radiation reaction terms (2.5PN, 3.5PN, 4.5PN)
double B_rrtot(const PNCoeffs& coeffs, const double& r_dot_sq, const double& v_dot_v, const double& GM_over_r, const double& c_val, int PNorder) {
    double sum = 0.0;
    for (int N = 1; N <= PNorder; ++N) {
        sum += B_rr(N, coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val);
    }
    return sum;
}

// Orbital elements definitions and transformations

// NormR = p/(1 + e*Cos[f])
double NormR(const double& p,const double& e,const double& alpha, const double& beta, const double& phi) {
    // 1 + e*Cos[f] = 1 + (α Cos[φ] + β Sin[φ])
    double cos_f = alpha * cos(phi) + beta * sin(phi);
    return p / (1.0 + e * cos_f);
}

// rDot = (Sqrt[G*M*p]*e/p)*Sin[f]
double rDot(const double& p, const double& e, const double& alpha, const double& beta, const double& phi) {
    // Sin[f] = α Sin[φ] - β Cos[φ]
    double sin_f = alpha * sin(phi) - beta * cos(phi);
    return (sqrt(G * M * p) * e / p) * sin_f;
}

// NormV = Sqrt[rDot^2 + (Sqrt[G*M*p]/NormR)^2]
double NormV(const double& p, const double& e, const double& alpha, const double& beta, const double& phi) {
    double r_dot = rDot(p, e, alpha, beta, phi);
    double norm_r = NormR(p, e, alpha, beta, phi);
    return sqrt(r_dot * r_dot + (sqrt(G * M * p) / norm_r) * (sqrt(G * M * p) / norm_r));
}



// Conservative terms (*Transformation of conservative A and B coefficients*)

double AplusB(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double norm_v = NormV(p, e, alpha, beta, phi);
    double r_dot_sq = r_dot * r_dot;
    double v_dot_v = norm_v * norm_v;
    double GM_over_r = G * M / norm_r;
    double a_tot = A_tot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
    double b_tot = B_tot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
    return a_tot + b_tot;
}


double Btot(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double norm_v = NormV(p, e, alpha, beta, phi);
    double r_dot_sq = r_dot * r_dot;
    double v_dot_v = norm_v * norm_v;
    double GM_over_r = G * M / norm_r;
    return B_tot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
}


double ScriptCapitalRcons(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double aplusb = AplusB(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    return G * M / (norm_r * norm_r) * aplusb;
}


double ScriptCapitalScons(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double btot = Btot(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    return G * M / (norm_r * norm_r * norm_r) * (sqrt(G * M * p) / r_dot) * btot;
}

// Radiation Reaction terms (*Transformation of radiation reaction A and B coefficients*)


double AplusBrr(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double norm_v = NormV(p, e, alpha, beta, phi);
    double r_dot_sq = r_dot * r_dot;
    double v_dot_v = norm_v * norm_v;
    double GM_over_r = G * M / norm_r;
    double a_rrtot = A_rrtot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
    double b_rrtot = B_rrtot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
    return a_rrtot + b_rrtot;
}

double Btotrr(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double norm_v = NormV(p, e, alpha, beta, phi);
    double r_dot_sq = r_dot * r_dot;
    double v_dot_v = norm_v * norm_v;
    double GM_over_r = G * M / norm_r;
    return B_rrtot(coeffs, r_dot_sq, v_dot_v, GM_over_r, c_val, PNorder);
}


double ScriptCapitalRrr(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double r_dot = rDot(p, e, alpha, beta, phi);
    double aplusbrr = AplusBrr(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    double epsilon = 1.0 / c_val;  // Since c = 1/ε in geometrized units, but ε is small, but in code c=1, so ε=1
    return (8.0 / 5.0) * eta * std::pow(epsilon, 3) * std::pow(G * M, 2) / std::pow(norm_r, 3) * r_dot * aplusbrr;
}


double ScriptCapitalSrr(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double btotrr = Btotrr(coeffs, p, e, alpha, beta, phi, c_val, PNorder + 3);  // Sum up to PNorder + 3
    double epsilon = 1.0 / c_val;
    double GM_over_normr_sq = G * M / (norm_r * norm_r);
    return (8.0 / 5.0) * eta * std::pow(epsilon, 3) * std::pow(GM_over_normr_sq, 2) * sqrt(G * M * p) * btotrr;
}

// Full equations of motion right-hand sides combining conservative and radiation reaction parts

double ScriptCapitalR(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    return ScriptCapitalRrr(coeffs, p, e, alpha, beta, phi, c_val, PNorder)
         + ScriptCapitalRcons(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
}

double ScriptCapitalS(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    return ScriptCapitalSrr(coeffs, p, e, alpha, beta, phi, c_val, PNorder)
         + ScriptCapitalScons(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
}

// dp/dphi = 2 NormR^3 / (G M) * ScriptCapitalS (eq.2.8 in the paper)
double QLT_p(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double scriptS = ScriptCapitalS(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    return 2.0 * norm_r * norm_r * norm_r / (G * M) * scriptS;
}

// dalpha/dphi = NormR^2 / (G M) * [ScriptCapitalR Sin[φ] + ScriptCapitalS (α + Cos[φ]) (1 + NormR/p)] (eq.2.8 in the paper)
double QLT_alpha(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double scriptR = ScriptCapitalR(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    double scriptS = ScriptCapitalS(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    return norm_r * norm_r / (G * M) * (
        scriptR * sin(phi)
      + scriptS * (alpha + cos(phi)) * (1.0 + norm_r / p)
      - scriptS * alpha
    );
}

// dbeta/dphi = NormR^2 / (G M) * [-ScriptCapitalR Cos[φ] + ScriptCapitalS (β + Sin[φ]) (1 + NormR/p)] (eq.2.8 in the paper)
double QLT_beta(const PNCoeffs& coeffs, const double& p, const double& e, const double& alpha, const double& beta, const double& phi, const double& c_val, int PNorder) {
    double norm_r = NormR(p, e, alpha, beta, phi);
    double scriptR = ScriptCapitalR(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    double scriptS = ScriptCapitalS(coeffs, p, e, alpha, beta, phi, c_val, PNorder);
    return norm_r * norm_r / (G * M) * (
        -scriptR * cos(phi)
      + scriptS * (beta + sin(phi)) * (1.0 + norm_r / p)
      - scriptS * beta
    );
}

// Orbit-averaged evolution helpers

using State3 = std::array<double, 3>;
static const int Ord = static_cast<int>(PNorderSetting);
static const int NonZeroOrd = 2;
static const std::array<int, 5> QOrders = {2, 4, 5, 7, 9};

// A simple numerical integration method to compute the average of a function over φ from 0 to 2π
template<typename Func>
double Avg(Func f, int samples = 2048) {
    double sum = 0.0;
    for (int i = 0; i < samples; ++i) {
        double phi = 2.0 * PI * (i + 0.5) / samples;
        sum += f(phi);
    }
    return sum / static_cast<double>(samples);
}

// We will use the same numerical integration method to compute the average of φ times a function over φ from 0 to 2π
using YpFunction = std::function<double(double, double, double, double)>;

struct YpSeries {
    std::array<std::vector<YpFunction>, 3> terms;
};

// A helper function to compute ε^order, which is used in the expansion of YpSeries terms. This allows us to easily adjust the order of the expansion by changing the value of ε.
inline double epsilonPow(int order, double epsilon) {
    return std::pow(epsilon, order);
}

// A helper function to compute the contribution of the YpSeries terms for a given index (idx) and φ. It sums up the contributions from the non-zero order terms, weighted by ε^order.
double Yexp(const YpSeries& yps, int idx, double pt, double alpha_t, double beta_t, double phi, double epsilon) {
    double sum = 0.0;
    for (int i = NonZeroOrd; i <= Ord - NonZeroOrd; ++i) {
        int termIndex = i - NonZeroOrd;
        if (termIndex >= 0 && termIndex < static_cast<int>(yps.terms[idx].size())) {
            sum += epsilonPow(i, epsilon) * yps.terms[idx][termIndex](pt, alpha_t, beta_t, phi);
        }
    }
    return sum;
}

// A helper function to compute the new state XindPhi by adding the contributions from the YpSeries terms to the current state Xvar. This function is used to compute the intermediate state at a given φ, which is then used in the Gseries function to compute the final Q expansion.
State3 XindPhi(const YpSeries& yps, const State3& Xvar, double pt, double alpha_t, double beta_t, double phi, double epsilon) {
    return {
        Xvar[0] + Yexp(yps, 0, pt, alpha_t, beta_t, phi, epsilon),
        Xvar[1] + Yexp(yps, 1, pt, alpha_t, beta_t, phi, epsilon),
        Xvar[2] + Yexp(yps, 2, pt, alpha_t, beta_t, phi, epsilon)
    };
}

// A helper function to create a PhiFunc that computes the contribution of a specific QSeries term for a given index (idx) and order. This function captures the necessary parameters and returns a function that can be evaluated at any φ to get the corresponding Q coefficient value.
using QFunc = std::function<double(const PNCoeffs&, double, double, double, double, int)>;
using PhiFunc = std::function<double(double)>;

struct QSeries {
    std::array<std::array<QFunc, 3>, 5> terms;
};

// A helper function to create a PhiFunc that computes the contribution of a specific QSeries term for a given index (idx) and order. This function captures the necessary parameters and returns a function that can be evaluated at any φ to get the corresponding Q coefficient value.
PhiFunc makeQCoefPhi(const QSeries& qseries, int idx, int order, const State3& Xvar, double pt, double alpha_t, double beta_t, double epsilon, const PNCoeffs& coeffs, int PNorder) {
    return [=](double phi) {
        double value = 0.0;
        for (int k = 0; k < static_cast<int>(QOrders.size()); ++k) {
            if (QOrders[k] == order) {
                value = qseries.terms[k][idx](coeffs, Xvar[0], Xvar[1], Xvar[2], phi, PNorder);
                break;
            }
        }
        return value;
    };
}

// A helper function to numerically integrate a PhiFunc from 0 to upperPhi using the midpoint rule. This function is used to compute the integral of the Q coefficients over φ, which is necessary for computing the YpSeries terms and the final averaged equations of motion.
double integratePhi(const PhiFunc& f, double upperPhi, int samples = 2048) {
    if (upperPhi == 0.0) {
        return 0.0;
    }
    double sum = 0.0;
    double h = upperPhi / samples;
    for (int i = 0; i < samples; ++i) {
        double x = (i + 0.5) * h;
        sum += f(x);
    }
    return sum * h;
}

// A helper function to compute the average of a PhiFunc over φ from 0 to 2π. This function uses the integratePhi function to compute the integral and then divides by 2π to get the average value.
double avgPhi(const PhiFunc& f, int samples = 2048) {
    return integratePhi(f, 2.0 * PI, samples) / (2.0 * PI);
}

// A helper function to compute the average of φ times a PhiFunc over φ from 0 to 2π. This function uses a numerical integration method similar to integratePhi, but it multiplies the function value by φ at each step before summing, and then divides by 2π to get the average value.
double avgPhiTimes(const PhiFunc& f, int samples = 2048) {
    double sum = 0.0;
    double h = 2.0 * PI / samples;
    for (int i = 0; i < samples; ++i) {
        double x = (i + 0.5) * h;
        sum += x * f(x);
    }
    return sum * h / (2.0 * PI);
}

// A helper function to reduce the order of a QSeries by zeroing out the terms that exceed the specified maxOrder. This function creates a new QSeries where the terms corresponding to orders greater than maxOrder are replaced with functions that return 0, effectively removing their contribution from the calculations.
QSeries ReduceOrder(const QSeries& qseries, int maxOrder) {
    QSeries reduced;
    for (int k = 0; k < static_cast<int>(QOrders.size()); ++k) {
        if (QOrders[k] <= maxOrder) {
            reduced.terms[k] = qseries.terms[k];
        } else {
            for (int idx = 0; idx < 3; ++idx) {
                reduced.terms[k][idx] = [](const PNCoeffs&, double, double, double, double, int) {
                    return 0.0;
                };
            }
        }
    }
    return reduced;
}

// A helper function to create a PhiFunc for a specific QSeries term based on the given parameters. This function captures the necessary parameters and returns a function that can be evaluated at any φ to get the corresponding Q coefficient value for the specified index and order.
PhiFunc GetCoeff(const QSeries& qseries, int idx, int order, const State3& Xvar, double pt, double alpha_t, double beta_t, double epsilon, const PNCoeffs& coeffs, int PNorder) {
    return makeQCoefPhi(qseries, idx, order, Xvar, pt, alpha_t, beta_t, epsilon, coeffs, PNorder);
}

// A helper function to create a PhiFunc that computes the contribution of a specific QSeries term for a given index (idx) and order. This function captures the necessary parameters and returns a function that can be evaluated at any φ to get the corresponding Y solution value, which is used in the YpSeries terms.
PhiFunc makeYSolPhi(const PhiFunc& qcoef) {
    double avg_q = avgPhi(qcoef);
    double avg_phi_q = avgPhiTimes(qcoef);
    return [=](double phi) {
        return integratePhi(qcoef, phi) - (phi + PI) * avg_q + avg_phi_q;
    };
}

// A helper function to compute the G series expansion for a given QSeries, PNCoeffs, intermediate state xphi, and other parameters. This function sums up the contributions from the QSeries terms weighted by ε^order to compute the final G series value, which is used in the Q expansion of the equations of motion.
State3 Gseries(const QSeries& qseries, const PNCoeffs& coeffs, const State3& xphi, double phi, double epsilon, int PNorder) {
    State3 result = {0.0, 0.0, 0.0};
    for (int k = 0; k < static_cast<int>(QOrders.size()); ++k) {
        double orderFactor = epsilonPow(QOrders[k], epsilon);
        for (int idx = 0; idx < 3; ++idx) {
            result[idx] += orderFactor * qseries.terms[k][idx](coeffs, xphi[0], xphi[1], xphi[2], phi, PNorder);
        }
    }
    return result;
}

// A helper function to compute the Q expansion for a given QSeries, YpSeries, current state Xvar, and other parameters. This function first computes the intermediate state xphi using the XindPhi function, and then uses the Gseries function to compute the final Q expansion value, which is used in the equations of motion.
State3 Qexpansion(const QSeries& qseries, const YpSeries& yps, const State3& Xvar, double pt, double alpha_t, double beta_t, double phi, double epsilon, const PNCoeffs& coeffs, int PNorder) {
    State3 xphi = XindPhi(yps, Xvar, pt, alpha_t, beta_t, phi, epsilon);
    return Gseries(qseries, coeffs, xphi, phi, epsilon, PNorder);
}

// A helper function to build the YpSeries based on the given QSeries, current state Xvar, and other parameters. This function iterates over the indices and orders of the YpSeries terms, creates the corresponding Q coefficient PhiFunc using makeQCoefPhi, and then creates the Y solution PhiFunc using makeYSolPhi. The resulting YpSeries is returned for use in the calculations of the equations of motion.
YpSeries buildYpSeries(const QSeries& qseries, const State3& Xvar, double pt, double alpha_t, double beta_t, double epsilon, const PNCoeffs& coeffs, int PNorder) {
    YpSeries yps;
    for (int idx = 0; idx < 3; ++idx) {
        for (int order = NonZeroOrd; order <= Ord - NonZeroOrd; ++order) {
            PhiFunc qcoef = makeQCoefPhi(qseries, idx, order, Xvar, pt, alpha_t, beta_t, epsilon, coeffs, PNorder);
            PhiFunc ysol = makeYSolPhi(qcoef);
            yps.terms[idx].push_back([=](double pt_in, double at, double bt, double phi_in) {
                return ysol(phi_in);
            });
        }
    }
    return yps;
}

// A helper function to compute the average of the final right-hand side of the equations of motion over φ from 0 to 2π. This function defines a lambda function rhsFunc that computes the Q expansion for a given φ, and then uses numerical integration to compute the average value of this function over φ, which gives the averaged equations of motion for the orbital elements.
State3 avgFinalRHS(const QSeries& qseries, const YpSeries& yps, const State3& Xvar, double pt, double alpha_t, double beta_t, double epsilon, const PNCoeffs& coeffs, int PNorder) {
    auto rhsFunc = [&](double phi) {
        State3 qval = Qexpansion(qseries, yps, Xvar, pt, alpha_t, beta_t, phi, epsilon, coeffs, PNorder);
        return qval;
    };
    State3 avgRHS = {0.0, 0.0, 0.0};
    int samples = 2048;
    double h = 2.0 * PI / samples;
    for (int i = 0; i < samples; ++i) {
        double phi = (i + 0.5) * h;
        State3 qval = rhsFunc(phi);
        avgRHS[0] += qval[0];
        avgRHS[1] += qval[1];
        avgRHS[2] += qval[2];
    }
    avgRHS[0] *= h / (2.0 * PI);
    avgRHS[1] *= h / (2.0 * PI);
    avgRHS[2] *= h / (2.0 * PI);
    return avgRHS;
}

// A helper function to generate all combinations of non-negative integers (i, j, k) such that i + j + k <= Max. This function is used to generate the indices for the QSeries and YpSeries terms based on the specified maximum order of the expansion.
static const int MaxDer = 2;

// Generate all combinations of non-negative integers (i, j, k) such that i + j + k <= Max. This is used to generate the indices for the QSeries and YpSeries terms based on the specified maximum order of the expansion.
std::vector<std::array<int, 3>> CombinationsIJK(int Max) {
    std::vector<std::array<int, 3>> results;
    results.reserve(Max * (Max + 1) * (Max + 2) / 6);
    for (int A = 1; A <= Max; ++A) {
        for (int i = 0; i <= A; ++i) {
            for (int j = 0; j <= A - i; ++j) {
                int k = A - i - j;
                if (k >= 0) {
                    results.push_back({i, j, k});
                }
            }
        }
    }
    return results;
}

// Solve the linear system Ax = b using Gaussian elimination with partial pivoting and return the solution vector x
std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = static_cast<int>(A.size());
    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int row = i + 1; row < n; ++row) {
            if (fabs(A[row][i]) > fabs(A[pivot][i])) {
                pivot = row;
            }
        }
        if (pivot != i) {
            swap(A[i], A[pivot]);
            swap(b[i], b[pivot]);
        }
        double diag = A[i][i];
        if (fabs(diag) < 1e-15) {
            continue;
        }
        for (int col = i; col < n; ++col) {
            A[i][col] /= diag;
        }
        b[i] /= diag;
        for (int row = 0; row < n; ++row) {
            if (row == i) {
                continue;
            }
            double factor = A[row][i];
            for (int col = i; col < n; ++col) {
                A[row][col] -= factor * A[i][col];
            }
            b[row] -= factor * b[i];
        }
    }
    return b;
}

// Fit a polynomial of degree maxOrder to the function f evaluated at small epsilon values and return the coefficients of the fitted polynomial
std::vector<double> getPolynomialCoefficients(const std::function<double(double)>& f, int maxOrder) {
    int size = maxOrder + 1;
    const double h = 1e-3;
    std::vector<std::vector<double>> A(size, std::vector<double>(size));
    std::vector<double> b(size);
    for (int row = 0; row < size; ++row) {
        double x = (row + 1) * h;
        double xp = 1.0;
        for (int col = 0; col < size; ++col) {
            A[row][col] = xp;
            xp *= x;
        }
        b[row] = f(x);
    }
    return solveLinearSystem(std::move(A), std::move(b));
}

// Evaluate the QLT expansion at a given epsilon by first calculating the YpSeries and then using it to compute the QSeries expansion, which gives the values of dp/dphi, dalpha/dphi, dbeta/dphi at that epsilon
State3 evaluateQLT(const PNCoeffs& coeffs, const State3& xvar, double phi, double epsilon, int PNorder) {
    double c_val = 1.0 / epsilon;
    double alpha_t = xvar[1];
    double beta_t = xvar[2];
    double e = sqrt(alpha_t * alpha_t + beta_t * beta_t);
    return {
        QLT_p(coeffs, xvar[0], e, alpha_t, beta_t, phi, c_val, PNorder),
        QLT_alpha(coeffs, xvar[0], e, alpha_t, beta_t, phi, c_val, PNorder),
        QLT_beta(coeffs, xvar[0], e, alpha_t, beta_t, phi, c_val, PNorder)
    };
}

// Get the coefficient of the QLT expansion for a specific qIndex and order by fitting a polynomial to the QLT values evaluated at small epsilon
double GetCoeffQLT(const PNCoeffs& coeffs, int qIndex, int order, const State3& xvar, double phi, int PNorder) {
    auto f = [&](double eps) {
        State3 q = evaluateQLT(coeffs, xvar, phi, eps, PNorder);
        return q[qIndex];
    };
    std::vector<double> coeffsVec = getPolynomialCoefficients(f, Ord);
    if (order < 0 || order >= static_cast<int>(coeffsVec.size())) {
        return 0.0;
    }
    return coeffsVec[order];
}

// Partial derivative of g with respect to p, alpha, beta according to the orders specified by dp, da, db
static double partialDerivative(const std::function<double(const State3&)>& g, const State3& x, int dp, int da, int db) {
    if (dp == 0 && da == 0 && db == 0) {
        return g(x);
    }
    const double h = 1e-6;
    if (dp > 0) {
        State3 xp = x;
        State3 xm = x;
        xp[0] += h;
        xm[0] -= h;
        return (partialDerivative(g, xp, dp - 1, da, db) - partialDerivative(g, xm, dp - 1, da, db)) / (2.0 * h);
    }
    if (da > 0) {
        State3 xp = x;
        State3 xm = x;
        xp[1] += h;
        xm[1] -= h;
        return (partialDerivative(g, xp, dp, da - 1, db) - partialDerivative(g, xm, dp, da - 1, db)) / (2.0 * h);
    }
    State3 xp = x;
    State3 xm = x;
    xp[2] += h;
    xm[2] -= h;
    return (partialDerivative(g, xp, dp, da, db - 1) - partialDerivative(g, xm, dp, da, db - 1)) / (2.0 * h);
}

double GetQDerivative(const PNCoeffs& coeffs, int qIndex, int qOrder, const std::array<int, 3>& derivOrder, const State3& xvar, double phi, int PNorder) {
    auto g = [&](const State3& v) {
        return GetCoeffQLT(coeffs, qIndex, qOrder, v, phi, PNorder);
    };
    return partialDerivative(g, xvar, derivOrder[0], derivOrder[1], derivOrder[2]);
}

std::array<double, 2> FinalFinalRHS(const State3& preFinalRHS, double alpha_t, double beta_t, double e) {
    return {
        preFinalRHS[0],
        (alpha_t / e) * preFinalRHS[1] + (beta_t / e) * preFinalRHS[2]
    };
}



// Tucker-Will results for evolution

// definition of the rescaling parameter x in Tucker-Will paper (eq. 2.17 in the paper)
double x_TW(const double& p) {
    return c * c * p / (G * M);
}

// Tucker-Will result of de/dtheta (eq. 2.18a in the paper)
// term1...2.5PN, term2...3.5PN, term3...4PN, term4...4.5PN
double de_TW_dtheta(const double& e, const double& x_TW) {
    double term1 = -((304.0 + 121.0 * e * e) / 15.0) * eta * e * std::pow(x_TW, -5.0 / 2.0);
    double factor1 = (1.0 / 30.0) * eta * e * std::pow(x_TW, - 7.0 / 2.0);
    double term2_1 = factor1 * ((1.0/28.0) * (144392.0 - 34768.0 * e * e - 2251.0 * e * e * e * e));
    double term2_2 = factor1 * eta * (1272.0 - 1829.0 * e * e - 538.0 * e * e * e * e);
    double factor2 = - (1.0 / 34560.0) * eta * PI * e * std::pow(x_TW, -4.0);
    double term3 = factor2 * (4538880.0 + 6876288.0 * e * e + 581208.0 * e * e * e * e + 623.0 * e * e * e * e * e * e);
    double factor3 = - (1.0 / 120.0) * eta * e * std::pow(x_TW, - 9.0 / 2.0);
    double term4_1 = factor3 * ((1.0 / 252.0) * (43837360.0 + 4258932.0 * e * e - 1211290.0 * e *e * e * e + 77535.0 * e * e * e * e * e * e));
    double term4_2 = factor3 * eta / 14.0 * (1239608.0 - 3232202.0 * e * e + 898433.0 * e * e * e * e + 13130.0 * e * e * e * e * e * e);
    double term4_3 = - factor3 * eta * eta * (9216.0 + 24353.0 * e * e + 45704.0 * e * e * e * e + 4304.0 * e * e * e * e * e * e);

    return term1 + term2_1 + term2_2 + term3 + term4_1 + term4_2 + term4_3;
}


// Tucker-Will result of dx/dtheta (eq. 2.18b in the paper)

double dx_TW_dtheta(const double& e, const double& x_TW) {
    double term1 = - (8.0 / 5.0) * eta * std::pow(x_TW, -3.0 / 2.0) * (8.0 + 7.0 * e * e);
    double factor1 = (1.0 / 15.0) * eta * std::pow(x_TW, - 5.0 / 2.0);
    double term2_1 = factor1 * ((1.0 / 14.0) * (22072.0 - 6064.0 * e * e - 1483.0 * e * e * e * e));
    double term2_2 = factor1 * 4.0 * eta * (36.0 - 127.0 * e * e - 79.0 * e * e * e * e);
    double factor2 = - (1.0 / 360.0) * eta * PI * std::pow(x_TW, - 3.0);
    double term3 = factor2 * (18432.0 + 55872.0 * e * e + 7056.0 * e * e * e * e * e - 49.0 * e * e * e * e * e * e);
    double factor3 = - (1.0 / 15.0) * eta * std::pow(x_TW, - 7.0 / 2.0);
    double term4_1 = factor3 / 756.0 * (8272600.0 + 777972.0 * e * e - 947991.0 * e * e * e * e - 4743.0 * e * e * e * e * e * e);
    double term4_2 = factor3 * eta / 84.0 * (232328.0 - 1581612.0 * e * e + 598485.0 * e * e * e * e + 6300.0 * e * e * e * e * e * e);
    double term4_3 = - factor3 * eta * eta * (384.0 + 1025.0 * e * e + 5276.0 * e * e * e * e + 632.0 * e * e * e * e * e * e);
    
    return term1 + term2_1 + term2_2 + term3 + term4_1 + term4_2 + term4_3;
}


int main() {
    return 0;
}
