#include <cmath>
#include <array>
#include <vector>

// Physical constants and parameters
struct PhysicalParams {
    double G;      // Gravitational constant
    double M;      // Total mass
    double eta;    // Symmetric mass ratio (eta = m1*m2/(m1+m2)^2)
    double phi;    // Phase angle
};

// State variables for binary evolution
struct BinaryState {
    double p;        // Orbital parameter (semi-latus rectum)
    double alpha;    // Spin parameter 1
    double beta;     // Spin parameter 2
};

// Container for RHS terms: {dp/dtheta, dalpha/dtheta, dbeta/dtheta}
using SecularRHS = std::array<double, 3>;

// ============================================================================
// 1 PN ORDER TERMS
// ============================================================================

SecularRHS secular_1PN(const BinaryState& state, const PhysicalParams& params) {
    const double& p = state.p;
    const double& alpha = state.alpha;
    const double& beta = state.beta;
    const double& eta = params.eta;
    const double GM = params.G * params.M;
    
    SecularRHS rhs;
    rhs[0] = 0.0;
    rhs[1] = -(3.0 * GM * beta) / p;
    rhs[2] = (3.0 * GM * alpha) / p;
    
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
    
    SecularRHS rhs;
    
    double bracket = -10.0 + beta * beta - 4.0 * eta + 10.0 * beta * beta * eta + 
                     alpha * alpha * (1.0 + 10.0 * eta);
    
    rhs[0] = 0.0;
    rhs[1] = -(3.0 * GM2 * beta * bracket) / (4.0 * p * p);
    rhs[2] = (3.0 * GM2 * alpha * bracket) / (4.0 * p * p);
    
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
        auto rhs = secular_2_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 4) {  // 3.5 PN (order 4 in array indexing)
        auto rhs = secular_3_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    if (max_PN_order >= 5) {  // 4.5 PN (order 5 in array indexing)
        auto rhs = secular_4_5PN(state, params);
        total_rhs[0] += rhs[0]; total_rhs[1] += rhs[1]; total_rhs[2] += rhs[2];
    }
    
    return total_rhs;
}

// ============================================================================
// NUMERICAL INTEGRATION SETUP
// ============================================================================

class SecularEquationIntegrator {
public:
    PhysicalParams params;
    int max_PN_order;
    
    SecularEquationIntegrator(double G, double M, double eta, int pn_order = 5)
        : max_PN_order(pn_order) {
        params.G = G;
        params.M = M;
        params.eta = eta;
        params.phi = 0.0;
    }
    
    // RK4 step function
    BinaryState rk4_step(const BinaryState& state, double dtheta) {
        // k1
        SecularRHS k1 = compute_secular_RHS(state, params, max_PN_order);
        
        // k2
        BinaryState state2 = state;
        state2.p = state.p + 0.5 * k1[0] * dtheta;
        state2.alpha = state.alpha + 0.5 * k1[1] * dtheta;
        state2.beta = state.beta + 0.5 * k1[2] * dtheta;
        SecularRHS k2 = compute_secular_RHS(state2, params, max_PN_order);
        
        // k3
        BinaryState state3 = state;
        state3.p = state.p + 0.5 * k2[0] * dtheta;
        state3.alpha = state.alpha + 0.5 * k2[1] * dtheta;
        state3.beta = state.beta + 0.5 * k2[2] * dtheta;
        SecularRHS k3 = compute_secular_RHS(state3, params, max_PN_order);
        
        // k4
        BinaryState state4 = state;
        state4.p = state.p + k3[0] * dtheta;
        state4.alpha = state.alpha + k3[1] * dtheta;
        state4.beta = state.beta + k3[2] * dtheta;
        SecularRHS k4 = compute_secular_RHS(state4, params, max_PN_order);
        
        // Combine
        BinaryState new_state = state;
        new_state.p += (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]) * dtheta / 6.0;
        new_state.alpha += (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]) * dtheta / 6.0;
        new_state.beta += (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]) * dtheta / 6.0;
        
        return new_state;
    }
    
    // Integrate from theta_initial to theta_final
    std::vector<BinaryState> integrate(const BinaryState& initial_state, 
                                       double theta_initial, 
                                       double theta_final, 
                                       double dtheta) {
        std::vector<BinaryState> trajectory;
        BinaryState current = initial_state;
        trajectory.push_back(current);
        
        double theta = theta_initial;
        if (dtheta > 0) {
            while (theta < theta_final) {
                current = rk4_step(current, dtheta);
                theta += dtheta;
                trajectory.push_back(current);
            }
        } else {
            while (theta > theta_final) {
                current = rk4_step(current, dtheta);
                theta += dtheta;
                trajectory.push_back(current);
            }
        }
        
        return trajectory;
    }
};
