#include "table_kinematics.h"
#include <stdexcept>

// ============================================================================
// Construction
// ============================================================================

TableKinematics::TableKinematics(const TableParams& params) : params_(params) {
    if (params.R_ground <= 0 || params.R_table <= 0 ||
        params.L1 <= 0 || params.L2 <= 0) {
        throw std::invalid_argument("All geometric parameters must be positive");
    }
    if (params.alpha_min >= params.alpha_max) {
        throw std::invalid_argument("alpha_min must be less than alpha_max");
    }
}

// ============================================================================
// Geometry
// ============================================================================

double TableKinematics::beta(int i) const {
    return i * 2.0 * M_PI / 3.0;
}

double TableKinematics::gamma(int i) const {
    return i * 2.0 * M_PI / 3.0 + params_.gamma_offset;
}

Eigen::Vector3d TableKinematics::ground_point(int i) const {
    double b = beta(i);
    return Eigen::Vector3d(params_.R_ground * std::cos(b),
                           params_.R_ground * std::sin(b),
                           0.0);
}

Eigen::Vector3d TableKinematics::radial_dir(int i) const {
    double b = beta(i);
    return Eigen::Vector3d(std::cos(b), std::sin(b), 0.0);
}

Eigen::Vector3d TableKinematics::tangent_dir(int i) const {
    double b = beta(i);
    return Eigen::Vector3d(-std::sin(b), std::cos(b), 0.0);
}

// ============================================================================
// Knee Positions (simplified: L1 in radial-vertical plane)
// ============================================================================

Eigen::Vector3d TableKinematics::knee_position(int i, double alpha_i) const {
    Eigen::Vector3d G = ground_point(i);
    Eigen::Vector3d r = radial_dir(i);
    Eigen::Vector3d z_hat(0.0, 0.0, 1.0);

    return G + params_.L1 * (std::cos(alpha_i) * r + std::sin(alpha_i) * z_hat);
}

// ============================================================================
// Table Attachment Points
// ============================================================================

Eigen::Vector3d TableKinematics::table_point_local(int i) const {
    double g = gamma(i);
    return Eigen::Vector3d(params_.R_table * std::cos(g),
                           params_.R_table * std::sin(g),
                           0.0);
}

Eigen::Matrix3d TableKinematics::table_rotation(double phi, double theta) const {
    // R = Ry(theta) * Rx(phi)
    double cp = std::cos(phi),   sp = std::sin(phi);
    double ct = std::cos(theta), st = std::sin(theta);

    Eigen::Matrix3d R;
    R <<  ct,      sp * st,    cp * st,
          0.0,     cp,        -sp,
         -st,      sp * ct,    cp * ct;
    return R;
}

Eigen::Vector3d TableKinematics::table_point_world(int i, const TablePose& pose) const {
    Eigen::Matrix3d R = table_rotation(pose.phi, pose.theta);
    Eigen::Vector3d p_local = table_point_local(i);
    Eigen::Vector3d center(0.0, 0.0, pose.z_c);
    return R * p_local + center;
}

// ============================================================================
// Inverse Kinematics (closed-form, per leg)
// ============================================================================

std::optional<double> TableKinematics::inverse_kinematics_leg(int i, const TablePose& pose) const {
    // Table attachment point in world frame
    Eigen::Vector3d P = table_point_world(i, pose);

    // Ground attachment
    Eigen::Vector3d G = ground_point(i);

    double px = P.x(), py = P.y(), pz = P.z();
    double b = beta(i);
    double cb = std::cos(b), sb = std::sin(b);

    // Projection of P onto radial direction at leg i
    double A = px * cb + py * sb;
    double B = pz;

    // ||P - K||^2 = L2^2 expanded:
    // P^2 + R_g^2 + L1^2 - L2^2 + 2*L1*(R_g - A)*cos(alpha) - 2*L1*B*sin(alpha) - 2*R_g*A = 0
    //
    // But we also have the tangential component of P:
    // P_tangential = -px*sin(b) + py*cos(b)
    // This contributes to ||P - K||^2 but is independent of alpha.

    // Full expansion of ||P - K||^2:
    // Let K = G + L1*(cos(a)*r + sin(a)*z)
    // d = P - K
    // d_radial = (A - R_g - L1*cos(a))    [radial component]
    // d_tang   = -px*sin(b) + py*cos(b)   [tangential component, independent of alpha]
    // d_z      = B - L1*sin(a)            [vertical component]
    //
    // ||d||^2 = d_radial^2 + d_tang^2 + d_z^2 = L2^2

    double P_tang = -px * sb + py * cb;
    double L2_eff_sq = params_.L2 * params_.L2 - P_tang * P_tang;

    if (L2_eff_sq < 0.0) {
        // Table point is too far in the tangential direction — unreachable
        return std::nullopt;
    }

    // Now solve in the radial-vertical plane:
    // (A - R_g - L1*cos(a))^2 + (B - L1*sin(a))^2 = L2_eff_sq
    //
    // Expand:
    // (A - R_g)^2 - 2*L1*(A - R_g)*cos(a) + L1^2*cos^2(a)
    //   + B^2 - 2*L1*B*sin(a) + L1^2*sin^2(a) = L2_eff_sq
    //
    // L1^2 + (A-R_g)^2 + B^2 - 2*L1*(A-R_g)*cos(a) - 2*L1*B*sin(a) = L2_eff_sq

    double Ar = A - params_.R_ground;
    double C = params_.L1 * params_.L1 + Ar * Ar + B * B - L2_eff_sq;

    // C - 2*L1*Ar*cos(a) - 2*L1*B*sin(a) = 0
    // => 2*L1*Ar*cos(a) + 2*L1*B*sin(a) = C
    // => a_coeff*cos(a) + b_coeff*sin(a) = d

    double a_coeff = 2.0 * params_.L1 * Ar;
    double b_coeff = 2.0 * params_.L1 * B;
    double d = C;

    double norm_ab = std::sqrt(a_coeff * a_coeff + b_coeff * b_coeff);

    if (norm_ab < 1e-12) {
        // Degenerate — L1 is zero or P is at G
        return std::nullopt;
    }

    double ratio = d / norm_ab;
    if (std::abs(ratio) > 1.0 + 1e-9) {
        // No solution — target unreachable
        return std::nullopt;
    }
    ratio = std::clamp(ratio, -1.0, 1.0);

    // alpha = atan2(b_coeff, a_coeff) ± acos(ratio)
    double base_angle = std::atan2(b_coeff, a_coeff);
    double delta = std::acos(ratio);

    // Two candidate solutions
    double alpha1 = base_angle + delta;
    double alpha2 = base_angle - delta;

    // Normalize to [0, 2*pi) for comparison
    auto normalize = [](double a) {
        while (a < 0) a += 2.0 * M_PI;
        while (a >= 2.0 * M_PI) a -= 2.0 * M_PI;
        return a;
    };

    alpha1 = normalize(alpha1);
    alpha2 = normalize(alpha2);

    // Select the solution within servo limits
    auto in_range = [&](double a) {
        return a >= params_.alpha_min - 1e-9 && a <= params_.alpha_max + 1e-9;
    };

    bool ok1 = in_range(alpha1);
    bool ok2 = in_range(alpha2);

    if (ok1 && ok2) {
        // Both valid — pick the one closer to home (45°)
        double home = M_PI / 4.0;
        return (std::abs(alpha1 - home) <= std::abs(alpha2 - home)) ? alpha1 : alpha2;
    } else if (ok1) {
        return alpha1;
    } else if (ok2) {
        return alpha2;
    }

    // Neither solution in range — infeasible
    return std::nullopt;
}

IKResult TableKinematics::inverse_kinematics(const TablePose& pose) const {
    IKResult result;
    result.feasible = true;

    for (int i = 0; i < 3; ++i) {
        auto alpha_opt = inverse_kinematics_leg(i, pose);
        if (alpha_opt) {
            result.alpha[i] = std::clamp(*alpha_opt, params_.alpha_min, params_.alpha_max);
        } else {
            result.alpha[i] = M_PI / 4.0; // fallback
            result.feasible = false;
        }
    }
    return result;
}

bool TableKinematics::is_feasible(const TablePose& pose) const {
    for (int i = 0; i < 3; ++i) {
        if (!inverse_kinematics_leg(i, pose)) return false;
    }
    return true;
}

// ============================================================================
// Forward Kinematics (Newton-Raphson)
// ============================================================================

Eigen::Vector3d TableKinematics::fk_residual(const std::array<double, 3>& alpha,
                                              const TablePose& pose) const {
    Eigen::Vector3d f;
    for (int i = 0; i < 3; ++i) {
        Eigen::Vector3d P = table_point_world(i, pose);
        Eigen::Vector3d K = knee_position(i, alpha[i]);
        Eigen::Vector3d d = P - K;
        f(i) = d.squaredNorm() - params_.L2 * params_.L2;
    }
    return f;
}

Eigen::Matrix3d TableKinematics::dR_dphi(double phi, double theta) const {
    double cp = std::cos(phi),   sp = std::sin(phi);
    double ct = std::cos(theta), st = std::sin(theta);

    Eigen::Matrix3d dR;
    dR <<  0.0,     cp * st,   -sp * st,
           0.0,    -sp,        -cp,
           0.0,     cp * ct,   -sp * ct;
    return dR;
}

Eigen::Matrix3d TableKinematics::dR_dtheta(double phi, double theta) const {
    double cp = std::cos(phi),   sp = std::sin(phi);
    double ct = std::cos(theta), st = std::sin(theta);

    Eigen::Matrix3d dR;
    dR << -st,      sp * ct,    cp * ct,
           0.0,     0.0,        0.0,
          -ct,     -sp * st,   -cp * st;
    return dR;
}

Eigen::Matrix3d TableKinematics::fk_jacobian(const std::array<double, 3>& alpha,
                                              const TablePose& pose) const {
    Eigen::Matrix3d J;

    Eigen::Matrix3d dRdp = dR_dphi(pose.phi, pose.theta);
    Eigen::Matrix3d dRdt = dR_dtheta(pose.phi, pose.theta);

    for (int i = 0; i < 3; ++i) {
        Eigen::Vector3d P = table_point_world(i, pose);
        Eigen::Vector3d K = knee_position(i, alpha[i]);
        Eigen::Vector3d d = P - K;
        Eigen::Vector3d p_local = table_point_local(i);

        // dP/dphi = dR/dphi * p_local
        Eigen::Vector3d dP_dphi = dRdp * p_local;
        // dP/dtheta = dR/dtheta * p_local
        Eigen::Vector3d dP_dtheta = dRdt * p_local;
        // dP/dz_c = [0, 0, 1]
        Eigen::Vector3d dP_dzc(0.0, 0.0, 1.0);

        // J(i, j) = 2 * d^T * dP/dx_j
        J(i, 0) = 2.0 * d.dot(dP_dphi);
        J(i, 1) = 2.0 * d.dot(dP_dtheta);
        J(i, 2) = 2.0 * d.dot(dP_dzc);
    }
    return J;
}

FKResult TableKinematics::forward_kinematics(const std::array<double, 3>& alpha,
                                              const TablePose& x0,
                                              int max_iter,
                                              double tol) const {
    FKResult result;
    result.pose = x0;
    result.converged = false;
    result.iterations = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::Vector3d f = fk_residual(alpha, result.pose);
        result.residual_norm = f.norm();
        result.iterations = iter + 1;

        if (result.residual_norm < tol) {
            result.converged = true;
            return result;
        }

        Eigen::Matrix3d J = fk_jacobian(alpha, result.pose);

        // Solve J * dx = -f
        Eigen::Vector3d dx = J.colPivHouseholderQr().solve(-f);

        // Backtracking line search: ensure residual decreases
        double step = 1.0;
        TablePose candidate = result.pose;
        for (int ls = 0; ls < 10; ++ls) {
            candidate.phi   = result.pose.phi   + step * dx(0);
            candidate.theta = result.pose.theta + step * dx(1);
            candidate.z_c   = result.pose.z_c   + step * dx(2);

            double new_residual = fk_residual(alpha, candidate).norm();
            if (new_residual < result.residual_norm) break;
            step *= 0.5;
        }

        result.pose = candidate;
    }

    // Final residual check
    Eigen::Vector3d f = fk_residual(alpha, result.pose);
    result.residual_norm = f.norm();
    result.converged = (result.residual_norm < tol);

    return result;
}

// ============================================================================
// Velocity Jacobian
// ============================================================================

Eigen::Vector3d TableKinematics::dKnee_dAlpha(int i, double alpha_i) const {
    // K_i = G_i + L1 * (cos(alpha) * r_hat + sin(alpha) * z_hat)
    // dK/dalpha = L1 * (-sin(alpha) * r_hat + cos(alpha) * z_hat)
    Eigen::Vector3d r = radial_dir(i);
    Eigen::Vector3d z_hat(0.0, 0.0, 1.0);
    return params_.L1 * (-std::sin(alpha_i) * r + std::cos(alpha_i) * z_hat);
}

Eigen::Matrix3d TableKinematics::constraint_jacobian_alpha(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    // J_alpha is diagonal: J_alpha(i,i) = df_i/d(alpha_i)
    // f_i = ||P_i - K_i||^2 - L2^2
    // df_i/dalpha_i = 2 * (P_i - K_i)^T * (-dK_i/dalpha_i)
    Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
    for (int i = 0; i < 3; ++i) {
        Eigen::Vector3d d = table_point_world(i, pose) - knee_position(i, alpha[i]);
        Eigen::Vector3d dK = dKnee_dAlpha(i, alpha[i]);
        J(i, i) = -2.0 * d.dot(dK);
    }
    return J;
}

Eigen::Matrix3d TableKinematics::velocity_jacobian(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    // Implicit differentiation of f(alpha, pose) = 0:
    //   J_pose * d(pose)/dt + J_alpha * d(alpha)/dt = 0
    //   d(pose)/dt = -J_pose^{-1} * J_alpha * d(alpha)/dt
    // So J_v = -J_pose^{-1} * J_alpha
    Eigen::Matrix3d Jp = fk_jacobian(alpha, pose);
    Eigen::Matrix3d Ja = constraint_jacobian_alpha(alpha, pose);
    return -Jp.colPivHouseholderQr().solve(Ja);
}

// ============================================================================
// Singularity Analysis
// ============================================================================

double TableKinematics::manipulability(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    Eigen::Matrix3d Jv = velocity_jacobian(alpha, pose);
    return std::sqrt(std::abs((Jv * Jv.transpose()).determinant()));
}

double TableKinematics::condition_number(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    Eigen::Matrix3d Jv = velocity_jacobian(alpha, pose);
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(Jv);
    auto sv = svd.singularValues();
    if (sv(2) < 1e-15) return std::numeric_limits<double>::infinity();
    return sv(0) / sv(2);
}

double TableKinematics::det_J_pose(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    return fk_jacobian(alpha, pose).determinant();
}

double TableKinematics::det_J_alpha(
    const std::array<double, 3>& alpha, const TablePose& pose) const {
    return constraint_jacobian_alpha(alpha, pose).determinant();
}

// ============================================================================
// Home Position
// ============================================================================

TablePose TableKinematics::home_pose(double alpha_home) const {
    // All servos at the same angle => table is level, centered on z-axis
    // Knee height: L1 * sin(alpha_home)
    // Knee radial position: R_g + L1 * cos(alpha_home)
    // Table point radial position: R_t (when phi=0, theta=0)
    // Horizontal distance knee to table point: |knee_radial - R_t|
    // Vertical distance: sqrt(L2^2 - horizontal^2)

    double knee_r = params_.R_ground + params_.L1 * std::cos(alpha_home);
    double knee_z = params_.L1 * std::sin(alpha_home);
    double dr = knee_r - params_.R_table; // radial distance (can be negative if table is wider)
    double dz_sq = params_.L2 * params_.L2 - dr * dr;

    if (dz_sq < 0.0) {
        // Geometry doesn't close — L2 too short. Return best guess.
        return {0.0, 0.0, knee_z};
    }

    double z_c = knee_z + std::sqrt(dz_sq);

    return {0.0, 0.0, z_c};
}
