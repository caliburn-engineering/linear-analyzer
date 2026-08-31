// src/auto_balance.cpp
#include "auto_balance.h"

#include <algorithm>
#include <cmath>

namespace caliburn {

namespace {
constexpr int kStates = 7;
constexpr int kLegs = 3;
}  // namespace

bool sameMechanism(const TableParams& a, const TableParams& b) {
    // A micron, not an epsilon.  The model's geometry arrives as `float`
    // sliders and the plate's is a `double` literal, so 0.300 and 0.300 differ
    // in the twelfth digit and an exact-ish comparison refuses two mechanisms
    // that are the same object.  A micron is below anything the sliders can
    // express and far above the round-trip.
    const double tol = 1e-6;
    return std::abs(a.R_ground - b.R_ground) < tol &&
           std::abs(a.R_table - b.R_table) < tol &&
           std::abs(a.L1 - b.L1) < tol &&
           std::abs(a.L2 - b.L2) < tol &&
           std::abs(a.gamma_offset - b.gamma_offset) < tol;
}

bool gainFitsCascade(const AutoBalanceDesign& d) {
    return d.K.rows() == kLegs && d.K.cols() == kStates;
}

Eigen::VectorXd cascadeState(const std::array<double, 3>& alpha_rad,
                             double home_leg_rad,
                             const Eigen::Vector4d& ball) {
    Eigen::VectorXd x(kStates);
    for (int i = 0; i < kLegs; ++i) x(i) = alpha_rad[i] - home_leg_rad;
    x.segment<4>(kLegs) = ball;
    return x;
}

LegCommand legCommand(const AutoBalanceDesign& d,
                      const std::array<double, 3>& alpha_rad,
                      const Eigen::Vector4d& ball,
                      double x_sp,
                      double y_sp) {
    LegCommand out{{d.home_leg_rad, d.home_leg_rad, d.home_leg_rad}, false};
    if (!gainFitsCascade(d)) return out;

    Eigen::VectorXd x_ref = Eigen::VectorXd::Zero(kStates);
    x_ref(3) = x_sp;
    x_ref(4) = y_sp;

    const Eigen::VectorXd u =
        -d.K * (cascadeState(alpha_rad, d.home_leg_rad, ball) - x_ref);

    for (int i = 0; i < kLegs; ++i) {
        const double raw = d.home_leg_rad + u(i);
        const double clamped = std::clamp(raw, d.alpha_min_rad, d.alpha_max_rad);
        if (clamped != raw) out.saturated = true;
        out.alpha_rad[i] = clamped;
    }
    return out;
}

std::array<double, 3> stepServos(const std::array<double, 3>& alpha_rad,
                                 const std::array<double, 3>& cmd_rad,
                                 double tau,
                                 double dt) {
    // tau <= 0 is a servo with no lag at all, which is what the plate had
    // before this existed.  It is reachable from the model panel's slider only
    // at its floor, but a degenerate tau must not divide.
    const double decay = (tau > 0.0) ? std::exp(-dt / tau) : 0.0;
    std::array<double, 3> next{};
    for (int i = 0; i < 3; ++i)
        next[i] = cmd_rad[i] + (alpha_rad[i] - cmd_rad[i]) * decay;
    return next;
}

Eigen::VectorXd defaultLqrStateWeights(int n) {
    Eigen::VectorXd q = Eigen::VectorXd::Ones(n);
    if (n == kStates) {
        q(3) = q(4) = 100.0;  // ball position is what the product is about
        q(5) = q(6) = 10.0;   // and enough velocity weight to arrive damped
    }
    return q;
}

Eigen::VectorXd defaultLqrInputWeights(int m) {
    return Eigen::VectorXd::Ones(m);
}

}  // namespace caliburn
