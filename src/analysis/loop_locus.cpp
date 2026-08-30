// src/analysis/loop_locus.cpp
#include "loop_locus.h"
#include "system_connect.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>

namespace caliburn {

std::vector<RootLocusPoint> computeLoopLocus(
    const LinearSystem& plant,
    const std::vector<Loop>& loops,
    LoopKind kind,
    int output_i, int input_j,
    double kappa_min, double kappa_max, int num_points) {
    std::vector<RootLocusPoint> result;
    if (loops.empty() || kappa_min <= 0.0 || kappa_max <= kappa_min ||
        num_points < 2)
        return result;

    // Same liveness predicate the Bode panel uses — one rule, two panels.
    // While inspecting an off-loop plant channel (the Ball-Balancer pairing
    // demo) the locus is correctly blank.
    bool live = false;
    for (const auto& l : loops)
        if (l.out == output_i && l.in == input_j) { live = true; break; }
    if (!live) return result;

    const double log_min = std::log10(kappa_min);
    const double log_max = std::log10(kappa_max);
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        const double kappa = std::pow(
            10.0, log_min + (log_max - log_min) * step / (num_points - 1));

        // kappa * C_k(s).  Scaling all three PID gains together is exactly
        // that, since kappa*(Kp + Ki/s + Kd*s/(tau_f*s+1)) leaves tau_f alone.
        //
        // kappa != 0 never changes WHICH gains are zero, so the minimal
        // realization emits the identical state set at every point — the
        // closed-loop dimension is constant by construction, which is what
        // makes matchPoles and the panel's fixed-width drawing loop safe.
        // kappa = 0 would collapse the loop's states and matchPoles would read
        // prev[i] out of bounds when the pole count grew again.
        std::vector<Loop> scaled = loops;
        for (auto& l : scaled) {
            if (l.out != output_i || l.in != input_j) continue;
            l.pid.Kp = static_cast<float>(l.pid.Kp * kappa);
            l.pid.Ki = static_cast<float>(l.pid.Ki * kappa);
            l.pid.Kd = static_cast<float>(l.pid.Kd * kappa);
            l.leadlag.Kc = static_cast<float>(l.leadlag.Kc * kappa);
        }

        // Literally the production path, so the locus cannot disagree with the
        // closed-loop pole plot drawn beside it.  buildLoopController's
        // contract is untouched.
        const LinearSystem ctrl = buildLoopController(scaled, kind, plant);
        const LinearSystem closed =
            feedbackConnect(seriesConnect(ctrl, plant));

        const int n = closed.states();
        std::vector<std::complex<double>> poles(n);
        if (n > 0) {
            Eigen::EigenSolver<Eigen::MatrixXd> es(closed.A, false);
            for (int i = 0; i < n; ++i) poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) matchPoles(poles, result.back().poles);
        result.push_back({kappa, poles});
    }
    return result;
}

double loopGainMargin(const std::vector<RootLocusPoint>& locus) {
    auto maxRe = [](const RootLocusPoint& pt) {
        double m = -std::numeric_limits<double>::infinity();
        for (const auto& p : pt.poles) m = std::max(m, p.real());
        return m;
    };
    for (std::size_t k = 1; k < locus.size(); ++k) {
        const double a = maxRe(locus[k - 1]);
        const double b = maxRe(locus[k]);
        if (a < 0.0 && b >= 0.0) {
            const double t = -a / (b - a);
            return locus[k - 1].gain + t * (locus[k].gain - locus[k - 1].gain);
        }
    }
    return -1.0;  // no crossing in range — never report a number
}

}  // namespace caliburn
