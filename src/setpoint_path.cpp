// src/setpoint_path.cpp
#include "setpoint_path.h"

#include <algorithm>
#include <cmath>

namespace caliburn {
namespace {

/// Corners of the regular polygon a shape describes, or 0 for the others.
int corners(PathShape s) {
    switch (s) {
        case PathShape::Triangle: return 3;
        case PathShape::Square:   return 4;
        default:                  return 0;
    }
}

/// Corner `k` of an `n`-gon on the circumradius.
///
/// The first corner is put at the top rather than at +x so that the square
/// reads as a square: on this plate the visitor is looking at a triad of legs
/// at 120 degrees, and a square rotated 45 degrees against it looks like a
/// mistake rather than a choice.
Eigen::Vector2d corner(int k, int n, double r) {
    const double a = M_PI / 2.0 + 2.0 * M_PI * k / n;
    return Eigen::Vector2d(r * std::cos(a), r * std::sin(a));
}

/// Position along the perimeter, as a distance and the edge it falls on.
struct Along {
    int edge;
    double frac;   ///< 0..1 along that edge
};

Along along(const SetpointPath& p, double t, int n) {
    const double side = 2.0 * p.radius_m * std::sin(M_PI / n);
    const double lap = (p.period_s > 0.0) ? std::fmod(t / p.period_s, 1.0) : 0.0;
    const double u = (lap < 0.0 ? lap + 1.0 : lap) * n;   // edges travelled
    int edge = static_cast<int>(u) % n;
    double frac = u - std::floor(u);
    (void)side;
    return {edge, frac};
}

}  // namespace

Eigen::Vector2d pathPoint(const SetpointPath& p, double t) {
    switch (p.shape) {
        case PathShape::Fixed:
            return Eigen::Vector2d::Zero();
        case PathShape::Circle: {
            const double w = (p.period_s > 0.0) ? 2.0 * M_PI * t / p.period_s : 0.0;
            return Eigen::Vector2d(p.radius_m * std::cos(w), p.radius_m * std::sin(w));
        }
        default: break;
    }
    const int n = corners(p.shape);
    const Along a = along(p, t, n);
    const Eigen::Vector2d from = corner(a.edge, n, p.radius_m);
    const Eigen::Vector2d to = corner(a.edge + 1, n, p.radius_m);
    return from + a.frac * (to - from);
}

Eigen::Vector2d pathVelocity(const SetpointPath& p, double t) {
    if (p.shape == PathShape::Fixed || p.period_s <= 0.0)
        return Eigen::Vector2d::Zero();
    if (p.shape == PathShape::Circle) {
        const double w = 2.0 * M_PI / p.period_s;
        return Eigen::Vector2d(-p.radius_m * w * std::sin(w * t),
                                p.radius_m * w * std::cos(w * t));
    }
    const int n = corners(p.shape);
    const Along a = along(p, t, n);
    const Eigen::Vector2d from = corner(a.edge, n, p.radius_m);
    const Eigen::Vector2d to = corner(a.edge + 1, n, p.radius_m);
    // One edge per period/n, covered at constant speed.
    return (to - from) * (n / p.period_s);
}

double pathLength(const SetpointPath& p) {
    switch (p.shape) {
        case PathShape::Fixed:  return 0.0;
        case PathShape::Circle: return 2.0 * M_PI * p.radius_m;
        default: break;
    }
    const int n = corners(p.shape);
    return n * 2.0 * p.radius_m * std::sin(M_PI / n);
}

double minPeriod(const SetpointPath& p) {
    const double len = pathLength(p);
    if (len <= 0.0) return 0.0;
    return std::max(len / kMaxSetpointSpeed, kMinLapSeconds);
}

void pathOutline(const SetpointPath& p, int samples,
                 Eigen::Matrix<double, 2, Eigen::Dynamic>* out) {
    if (!out) return;
    if (p.shape == PathShape::Fixed) { out->resize(2, 0); return; }
    if (p.shape == PathShape::Circle) {
        out->resize(2, samples + 1);
        for (int i = 0; i <= samples; ++i) {
            const double a = 2.0 * M_PI * i / samples;
            out->col(i) = Eigen::Vector2d(p.radius_m * std::cos(a),
                                          p.radius_m * std::sin(a));
        }
        return;
    }
    const int n = corners(p.shape);
    out->resize(2, n + 1);
    for (int k = 0; k <= n; ++k) out->col(k) = corner(k, n, p.radius_m);
}

}  // namespace caliburn
