// src/attract_mode.cpp
#include "attract_mode.h"

namespace caliburn {

Eigen::Vector4d attractStart(const AttractSchedule& s) {
    return Eigen::Vector4d(s.start_x, s.start_y, 0.0, 0.0);
}

int attractKicksBy(const AttractSchedule& s, double t) {
    if (t < s.first_kick_s) return 0;
    // A zero or negative period is one kick.  Not defensiveness about a value
    // nothing sets today: the division below would be by zero, and casting the
    // resulting infinity to `int` is undefined behaviour rather than a large
    // number — the kind of bug that survives every test and appears once, in a
    // release build, on somebody else's machine.
    if (!(s.period_s > 0.0)) return 1;
    return 1 + static_cast<int>((t - s.first_kick_s) / s.period_s);
}

Eigen::Vector4d attractKick(const AttractSchedule& s, int n) {
    const double dir = s.turn_rad * n;
    return Eigen::Vector4d(0.0, 0.0,
                           s.speed * std::cos(dir),
                           s.speed * std::sin(dir));
}

}  // namespace caliburn
