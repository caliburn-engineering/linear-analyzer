// src/setpoint_path.h
#pragma once

#include <Eigen/Core>

namespace caliburn {

/// The shapes the ball can be asked to trace.
///
/// `Fixed` is the setpoint the balance loop has always had — a point the
/// visitor drags.  The rest move it, which asks a different question of the
/// controller: not "can you hold a ball still" but "can you make it go where
/// I say".  A circle the ball actually follows is self-evidently a working
/// controller to someone who does not read a pole-zero map, and the corner it
/// rounds off on the square is self-evidently a bandwidth limit.
///
/// See [#24](https://github.com/caliburn-engineering/caliburn/issues/24).
enum class PathShape { Fixed, Circle, Square, Triangle };

/// A closed path for the ball setpoint to run round, in the plate's own frame.
struct SetpointPath {
    PathShape shape = PathShape::Fixed;

    /// Circumradius: the distance from the centre to a corner, and for the
    /// circle simply its radius.  Sized from the corner rather than the edge
    /// so that every shape at the same setting reaches equally far out, which
    /// is what a visitor comparing them expects.
    double radius_m = 0.12;

    /// Seconds for one lap.  Time rather than speed because it is what makes
    /// the shapes comparable: at the same period a triangle and a circle come
    /// round together, and the triangle is simply travelling faster to do it.
    double period_s = 10.0;
};

/// The fastest the setpoint may be driven, and the furthest out it may go.
///
/// Both measured, both against the shipped tuning with velocity feedforward.
/// The ball trails the setpoint and overshoots the corners, so the path's own
/// radius is not the radius the BALL reaches: at 180 mm and 283 mm/s the ball
/// swings to 195 mm, and the plate loses it entirely beyond about 400 mm/s —
/// a 250 mm circle at a two-second lap throws it clean off.
///
/// So the sliders are bounded rather than left to find that out.  A demo whose
/// controls include a setting that breaks it is not offering a choice, it is
/// offering a trap.
inline constexpr double kMaxSetpointSpeed = 0.25;   ///< [m/s]
inline constexpr double kMaxPathRadius = 0.18;      ///< [m]

/// And a floor on the lap time whatever the size.
///
/// A speed cap alone is the wrong bound at small radii: 250 mm/s round a 20 mm
/// circle is a half-second lap, which is an angular rate the plate cannot
/// follow however short the distance.  Two constraints because there are two
/// ways to ask for something impossible — go too far too fast, or go round too
/// often — and neither implies the other.
inline constexpr double kMinLapSeconds = 2.0;

/// The shortest lap that keeps a path under `kMaxSetpointSpeed`.
///
/// The larger of the two bounds: `length / kMaxSetpointSpeed`, which binds on
/// the big paths, and `kMinLapSeconds`, which binds on the small ones.
double minPeriod(const SetpointPath& p);

/// How far round the lap the setpoint is: 0 at the start, 1 back at the start.
///
/// **Phase, not time, and this is the whole of it.**  The setpoint used to be
/// evaluated at `t / period_s`, which reads as a pure function of the clock
/// and is not one: changing `period_s` moves that quantity by `t dT / T^2`,
/// and `t` is the entire time the simulation has been running.  Measured, at
/// 100 s, a lap change of 10.0 -> 9.5 s moved the setpoint **170 degrees** —
/// so a nudge of the slider teleported the target across the plate and the
/// loop hauled the ball after it.  The jump grew with run time and wrapped, so
/// any nudge could land anywhere.
///
/// Accumulating the phase instead makes a lap change alter only the RATE from
/// that moment on, which is what the slider says it does.  See #24.
///
/// A degenerate period does not advance rather than dividing by zero.
double advancePhase(double phase, double dt, double period_s);

/// Where the setpoint is at phase `u`.
///
/// The polygons are traversed at constant SPEED, not constant angle: a corner
/// is a change of direction, and slowing into it would hide exactly the
/// behaviour the cornered shapes exist to show.
///
/// Linear in `radius_m` for every shape, which is what makes the size slider
/// safe to drag: the setpoint slides straight out along the ray it was already
/// on, same angle, bigger shape.
Eigen::Vector2d pathPoint(const SetpointPath& p, double u);

/// How fast the setpoint is moving at phase `u`, and where it is going.
///
/// The one place `period_s` still enters directly, because it must: the phase
/// says WHERE round the lap, and the period says how fast that is being
/// walked.  So the reference velocity DOES step when the lap slider moves,
/// while the position does not — which is correct, and is exactly what the
/// slider was asking for.
///
/// Undefined for an instant at each corner, where the path's velocity is
/// genuinely discontinuous; the value returned there is the edge being left.
/// That is honest — a corner IS a step in the reference velocity, and it is
/// the reason the ball rounds one.
Eigen::Vector2d pathVelocity(const SetpointPath& p, double u);

/// The path's total length, for drawing it and for reasoning about speed.
double pathLength(const SetpointPath& p);

/// Points around one lap, for drawing.  A circle gets `samples` of them; a
/// polygon gets its corners, because a corner drawn as a chord of samples is
/// a corner drawn wrong.
void pathOutline(const SetpointPath& p, int samples,
                 Eigen::Matrix<double, 2, Eigen::Dynamic>* out);

}  // namespace caliburn
