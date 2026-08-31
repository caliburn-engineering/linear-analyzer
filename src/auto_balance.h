// src/auto_balance.h
#pragma once

#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>
#include <cmath>

namespace caliburn {

/// A designed gain, plus everything about the operating point it was designed
/// about.  The plate simulates a *nonlinear* mechanism; the gain came from a
/// linearisation.  Every field here is part of the contract between them, and
/// the loop is only honest while all of them agree.
struct AutoBalanceDesign {
    /// 3 x 7 state feedback, u = -K(x - x_ref), u in leg-command deviations.
    Eigen::MatrixXd K;

    /// The leg angle the plant was linearised about — the `a0` parameter of
    /// the cascade preset.  States 0..2 are deviations *from this*, so a
    /// mismatch does not fail, it quietly regulates to the wrong pose.
    double home_leg_rad = M_PI / 4.0;

    /// First-order leg lag, the `tau` of the same preset.  The simulator had
    /// no servo dynamics at all before the loop was closed: the slider WAS the
    /// leg angle.  Closing u = -Kx around that is a unit-delay algebraic loop
    /// on the leg states, so the lag the model claims has to actually exist.
    double servo_tau = 0.05;

    /// Servo travel.  A command outside it is clamped, not refused — the plate
    /// really does stop there, and hiding that would make an over-aggressive
    /// tuning look better than it is.
    double alpha_min_rad = 10.0 * M_PI / 180.0;
    double alpha_max_rad = 80.0 * M_PI / 180.0;

    /// The plant the gain was designed against.  The model panel's physical
    /// sliders move it and the simulated plate does not follow, so a gain
    /// designed for longer legs — or for the moon — is not wrong-shaped, it is
    /// simply wrong, and nothing downstream could tell.  Default-constructed
    /// to zeros, so a design nobody filled in is refused rather than trusted.
    TableParams mechanism{};
    double gravity = 0.0;  ///< [m/s^2]; zero is "unset", never a real plant
};

/// Do two designs describe the same plate?  Geometry and gravity: the shape
/// the legs make and the acceleration acting on the ball are independent, and
/// getting either wrong is the same silent failure.  `alpha_min`/`alpha_max`
/// are travel limits rather than shape, and are not compared.
bool samePlant(const TableParams& a, double g_a,
               const TableParams& b, double g_b);

/// True when `d.K` has the shape the cascade plant's state feedback must have.
/// Shape is necessary, not sufficient — the caller also has to know the plant
/// IS the cascade, which no matrix dimension can tell it.
bool gainFitsCascade(const AutoBalanceDesign& d);

/// Assemble the cascade plant's state vector out of what the simulator holds:
///
///     [da1, da2, da3, x, y, x', y']
///
/// Leg deviations from the linearisation point first, then the ball's position
/// and velocity in the plate's own frame — the frame `RollingBallDynamics`
/// already integrates in, so nothing is transformed here.
Eigen::VectorXd cascadeState(const std::array<double, 3>& alpha_rad,
                             double home_leg_rad,
                             const Eigen::Vector4d& ball);

struct LegCommand {
    std::array<double, 3> alpha_rad;
    bool saturated;  ///< at least one leg hit its travel limit
};

/// The loop, evaluated once: u = -K(x - x_ref), commanded as `home + u`.
///
/// `x_ref` is `[0, 0, 0, x_sp, y_sp, 0, 0]` — and it needs no feedforward
/// term, because a ball at rest ANYWHERE on a flat plate is an equilibrium of
/// this plant.  The reference state is reachable with zero leg deviation, so
/// the regulator is already a tracker.  That is a property of the ball, not a
/// convenience: it is why the setpoint is a position and never a velocity.
///
/// A gain of the wrong shape returns the home pose and no saturation, so a
/// caller that forgets to check `gainFitsCascade` gets a flat plate rather
/// than an out-of-bounds read.
LegCommand legCommand(const AutoBalanceDesign& d,
                      const std::array<double, 3>& alpha_rad,
                      const Eigen::Vector4d& ball,
                      double x_sp,
                      double y_sp);

/// Advance the three first-order servos one step, integrated exactly:
///
///     alpha <- cmd + (alpha - cmd) * exp(-dt / tau)
///
/// Exact rather than forward-Euler because the plate runs at a fixed 60 Hz
/// against tau = 0.05 s — three steps per time constant, where Euler already
/// overshoots, and where a user dragging tau below 1/30 s would make Euler
/// oscillate and then diverge.  The exponential form cannot, at any dt.
std::array<double, 3> stepServos(const std::array<double, 3>& alpha_rad,
                                 const std::array<double, 3>& cmd_rad,
                                 double tau,
                                 double dt);

/// The weights the LQR designer opens on.
///
/// Unit weights are the textbook default, and on this plant they do not
/// balance the ball.  The rolling model carries Coulomb resistance, so a plate
/// tilted by less than asin(c_rr) = 0.573 deg cannot start the ball moving at
/// all — and a state feedback with no integral term has no way out of that
/// dead band.  It parks the ball wherever the tilt it is asking for falls
/// inside it, which under unit weights is 37 mm off centre: a loop that looks
/// broken while behaving exactly as designed.
///
/// The residual is inversely proportional to the position gain, so the fix is
/// a default that actually weights ball position — 1 mm under these.  Pinned
/// by `test_auto_balance`, both the good case and the dead band itself.
///
/// Keyed on the plant's dimensions, the only thing this layer can see.  A
/// different 7-state, 3-input plant gets the ball-balancer's opening tuning,
/// which is a poor guess rather than a wrong answer — the sliders are right
/// there, and every other shape still opens on unit weights.
Eigen::VectorXd defaultLqrStateWeights(int n);
Eigen::VectorXd defaultLqrInputWeights(int m);

}  // namespace caliburn
