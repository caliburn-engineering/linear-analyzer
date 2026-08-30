// src/ball_sim.h
#pragma once

#include "rolling_dynamics.h"
#include "table_kinematics.h"

#include <Eigen/Core>

namespace caliburn {

/// Tilt angles as the rolling-ball model names them: `alpha` is the tilt that
/// drives the ball along y, `beta` the tilt that drives it along x.
struct BallTilt {
    double alpha;
    double beta;
};

/// Read the rolling model's tilt inputs off a table pose.
///
/// The two halves disagree, and this function is where they are reconciled:
/// `TableKinematics` uses R = Ry(theta) * Rx(phi), so the plate's local +Y edge
/// rises with phi while its local +X edge falls with theta.  The rolling model
/// accelerates the ball by +g*sin(alpha) along y and +g*sin(beta) along x.  So
/// theta carries straight through and phi has to flip.
BallTilt ballTiltFromPose(const TablePose& pose);

/// Advance the ball one RK4 step over a plate held at `pose`.
/// State is [x, y, vx, vy] in the plate's own frame.
Eigen::Vector4d stepBall(const RollingBallDynamics& dynamics,
                         const Eigen::Vector4d& state,
                         const TablePose& pose,
                         double dt);

/// True while the ball's contact patch is still over a disc plate.
/// The plate is a disc, not the square that `RollingBallDynamics::on_plate`
/// assumes, so the bound is radial and inset by the ball's own radius.
bool ballOnPlate(const Eigen::Vector4d& state,
                 double plate_radius,
                 double ball_radius);

}  // namespace caliburn
