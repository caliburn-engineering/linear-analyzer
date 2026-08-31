// src/attract_mode.h
#pragma once

#include <Eigen/Core>
#include <cmath>

namespace caliburn {

/// The demo driving itself for a visitor who has not arrived yet.
///
/// A motionless canvas is indistinguishable from a broken build, and a ball
/// sitting balanced reads as a still image of one.  What says "control system"
/// rather than "animation" is *recovery*: something disturbs the ball, and the
/// loop puts it back.  So the application opens with the loop already closed,
/// the ball already off centre, and a kick arriving every few seconds until
/// somebody touches the page.
///
/// The schedule is a pure function of elapsed time, and nothing here reaches
/// for a random number.  That is deliberate: the acceptance criteria are
/// claims about what a visitor sees — the disturbance is visible, the recovery
/// is visible, the ball never leaves the plate — and a schedule that differed
/// run to run could only be checked by watching it.  Written this way,
/// `test_attract_mode` runs the whole minute against the nonlinear plate.
///
/// The time it is a function of is the **simulation's**, not the wall's.  The
/// plate advances by a fixed 1/60 per frame whatever the frame rate, so on a
/// device that cannot hold 60 fps the whole demo plays slow — and pacing the
/// kicks by a wall clock would then land the next one on a recovery that has
/// not finished, which is precisely the reading the ticket asks for and would
/// destroy.  A slow device gets the demo late rather than incoherent.  That is
/// also why the first kick comes well inside the settling time's own margin
/// (see `first_kick_s`) rather than being spaced like the rest.
///
/// See [#17](https://github.com/caliburn-engineering/caliburn/issues/17).
struct AttractSchedule {
    /// Where the ball is put at t = 0, at rest, in the plate frame.  The
    /// opening motion is a recovery from a displacement rather than a kick:
    /// at t = 0 there is nothing on screen yet for a kick to be a change
    /// *from*, and a ball that starts off centre is already moving on the
    /// first frame that draws.
    double start_x = 0.06;   ///< [m]
    double start_y = -0.04;  ///< [m]

    /// The first kick, and the gap between kicks.  Both clear the ~1.8 s the
    /// shipped tuning takes to settle the opening displacement, so a recovery
    /// always finishes before the next disturbance arrives — a kick landing
    /// mid-recovery reads as noise rather than as a loop rejecting anything,
    /// and the still moment between them is what makes the next one legible.
    ///
    /// The first is deliberately tighter than the rest.  It is the only one a
    /// visitor is guaranteed to be present for, and on a device running at
    /// half frame rate it arrives at five wall-seconds rather than eight —
    /// inside the ten seconds the ticket says a stranger will give the page.
    /// The later kicks have no such deadline and get the full margin.
    double first_kick_s = 2.5;
    double period_s = 4.0;

    /// Kick magnitude, as a velocity impulse.  Chosen at the magnitude of the
    /// nudge `test_default_weights_reject_a_nudge` already runs — 0.43 m/s —
    /// and measured here at 43 to 85 mm of excursion depending on direction,
    /// against a plate of 300 mm radius.  Unmissable in the 3D view, and
    /// nowhere near the edge.
    double speed = 0.42;  ///< [m/s]

    /// The advance in direction from one kick to the next: the golden angle,
    /// so successive kicks point somewhere new for a very long time.  A demo
    /// that alternates between two directions is an animation loop, and a
    /// visitor spots the repeat in about fifteen seconds.
    ///
    /// Kick 0 goes along +x.  There is no first-direction field, because the
    /// plate has no preferred axis to set one against — the three legs sit at
    /// 120 degrees and the kick sequence has to cover the circle either way.
    double turn_rad = 137.5077640500378 * M_PI / 180.0;
};

/// The ball state the demo opens on: displaced, at rest.
Eigen::Vector4d attractStart(const AttractSchedule& s);

/// How many kicks are due by time `t`.
///
/// A count, not an "is one due now" predicate: the caller keeps the number it
/// has applied and compares, so no arithmetic on either side can deliver the
/// same kick twice or let a paused simulation lose one.
int attractKicksBy(const AttractSchedule& s, double t);

/// Kick `n` (0-based), as a delta on the ball state — velocity only.
///
/// A kick ADDS to whatever the ball is already doing.  That is what makes it a
/// disturbance rather than a reset: a position teleport reads as a glitch in
/// the render, an impulse reads as something having hit the ball.
Eigen::Vector4d attractKick(const AttractSchedule& s, int n);

}  // namespace caliburn
