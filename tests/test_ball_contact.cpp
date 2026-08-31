// tests/test_ball_contact.cpp
//
// Whether a ball stays on the plate is a question about the surface's
// ACCELERATION, and the rolling model cannot ask it: `RollingBallDynamics`
// carries [x, y, vx, vy] in the plate frame, has no normal force, and glues
// the ball down forever.  The plate heaves hard — `z_c = 0.30 sin(alpha)` for
// this geometry, so a 20 degree leg swing moves the table 74 mm vertically and
// the loop does that in about a tenth of a second — so the glue is doing real
// work.
//
// These tests build the normal force up from cases where the answer is known
// by hand, then point it at the shipped tuning.  See issue #23.
#include "ball_contact.h"

#include "cascade_fixture.h"
#include "test_helpers.h"

#include <cmath>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kG = 9.81;
constexpr double kR = kFixtureBallRadius;

TableParams plate() {
    return cascadeMechanism(cascadeModel(getBuiltinModels()).params);
}

/// A plate that is not moving at all, at a given tilt.
PlateMotion still(const TableKinematics& tk, double phi, double theta, double z_c) {
    PlateMotion m;
    m.R = tk.table_rotation(phi, theta);
    m.c = Eigen::Vector3d(0, 0, z_c);
    return m;
}

// A ball sitting on a level, motionless plate is held by exactly its weight.
// The whole calculation reduces to this, and if it does not, nothing further
// down is worth reading.
void test_a_still_level_plate_holds_the_ball_at_one_g() {
    const TableKinematics tk(plate());
    const PlateMotion m = still(tk, 0.0, 0.0, 0.2121);
    const double N = normalAccel(m, m, 1.0 / 60.0,
                                 Eigen::Vector3d(0.0, 0.0, kR),
                                 Eigen::Vector2d::Zero(), kG);
    ASSERT_NEAR(N, kG, 1e-12);
}

// Tilt it and the normal takes only gravity's share along the normal, which is
// the cosine.  The rest is what rolls the ball, and belongs to the roller.
void test_a_tilted_still_plate_holds_the_cosine_share() {
    const TableKinematics tk(plate());
    for (double deg : {5.0, 15.0, 30.0}) {
        const PlateMotion m = still(tk, deg * kDeg, 0.0, 0.2121);
        const double N = normalAccel(m, m, 1.0 / 60.0,
                                     Eigen::Vector3d(0.0, 0.0, kR),
                                     Eigen::Vector2d::Zero(), kG);
        ASSERT_NEAR(N, kG * std::cos(deg * kDeg), 1e-12);
    }
}

// The case the whole ticket is about: drop the plate away from under the ball.
// At exactly one g of downward heave the ball is weightless — the plate is
// falling with it and touching it with nothing.  Below that, contact is over.
void test_heave_cancels_the_normal_force_at_one_g() {
    const TableKinematics tk(plate());
    const double dt = 1.0 / 60.0;

    auto heaving = [&](double zc_ddot) {
        PlateMotion now = still(tk, 0, 0, 0.2121), prev = now;
        // c_ddot comes from differencing c_dot, so stage the two frames.
        now.c_dot = Eigen::Vector3d(0, 0, zc_ddot * dt);
        return normalAccel(now, prev, dt, Eigen::Vector3d(0, 0, kR),
                           Eigen::Vector2d::Zero(), kG);
    };

    ASSERT_NEAR(heaving(0.0), kG, 1e-12);
    ASSERT_NEAR(heaving(-kG), 0.0, 1e-9);        // weightless, exactly
    ASSERT_TRUE(heaving(-1.5 * kG) < 0.0);       // separated
    ASSERT_NEAR(heaving(kG), 2.0 * kG, 1e-9);    // and pushed twice as hard
}

// Steady rotation is centripetal only, and about a horizontal axis through the
// table centre that pulls straight down on the ball's own radius — the same
// for a ball anywhere on the plate, because a sideways offset adds nothing to
// the vertical component.  Worth pinning precisely: it is the term most easily
// confused with the next one.
void test_steady_rotation_pulls_on_the_balls_own_radius() {
    const TableKinematics tk(plate());
    const double dt = 1.0 / 60.0;

    PlateMotion now = still(tk, 0, 0, 0.2121), prev = now;
    now.omega = Eigen::Vector3d(2.0, 0.0, 0.0);
    prev.omega = now.omega;

    const double centred = normalAccel(now, prev, dt, Eigen::Vector3d(0, 0, kR),
                                       Eigen::Vector2d::Zero(), kG);
    const double off = normalAccel(now, prev, dt, Eigen::Vector3d(0, 0.10, kR),
                                   Eigen::Vector2d::Zero(), kG);
    ASSERT_NEAR(centred, kG - 2.0 * 2.0 * kR, 1e-9);
    ASSERT_NEAR(off, centred, 1e-12);
}

// Angular ACCELERATION is the one that cares where the ball is: the plate
// tipping about its centre lifts one side and drops the other, so two balls on
// opposite sides are held by different forces at the same instant.  This is the
// term that takes the ball off the plate when the loop slams the legs over.
void test_angular_acceleration_lifts_one_side_and_drops_the_other() {
    const TableKinematics tk(plate());
    const double dt = 1.0 / 60.0;
    const double alpha_dd = 40.0;   // [rad/s^2] about x

    PlateMotion now = still(tk, 0, 0, 0.2121), prev = now;
    now.omega = Eigen::Vector3d(alpha_dd * dt, 0.0, 0.0);   // prev.omega is zero

    const double centred = normalAccel(now, prev, dt, Eigen::Vector3d(0, 0, kR),
                                       Eigen::Vector2d::Zero(), kG);
    const double plus_y = normalAccel(now, prev, dt, Eigen::Vector3d(0, 0.10, kR),
                                      Eigen::Vector2d::Zero(), kG);
    const double minus_y = normalAccel(now, prev, dt, Eigen::Vector3d(0, -0.10, kR),
                                       Eigen::Vector2d::Zero(), kG);

    // omega_dot x r has vertical component alpha_dd * y, so the two sides
    // straddle the centred case by that much, symmetrically.
    ASSERT_NEAR(plus_y - centred, alpha_dd * 0.10, 1e-6);
    ASSERT_NEAR(centred - minus_y, alpha_dd * 0.10, 1e-6);
    // 40 rad/s^2 at 100 mm is 4 m/s^2 — nearly half a g, from tipping alone.
    ASSERT_TRUE(minus_y < kG - 3.0);
}

// Coriolis: a ball ROLLING on a tilting plate is held by a different force
// than one sitting still on it, at the same instant and the same place.
void test_a_rolling_ball_is_held_differently_from_a_still_one() {
    const TableKinematics tk(plate());
    const double dt = 1.0 / 60.0;

    PlateMotion now = still(tk, 0, 0, 0.2121), prev = now;
    now.omega = Eigen::Vector3d(1.5, 0.0, 0.0);   // tilting about x
    prev.omega = now.omega;

    const Eigen::Vector3d s(0.05, 0.0, kR);
    const double at_rest = normalAccel(now, prev, dt, s, Eigen::Vector2d(0, 0), kG);
    const double rolling = normalAccel(now, prev, dt, s, Eigen::Vector2d(0, 0.4), kG);
    ASSERT_TRUE(std::abs(rolling - at_rest) > 0.1);
    // Rolling along +y under a +x rotation is carried downward, so it presses
    // less hard.  Reverse the roll and it presses harder by the same amount.
    const double other_way = normalAccel(now, prev, dt, s, Eigen::Vector2d(0, -0.4), kG);
    ASSERT_NEAR(rolling + other_way, 2.0 * at_rest, 1e-9);
}

// ---------------------------------------------------------------------------
// Leaving, flying, and landing
// ---------------------------------------------------------------------------

RollingBallDynamics roller() {
    const TableParams tp = plate();
    return RollingBallDynamics(kPlateBall, PlateParams{tp.R_table, kG});
}

// A ball on a plate that is behaving itself never leaves, and the rolling path
// is bit-for-bit the one `stepBall` gives — the contact layer adds a question,
// not a different answer.
//
// `stepBall` has to be told the same normal force the contact layer computed,
// which since #23 is the plate's cosine share of gravity and not `g`.  On a
// motionless plate `normalAccel` reduces to exactly `quasiStaticNormalAccel`,
// so this stays an equality rather than becoming a tolerance.  Handing the bare
// roller `g` here instead would compare the contact model against a ball on a
// LEVEL plate, which at 2 deg of tilt disagrees in the fourth decimal — the
// divergence this assertion exists to catch.
void test_a_settled_plate_never_lets_go() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{2.0 * kDeg, -1.0 * kDeg, 0.2121};
    const PlateMotion m = still(tk, pose.phi, pose.theta, pose.z_c);

    BallState b;
    b.rolling << 0.05, -0.02, 0.0, 0.0;
    Eigen::Vector4d bare = b.rolling;

    for (int k = 0; k < 240; ++k) {
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, 1.0 / 60.0);
        bare = stepBall(dyn, bare, pose, quasiStaticNormalAccel(m, kG), 1.0 / 60.0);
        ASSERT_TRUE(!b.airborne);
    }
    for (int i = 0; i < 4; ++i) ASSERT_NEAR(b.rolling(i), bare(i), 1e-15);
}

// Drop the plate out from under it and the ball goes ballistic — and the arc
// it follows is the one gravity alone gives, checked against the closed form
// rather than against itself.
//
// The plate has to actually recede, not merely be given a downward velocity:
// a ball in free fall drops at g, so a plate that wants to leave it behind has
// to drop faster.  At 3g it does, and the gap opens.
void test_a_dropped_plate_launches_the_ball_on_a_parabola() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;
    const double a_plate = -3.0 * kG;

    double zc = 0.2121, vzc = 0.0;
    PlateMotion prev = still(tk, 0, 0, zc);

    BallState b;
    b.rolling << 0.0, 0.0, 0.30, 0.0;              // rolling along +x

    Eigen::Vector3d p0, v0;
    int launched_at = -1;
    for (int k = 0; k < 40; ++k) {
        vzc += a_plate * dt;
        zc += vzc * dt;
        const TablePose pose{0.0, 0.0, zc};
        PlateMotion now = still(tk, 0, 0, zc);
        now.c_dot = Eigen::Vector3d(0, 0, vzc);

        const bool was_rolling = !b.airborne;
        b = stepBallContact(dyn, b, now, prev, pose, kR, kG, dt);
        prev = now;

        if (was_rolling && b.airborne) {
            launched_at = k;
            // It left carrying its rolling speed, and the plate's motion.
            ASSERT_NEAR(b.flight_v(0), 0.30, 1e-12);
            ASSERT_TRUE(b.flight_v(2) < 0.0);
            // Rewind the one step it took after leaving, to get the launch.
            v0 = b.flight_v - Eigen::Vector3d(0, 0, -kG) * dt;
            p0 = b.flight_p - v0 * dt - 0.5 * Eigen::Vector3d(0, 0, -kG) * dt * dt;
        }
    }
    ASSERT_TRUE(launched_at == 0);      // the very first frame, at 3g
    ASSERT_TRUE(b.airborne);            // and the plate never caught it again

    // s = ut + at^2/2, in the world, from the launch.
    const double t = (40 - launched_at) * dt;
    ASSERT_NEAR(b.flight_p(0), p0(0) + v0(0) * t, 1e-12);
    ASSERT_NEAR(b.flight_p(2), p0(2) + v0(2) * t - 0.5 * kG * t * t, 1e-12);
}

// And it comes back.  Landing is inelastic: the approach speed along the
// normal is absorbed, the sideways carry is kept, and the ball is rolling
// again from where it touched down rather than from where it took off.
void test_the_ball_lands_and_rolls_on() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const TablePose pose{0.0, 0.0, 0.2121};
    const double dt = 1.0 / 60.0;
    const PlateMotion m = still(tk, 0, 0, pose.z_c);

    BallState b;
    b.airborne = true;
    b.flight_p = Eigen::Vector3d(0.04, 0.0, pose.z_c + kR + 0.05);  // 50 mm up
    b.flight_v = Eigen::Vector3d(0.20, 0.0, 0.0);                   // drifting +x

    int frames = 0;
    while (b.airborne && frames < 600) {
        b = stepBallContact(dyn, b, m, m, pose, kR, kG, dt);
        ++frames;
    }
    ASSERT_TRUE(!b.airborne);
    // Free fall from 50 mm is about 0.10 s; one frame either side is fine.
    ASSERT_NEAR(frames * dt, std::sqrt(2.0 * 0.05 / kG), 0.02);
    // Carried downrange while it fell, and still moving that way.
    ASSERT_TRUE(b.rolling(0) > 0.04);
    ASSERT_NEAR(b.rolling(2), 0.20, 1e-9);
    // The drop is gone, not turned into a bounce.
    ASSERT_NEAR(b.rolling(3), 0.0, 1e-12);
}

// The plate-frame view is the same physical ball in both phases.  Converting a
// rolling ball out to the world and straight back must return exactly what it
// started as — including the velocity, which is the part that can go wrong:
// the plate frame is moving, so the world velocity of a ball rolling on a
// tilting plate is not its rolling velocity, and the way back has to remove
// again precisely what the way out added.
void test_the_two_frames_describe_the_same_ball() {
    const TableKinematics tk(plate());

    // A plate that is doing everything at once: tilted, heaving, and rotating.
    PlateMotion m = still(tk, 4.0 * kDeg, -3.0 * kDeg, 0.2121);
    m.c_dot = Eigen::Vector3d(0.0, 0.0, 0.35);
    m.omega = Eigen::Vector3d(0.8, -0.5, 0.2);

    BallState rolling;
    rolling.rolling << 0.06, -0.03, 0.10, 0.05;

    BallState flying;
    flying.airborne = true;
    worldOf(rolling, m, kR, &flying.flight_p, &flying.flight_v);

    const Eigen::Matrix<double, 6, 1> a = plateFrame(rolling, m, kR);
    const Eigen::Matrix<double, 6, 1> b = plateFrame(flying, m, kR);
    for (int i = 0; i < 6; ++i) ASSERT_NEAR(b(i), a(i), 1e-12);

    // And the world velocity really is more than the rolling velocity — if it
    // were not, the round trip above would be proving nothing.
    ASSERT_TRUE((flying.flight_v - m.R * Eigen::Vector3d(0.10, 0.05, 0.0)).norm() > 0.1);
}

}  // namespace

// ---------------------------------------------------------------------------
// The trust gate
// ---------------------------------------------------------------------------
//
// A separation is a claim about the plate's VELOCITY, and those rates come
// through `-J_pose^-1 J_alpha`, which near a singularity amplifies without
// bound.  `rates_trustworthy` is what stops the contact model acting on rates
// the mechanism's own condition number says are meaningless — the guard
// CONTEXT.md gives a pull-quote to, and which nothing tested.

// A plate diving hard enough to leave the ball behind — but whose rates nobody
// believes — keeps the ball.  Refusing to act is the whole point of the guard:
// the alternative is a hop that is arithmetic rather than physics.
void test_untrustworthy_rates_never_separate_the_ball() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;

    // Falling at 3g, which for a trusted plate is a launch (see the parabola
    // test above), and the only difference here is the flag.
    PlateMotion prev = still(tk, 0, 0, 0.2121);
    PlateMotion now = still(tk, 0, 0, 0.2121 - 0.5 * 3.0 * kG * dt * dt);
    now.c_dot = Eigen::Vector3d(0, 0, -3.0 * kG * dt);

    BallState trusted;
    trusted.rolling << 0.03, 0.0, 0.0, 0.0;
    BallState doubted = trusted;

    now.rates_trustworthy = true;
    prev.rates_trustworthy = true;
    trusted = stepBallContact(dyn, trusted, now, prev, TablePose{0, 0, 0.2121},
                              kR, kG, dt);
    ASSERT_TRUE(trusted.airborne);

    now.rates_trustworthy = false;
    doubted = stepBallContact(dyn, doubted, now, prev, TablePose{0, 0, 0.2121},
                              kR, kG, dt);
    ASSERT_TRUE(!doubted.airborne);
}

// And the half of the gate that was missing: `normalAccel` reaches an
// acceleration by differencing `now` against `prev`, so the FIRST believable
// frame after a singular stretch is differenced against a garbage one.  Trust
// has to cover both ends of that difference or the guard leaks exactly the hop
// it exists to refuse.
void test_a_trustworthy_frame_after_an_untrusted_one_still_declines() {
    const TableKinematics tk(plate());
    const RollingBallDynamics dyn = roller();
    const double dt = 1.0 / 60.0;

    // `prev` carries a wild rate — the arithmetic a near-singular Jacobian
    // hands back — and `now` is a plate sitting still.  Differenced, that reads
    // as an enormous upward acceleration reversing, and the ball is reported
    // to leave a motionless plate.
    PlateMotion prev = still(tk, 0, 0, 0.2121);
    prev.c_dot = Eigen::Vector3d(0, 0, 12.0);
    prev.rates_trustworthy = false;

    PlateMotion now = still(tk, 0, 0, 0.2121);
    now.rates_trustworthy = true;

    BallState b;
    b.rolling << 0.03, 0.0, 0.0, 0.0;

    // The rates alone would say separation; the gate says no.
    const double N = normalAccel(now, prev, dt,
                                 Eigen::Vector3d(0.03, 0.0, kR),
                                 Eigen::Vector2d::Zero(), kG);
    ASSERT_TRUE(N < 0.0);

    b = stepBallContact(dyn, b, now, prev, TablePose{0, 0, 0.2121}, kR, kG, dt);
    ASSERT_TRUE(!b.airborne);
}

// The threshold is the application's own "Poor" line, not a second opinion.
void test_the_trust_threshold_is_the_condition_number_the_app_shows() {
    ASSERT_NEAR(kRatesUntrustworthyAbove, 20.0, 1e-15);
}

int main() {
    test_a_still_level_plate_holds_the_ball_at_one_g();
    test_a_tilted_still_plate_holds_the_cosine_share();
    test_heave_cancels_the_normal_force_at_one_g();
    test_steady_rotation_pulls_on_the_balls_own_radius();
    test_angular_acceleration_lifts_one_side_and_drops_the_other();
    test_a_rolling_ball_is_held_differently_from_a_still_one();
    test_a_settled_plate_never_lets_go();
    test_a_dropped_plate_launches_the_ball_on_a_parabola();
    test_the_ball_lands_and_rolls_on();
    test_the_two_frames_describe_the_same_ball();
    test_untrustworthy_rates_never_separate_the_ball();
    test_a_trustworthy_frame_after_an_untrusted_one_still_declines();
    test_the_trust_threshold_is_the_condition_number_the_app_shows();
    std::printf("test_ball_contact: all passed\n");
    return 0;
}
