// tests/test_attract_mode.cpp
//
// Attract mode is the demo running itself for a visitor who has not arrived
// yet, and every one of its acceptance criteria is phrased as something a
// person would see: "motion is apparent", "a disturbance is visible", "the
// recovery is visible".  None of that can be checked by looking at a canvas
// once and calling it done — the tuning, the plate geometry and the rolling
// model all move underneath it.
//
// So the schedule is a pure function of elapsed time, and this file runs it
// against the same nonlinear plate the browser does.  "Visible" becomes
// millimetres, "recovers" becomes a settling time measured before the next
// kick lands, and "the visitor never sees a broken demo" becomes: the ball
// does not leave the plate in a minute of unattended running.
#include "attract_mode.h"

#include "auto_balance.h"
#include "ball_contact.h"
#include "ball_sim.h"
#include "cascade_fixture.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>
#include <vector>

using namespace caliburn;

namespace {

constexpr double kBallRadius = kFixtureBallRadius;

// "Home again", for a plant whose rolling friction makes the centre a REGION
// rather than a point.  A state feedback has no integral term, so it parks
// wherever the tilt it is asking for falls inside the dead band and then
// creeps: measured 4.8 mm just before the next kick lands, and down to about
// 1 mm if left for twenty seconds, which is what
// `test_default_weights_balance_the_ball` reads instead.
//
// The bar is twice the worst measurement rather than the measurement itself.
// Pinned at 4.8 mm this test would be a transcript of one LQR solve; at 10 mm
// it still rejects every failure that matters — a kick not recovered, a ball
// left drifting, a tuning that stalls it half the plate out — and 10 mm on a
// 300 mm radius is a handful of pixels, which is the sense in which the ticket
// asks for the recovery to be visible.
constexpr double kHome = 0.010;  // [m]

// ---------------------------------------------------------------------------
// The schedule itself
// ---------------------------------------------------------------------------

// The opening is a displacement at rest, not a kick: at t = 0 there is nothing
// on screen yet for a kick to be a change *from*, and a ball that starts
// already off centre is moving on the first frame that draws.
void test_the_demo_opens_displaced_and_at_rest() {
    const AttractSchedule s;
    const Eigen::Vector4d b = attractStart(s);

    ASSERT_NEAR(b(0), s.start_x, 1e-15);
    ASSERT_NEAR(b(1), s.start_y, 1e-15);
    ASSERT_NEAR(b(2), 0.0, 1e-15);
    ASSERT_NEAR(b(3), 0.0, 1e-15);
    // Far enough out to read as displaced rather than as centred.
    ASSERT_TRUE(std::hypot(s.start_x, s.start_y) > 0.05);
}

// Counted, never triggered on an interval test: the caller keeps the number it
// has actually delivered and applies the difference, so a long frame or a
// paused simulation cannot silently swallow a kick.
void test_kicks_are_counted_from_the_start_and_the_period() {
    AttractSchedule s;
    s.first_kick_s = 4.0;
    s.period_s = 4.0;

    ASSERT_EQ(attractKicksBy(s, 0.0), 0);
    ASSERT_EQ(attractKicksBy(s, 3.999), 0);
    ASSERT_EQ(attractKicksBy(s, 4.0), 1);
    ASSERT_EQ(attractKicksBy(s, 7.999), 1);
    ASSERT_EQ(attractKicksBy(s, 8.0), 2);
    ASSERT_EQ(attractKicksBy(s, 20.0), 5);

    // A frame long enough to span two kicks reports both, so the caller's
    // catch-up loop delivers both rather than losing one.
    ASSERT_EQ(attractKicksBy(s, 12.5) - attractKicksBy(s, 4.5), 2);

    // A degenerate period is one kick, not an infinity of them.
    s.period_s = 0.0;
    ASSERT_EQ(attractKicksBy(s, 1000.0), 1);
}

// A kick adds to whatever the ball is already doing.  That is what makes it a
// disturbance rather than a reset — a teleport reads as a glitch, a velocity
// impulse reads as something having hit the ball.
void test_a_kick_is_velocity_only() {
    const AttractSchedule s;
    for (int n = 0; n < 8; ++n) {
        const Eigen::Vector4d k = attractKick(s, n);
        ASSERT_NEAR(k(0), 0.0, 1e-15);
        ASSERT_NEAR(k(1), 0.0, 1e-15);
        ASSERT_NEAR(std::hypot(k(2), k(3)), s.speed, 1e-12);
    }
}

// The golden angle, and the reason for it: a demo whose kicks alternate
// between two directions is an animation loop, and a visitor spots the repeat
// in about fifteen seconds.  Twelve kicks is nearly a minute of running, and
// no two of them point within ten degrees of each other.
void test_successive_kicks_point_somewhere_new() {
    const AttractSchedule s;
    std::vector<double> dir;
    for (int n = 0; n < 12; ++n) {
        const Eigen::Vector4d k = attractKick(s, n);
        dir.push_back(std::atan2(k(3), k(2)));
    }
    for (std::size_t i = 0; i < dir.size(); ++i) {
        for (std::size_t j = i + 1; j < dir.size(); ++j) {
            double d = std::fmod(std::abs(dir[i] - dir[j]), 2.0 * M_PI);
            if (d > M_PI) d = 2.0 * M_PI - d;
            ASSERT_TRUE(d > 10.0 * M_PI / 180.0);
        }
    }
}

// ---------------------------------------------------------------------------
// The demo, run against the nonlinear plate
// ---------------------------------------------------------------------------

struct Sample {
    double t;
    double radius;  // [m] from the centre, which is also the setpoint
};

struct DemoRun {
    std::vector<Sample> trace;
    bool left_plate = false;
    int kicks = 0;
    int airborne_frames = 0;   ///< the ball can leave upward now, see #23
    double peak_height = 0.0;  ///< [m] above the plate surface
};

// PlateView::step's calls, in PlateView::step's order: legCommand ->
// stepServosOnPlate -> solve_pose -> attract kick -> stepBall.  The kick
// lands after the controller has read the ball, which is the causal order a
// disturbance actually has — the loop sees it on the next frame, not the one
// it happened on.
DemoRun runDemo(const ModelEntry& e, const Eigen::MatrixXd& K,
                const AttractSchedule& s, double duration) {
    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);

    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;

    const RollingBallDynamics dynamics(
        kPlateBall, PlateParams{tp.R_table, cascadeGravity(e.params)});

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    BallState ball;
    ball.rolling = attractStart(s);
    TablePose pose = tk.home_pose(d.home_leg_rad);
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;

    const double dt = 1.0 / 60.0;
    const int steps = static_cast<int>(duration / dt);
    const double g = cascadeGravity(e.params);

    DemoRun r;
    double t = 0.0;
    r.trace.push_back({0.0, std::hypot(ball.rolling(0), ball.rolling(1))});
    for (int k = 0; k < steps; ++k) {
        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        const LegCommand c = legCommand(tk, d, alpha, seen, 0.0, 0.0);
        t += dt;
        std::array<double, 3> adot{};
        for (int i = 0; i < 3; ++i)
            adot[i] = (c.alpha_rad[i] - alpha[i]) / d.servo_tau;
        alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);

        const FKResult fk = tk.solve_pose(alpha, pose);
        if (fk.converged) pose = fk.pose;
        pm_prev = pm;
        pm = plateMotion(tk, pose, alpha, adot);

        for (const int due = attractKicksBy(s, t); r.kicks < due; ++r.kicks)
            if (!ball.airborne) ball.rolling += attractKick(s, r.kicks);

        ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                               kBallRadius, g, dt);
        if (ball.airborne) ++r.airborne_frames;

        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
        r.peak_height = std::max(r.peak_height, now(2) - kBallRadius);
        r.trace.push_back({t, std::hypot(now(0), now(1))});
        if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                         tp.R_table, kBallRadius)) {
            r.left_plate = true;
            return r;
        }
    }
    return r;
}

double radiusAt(const DemoRun& r, double t) {
    double best = 0.0;
    for (const Sample& s : r.trace) {
        if (s.t > t) break;
        best = s.radius;
    }
    return best;
}

double peakBetween(const DemoRun& r, double t0, double t1) {
    double peak = 0.0;
    for (const Sample& s : r.trace)
        if (s.t >= t0 && s.t <= t1) peak = std::max(peak, s.radius);
    return peak;
}

// "Motion is apparent within the first few seconds."  As a number: the ball
// starts 72 mm out and has crossed most of that within one second — measured
// 8.9 mm, against a bar set at half the opening distance.  A still image
// cannot do that, and neither can a demo that waits for a click.
void test_the_opening_is_moving_within_a_second() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const AttractSchedule s;

    const DemoRun r = runDemo(e, defaultGain(e), s, 4.0);
    ASSERT_TRUE(!r.left_plate);

    const double start = radiusAt(r, 0.0);
    ASSERT_TRUE(start > 0.05);
    ASSERT_TRUE(radiusAt(r, 1.0) < 0.5 * start);
}

// "Controller is already active and stabilising at load", and the recovery
// finishes before the first kick lands — a kick arriving mid-recovery reads as
// noise, not as a loop rejecting a disturbance.
void test_the_opening_displacement_is_recovered_before_the_first_kick() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const AttractSchedule s;

    const DemoRun r = runDemo(e, defaultGain(e), s, s.first_kick_s);
    ASSERT_TRUE(!r.left_plate);
    ASSERT_EQ(r.kicks, 0);
    ASSERT_TRUE(radiusAt(r, s.first_kick_s - 0.1) < kHome);
}

// "A disturbance is visible, and recovery from it is visible."  Every kick in
// a minute of running, individually: it throws the ball at least 18 mm out —
// measured 24 to 30 mm — and the ball is home again before the next one
// arrives.
//
// The bar came down with the kick when #22 showed the old magnitude was
// outside what this loop can actually reject.  30 mm is a tenth of the plate's
// radius, which is plainly visible motion on a canvas a few hundred pixels
// wide — and, unlike the 85 mm this used to claim, it is a disturbance the
// controller survives from every direction rather than from most of them.
void test_every_kick_is_visible_and_recovered_before_the_next() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const AttractSchedule s;

    const DemoRun r = runDemo(e, defaultGain(e), s, 60.0);
    ASSERT_TRUE(!r.left_plate);
    ASSERT_TRUE(r.kicks >= 10);

    for (int n = 0; n < r.kicks; ++n) {
        const double t_kick = s.first_kick_s + n * s.period_s;
        const double t_next = t_kick + s.period_s;
        ASSERT_TRUE(peakBetween(r, t_kick, t_next) > 0.018);
        ASSERT_TRUE(radiusAt(r, t_next - 0.1) < kHome);
    }
}

// A minute unattended, and the ball is still on the plate.  The failure this
// pins is the one that would be worst in front of a stranger: an auto-reset
// blink, or a demo that quietly stops having a ball in it.
void test_a_minute_unattended_never_loses_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), AttractSchedule{}, 60.0);
    ASSERT_TRUE(!r.left_plate);
    // And never even came close: the ball leaves at R_table - r_ball.
    const TableParams tp = cascadeMechanism(e.params);
    ASSERT_TRUE(peakBetween(r, 0.0, 60.0) < 0.5 * (tp.R_table - kBallRadius));
}

// Issue #22, as the two checks that would have caught it.
//
// The old one-minute run passed while the demo was broken, because it saw
// only the fifteen directions the golden angle happens to produce in sixty
// seconds and the first failure was at 59.3 s.  These sweep the disturbance
// instead of sampling it, and run long enough for a permanent failure to have
// somewhere to show itself.

// One kick, swept over every direction, from the state the schedule actually
// kicks from — the loop settled, the ball parked a few millimetres off centre
// in the friction dead band, the legs slightly off home.  Kicking a pristine
// plate at exact centre would be an easier test than the thing it stands in
// for: #22's real failure had the ball at +3.4, -3.3 mm with the legs already
// displaced.
double sweepOneDirection(const ModelEntry& e, const Eigen::MatrixXd& K,
                         const AttractSchedule& s, double theta) {
    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;
    const RollingBallDynamics dynamics(
        kPlateBall, PlateParams{tp.R_table, cascadeGravity(e.params)});

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    TablePose pose = tk.home_pose(d.home_leg_rad);
    BallState ball;
    ball.rolling = attractStart(s);
    PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;
    const double g = cascadeGravity(e.params);

    const double dt = 1.0 / 60.0;
    double peak = 0.0;
    // Settle out of the opening displacement first, exactly as the demo does,
    // then kick from wherever that left everything.
    const int settle = static_cast<int>(s.first_kick_s / dt);
    const int total = settle + static_cast<int>(12.0 / dt);
    for (int k = 0; k < total; ++k) {
        if (k == settle && !ball.airborne) {
            ball.rolling(2) += s.speed * std::cos(theta);
            ball.rolling(3) += s.speed * std::sin(theta);
        }
        const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
        Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
        if (ball.airborne) {
            const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
            seen << land(0), land(1), 0.0, 0.0;
        }
        const LegCommand c = legCommand(tk, d, alpha, seen, 0.0, 0.0);
        std::array<double, 3> adot{};
        for (int i = 0; i < 3; ++i)
            adot[i] = (c.alpha_rad[i] - alpha[i]) / d.servo_tau;
        alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);
        const FKResult fk = tk.solve_pose(alpha, pose);
        // The plate always has an assembly now.  A stronger statement than
        // "the ball stayed on", and the one the fix actually makes.
        ASSERT_TRUE(fk.converged);
        pose = fk.pose;
        pm_prev = pm;
        pm = plateMotion(tk, pose, alpha, adot);
        ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                               kBallRadius, g, dt);
        const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
        if (k >= settle) peak = std::max(peak, std::hypot(now(0), now(1)));
        ASSERT_TRUE(ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                                tp.R_table, kBallRadius));
    }
    const Eigen::Matrix<double, 6, 1> end = plateFrame(ball, pm, kBallRadius);
    ASSERT_TRUE(std::hypot(end(0), end(1)) < kHome);   // and came home
    return peak;
}

// Every direction, not the handful the schedule reaches.  A kick the loop
// cannot reject is a kick the demo must not deliver, wherever it points.
void test_the_kick_is_rejected_from_every_direction() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);
    const Eigen::MatrixXd K = defaultGain(e);

    double worst_peak = 0.0;
    for (int i = 0; i < 720; ++i)
        worst_peak = std::max(worst_peak,
                              sweepOneDirection(e, K, AttractSchedule{},
                                                i * M_PI / 360.0));

    // Visible, and nowhere near the edge: 30 mm against the 280 mm the plate
    // has.  Measured; the assertion is loose around it on purpose.
    ASSERT_TRUE(worst_peak > 0.018);
    // Wider than it was before #23: a ball that hops carries its speed without
    // rolling friction, so it lands further out.  Still a quarter of what the
    // plate has.
    ASSERT_TRUE(worst_peak < 0.09);
}

// "A tuning that genuinely saturates the legs still recovers."
//
// Q on ball position at 200 saturates the legs and is clipped against the
// workspace from every direction, and still brings the ball back from all 180
// tested without losing it — measured 0 of 720 at this weight.  That is the
// workspace retreat doing its job under a tuning harder than the shipped one.
void test_a_saturating_tuning_still_recovers() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q(7);
    q << 1, 1, 1, 200, 200, 2, 2;
    const Eigen::MatrixXd K = gainFor(e, q, defaultLqrInputWeights(3));

    for (int i = 0; i < 180; ++i)
        sweepOneDirection(e, K, AttractSchedule{}, i * M_PI / 90.0);
}

// And the fact that only exists now the ball can leave: an over-aggressive
// tuning does not merely overshoot, it THROWS THE BALL OFF THE PLATE.
//
// Q at 2000 slams the legs hard enough that the plate is moving out from under
// the ball rather than tilting under it: measured, 54 of 720 kick directions
// end with the ball gone, and one launch reached 689 mm of altitude.  Before
// #23 the same tuning looked merely fast — 0.35 s to settle — because a ball
// glued to the plate cannot be thrown off it.
//
// This is the honest ceiling on #19's aggressive preset, and a far better
// demonstration of over-aggressive control than a settling time.
void test_an_over_aggressive_tuning_throws_the_ball_off() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q(7);
    q << 1, 1, 1, 2000, 2000, 2, 2;
    const Eigen::MatrixXd K = gainFor(e, q, defaultLqrInputWeights(3));

    const TableParams tp = cascadeMechanism(e.params);
    const TableKinematics tk(tp);
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = cascadeHomeLegAngle(e.params);
    d.servo_tau = cascadeServoTau(e.params);
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;
    const double g = cascadeGravity(e.params);
    const RollingBallDynamics dynamics(kPlateBall, PlateParams{tp.R_table, g});
    const AttractSchedule s;
    const double dt = 1.0 / 60.0;

    int lost = 0;
    double peak_hop = 0.0;
    for (int i = 0; i < 90; ++i) {
        const double theta = i * M_PI / 45.0;
        std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
        TablePose pose = tk.home_pose(d.home_leg_rad);
        BallState ball;
        ball.rolling = attractStart(s);
        PlateMotion pm = plateMotion(tk, pose, alpha, {0.0, 0.0, 0.0}), pm_prev = pm;

        const int settle = static_cast<int>(s.first_kick_s / dt);
        for (int k = 0; k < settle + static_cast<int>(6.0 / dt); ++k) {
            if (k == settle && !ball.airborne) {
                ball.rolling(2) += s.speed * std::cos(theta);
                ball.rolling(3) += s.speed * std::sin(theta);
            }
            const Eigen::Matrix<double, 6, 1> bp = plateFrame(ball, pm, kBallRadius);
            Eigen::Vector4d seen(bp(0), bp(1), bp(3), bp(4));
            if (ball.airborne) {
                const Eigen::Vector2d land = predictedLanding(bp, kBallRadius, g);
                seen << land(0), land(1), 0.0, 0.0;
            }
            const LegCommand c = legCommand(tk, d, alpha, seen, 0.0, 0.0);
            std::array<double, 3> adot{};
            for (int j = 0; j < 3; ++j)
                adot[j] = (c.alpha_rad[j] - alpha[j]) / d.servo_tau;
            alpha = stepServosOnPlate(tk, alpha, c.alpha_rad, d.servo_tau, dt);
            const FKResult fk = tk.solve_pose(alpha, pose);
            if (fk.converged) pose = fk.pose;
            pm_prev = pm;
            pm = plateMotion(tk, pose, alpha, adot);
            ball = stepBallContact(dynamics, ball, pm, pm_prev, pose,
                                   kBallRadius, g, dt);
            const Eigen::Matrix<double, 6, 1> now = plateFrame(ball, pm, kBallRadius);
            peak_hop = std::max(peak_hop, now(2) - kBallRadius);
            if (!ballOnPlate(Eigen::Vector4d(now(0), now(1), now(3), now(4)),
                             tp.R_table, kBallRadius)) { ++lost; break; }
        }
    }
    ASSERT_TRUE(lost > 0);            // it really does lose the ball
    ASSERT_TRUE(peak_hop > 0.05);     // and really does launch it, metres not mm
}

// Ten minutes unattended.  The failure this pins is the one that shipped: the
// ball left at 59.3 s and then left on every kick after it, for as long as the
// page stayed open, because the plate had folded onto an assembly it could
// never climb out of.
void test_ten_minutes_unattended_never_loses_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const DemoRun r = runDemo(e, defaultGain(e), AttractSchedule{}, 600.0);
    ASSERT_TRUE(!r.left_plate);
    ASSERT_TRUE(r.kicks >= 140);

    const TableParams tp = cascadeMechanism(e.params);
    ASSERT_TRUE(peakBetween(r, 0.0, 600.0) < 0.5 * (tp.R_table - kBallRadius));

    // Still balancing at the end, not merely still on the plate.
    ASSERT_TRUE(radiusAt(r, 600.0) < kHome);
}

}  // namespace

int main() {
    test_the_demo_opens_displaced_and_at_rest();
    test_kicks_are_counted_from_the_start_and_the_period();
    test_a_kick_is_velocity_only();
    test_successive_kicks_point_somewhere_new();
    test_the_opening_is_moving_within_a_second();
    test_the_opening_displacement_is_recovered_before_the_first_kick();
    test_every_kick_is_visible_and_recovered_before_the_next();
    test_a_minute_unattended_never_loses_the_ball();
    test_the_kick_is_rejected_from_every_direction();
    test_a_saturating_tuning_still_recovers();
    test_an_over_aggressive_tuning_throws_the_ball_off();
    test_ten_minutes_unattended_never_loses_the_ball();
    std::printf("test_attract_mode: all passed\n");
    return 0;
}
