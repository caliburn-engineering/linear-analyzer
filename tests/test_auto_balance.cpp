// tests/test_auto_balance.cpp
//
// Closing the loop is the one place where the two halves of this application
// have to agree numerically rather than merely coexist.  The gain comes from a
// LINEAR model; it drives the NONLINEAR simulator — a Newton-solved mechanism
// and a rolling ball, neither of which the design ever saw.
//
// So the tests come in two layers.  The first pins the wiring: which state
// goes in which slot, what the setpoint means, where the servo travel clips.
// The second runs the real thing and asserts what the ticket actually asks
// for — that the ball balances under a sensible tuning and visibly does not
// under a poor one.  A criterion phrased as "visibly" is still a measurable
// claim, and measuring it here is what keeps it from being re-judged by eye
// every time the model moves.
#include "auto_balance.h"

#include "analysis/lqr.h"
#include "analysis/model_library.h"
#include "ball_sim.h"
#include "rolling_dynamics.h"
#include "table_kinematics.h"
#include "test_helpers.h"

#include <cmath>
#include <string>

using namespace caliburn;

namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr double kHome = 45.0 * kDeg;
constexpr double kBallRadius = 0.02;
constexpr double kBallMass = 0.05;
constexpr double kBallFriction = 0.01;

AutoBalanceDesign designWith(const Eigen::MatrixXd& K) {
    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = kHome;
    d.servo_tau = 0.05;
    return d;
}

// ---------------------------------------------------------------------------
// Layer 1: the wiring
// ---------------------------------------------------------------------------

void test_state_vector_layout() {
    const std::array<double, 3> alpha = {kHome + 0.1, kHome - 0.2, kHome};
    const Eigen::Vector4d ball(0.03, -0.04, 0.5, -0.6);

    const Eigen::VectorXd x = cascadeState(alpha, kHome, ball);

    ASSERT_EQ((int)x.size(), 7);
    ASSERT_NEAR(x(0), 0.1, 1e-12);
    ASSERT_NEAR(x(1), -0.2, 1e-12);
    ASSERT_NEAR(x(2), 0.0, 1e-12);
    ASSERT_NEAR(x(3), 0.03, 1e-12);   // ball x
    ASSERT_NEAR(x(4), -0.04, 1e-12);  // ball y
    ASSERT_NEAR(x(5), 0.5, 1e-12);    // ball x'
    ASSERT_NEAR(x(6), -0.6, 1e-12);   // ball y'
}

// The deviation is from the DESIGN's home angle, not from any fixed 45 degrees.
// Get this wrong and the loop still converges — to a tilted plate.
void test_state_deviation_follows_the_design_home() {
    const std::array<double, 3> alpha = {60.0 * kDeg, 60.0 * kDeg, 60.0 * kDeg};
    const Eigen::VectorXd x = cascadeState(alpha, 60.0 * kDeg, Eigen::Vector4d::Zero());
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(x(i), 0.0, 1e-12);
}

void test_shape_gate() {
    AutoBalanceDesign d = designWith(Eigen::MatrixXd::Zero(3, 7));
    ASSERT_TRUE(gainFitsCascade(d));
    d.K = Eigen::MatrixXd::Zero(1, 4);
    ASSERT_TRUE(!gainFitsCascade(d));
    d.K = Eigen::MatrixXd();
    ASSERT_TRUE(!gainFitsCascade(d));
}

// A ball at rest at the setpoint, legs at home: nothing to do.  This is the
// claim that makes a position setpoint need no feedforward.
// The model panel carries geometry as `float` sliders; the plate holds
// `double` literals.  0.300f widened is 0.3000000119, so a comparison tight
// enough to be called exact refuses two mechanisms that are the same object —
// and the balance loop then never becomes available at all.
void test_mechanism_identity_survives_the_float_round_trip() {
    TableParams plate;
    plate.R_ground = 0.300;
    plate.R_table = 0.300;
    plate.L1 = 0.150;
    plate.L2 = 0.150;

    TableParams from_sliders;
    from_sliders.R_ground = static_cast<double>(0.300f);
    from_sliders.R_table = static_cast<double>(0.300f);
    from_sliders.L1 = static_cast<double>(0.150f);
    from_sliders.L2 = static_cast<double>(0.150f);
    ASSERT_TRUE(sameMechanism(plate, from_sliders));

    // A millimetre of leg, however, is a different mechanism.
    from_sliders.L1 = 0.151;
    ASSERT_TRUE(!sameMechanism(plate, from_sliders));

    // A design nobody filled in is all zeros, and matches nothing.
    ASSERT_TRUE(!sameMechanism(plate, TableParams{}));
}

void test_zero_error_commands_the_home_pose() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Random(3, 7) * 10.0;
    const AutoBalanceDesign d = designWith(K);

    const std::array<double, 3> alpha = {kHome, kHome, kHome};
    const Eigen::Vector4d ball(0.07, -0.02, 0.0, 0.0);

    const LegCommand c = legCommand(d, alpha, ball, 0.07, -0.02);
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(c.alpha_rad[i], kHome, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

void test_command_is_home_minus_k_error() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 2.0;   // leg 0 reacts to ball x
    K(1, 4) = -3.0;  // leg 1 reacts to ball y
    K(2, 0) = 0.5;   // leg 2 reacts to its own neighbour's deviation
    const AutoBalanceDesign d = designWith(K);

    const std::array<double, 3> alpha = {kHome + 0.04, kHome, kHome};
    const Eigen::Vector4d ball(0.01, 0.02, 0.0, 0.0);

    const LegCommand c = legCommand(d, alpha, ball, 0.0, 0.0);
    ASSERT_NEAR(c.alpha_rad[0], kHome - 2.0 * 0.01, 1e-12);
    ASSERT_NEAR(c.alpha_rad[1], kHome + 3.0 * 0.02, 1e-12);
    ASSERT_NEAR(c.alpha_rad[2], kHome - 0.5 * 0.04, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

// The setpoint enters as a reference STATE, so displacing both the ball and the
// setpoint together is a no-op while displacing only one is not.
void test_setpoint_moves_the_target() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 1.0;
    const AutoBalanceDesign d = designWith(K);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    const LegCommand at = legCommand(d, alpha, Eigen::Vector4d(0.05, 0, 0, 0), 0.05, 0.0);
    ASSERT_NEAR(at.alpha_rad[0], kHome, 1e-12);

    const LegCommand off = legCommand(d, alpha, Eigen::Vector4d(0.05, 0, 0, 0), 0.0, 0.0);
    ASSERT_NEAR(off.alpha_rad[0], kHome - 0.05, 1e-12);
}

void test_command_clamps_to_servo_travel() {
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(3, 7);
    K(0, 3) = 1e4;
    K(1, 3) = -1e4;
    const AutoBalanceDesign d = designWith(K);
    const std::array<double, 3> alpha = {kHome, kHome, kHome};

    const LegCommand c = legCommand(d, alpha, Eigen::Vector4d(0.05, 0, 0, 0), 0.0, 0.0);
    ASSERT_NEAR(c.alpha_rad[0], d.alpha_min_rad, 1e-12);
    ASSERT_NEAR(c.alpha_rad[1], d.alpha_max_rad, 1e-12);
    ASSERT_NEAR(c.alpha_rad[2], kHome, 1e-12);
    ASSERT_TRUE(c.saturated);
}

void test_wrong_shape_commands_the_home_pose() {
    const AutoBalanceDesign d = designWith(Eigen::MatrixXd::Ones(2, 3));
    const std::array<double, 3> alpha = {kHome + 0.2, kHome, kHome};
    const LegCommand c = legCommand(d, alpha, Eigen::Vector4d(0.1, 0.1, 1, 1), 0, 0);
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(c.alpha_rad[i], kHome, 1e-12);
    ASSERT_TRUE(!c.saturated);
}

void test_servo_lag_is_first_order() {
    const double tau = 0.05;
    const std::array<double, 3> cmd = {1.0, 1.0, 1.0};
    std::array<double, 3> a = {0.0, 0.0, 0.0};

    // One time constant, reached in one step and in many: the exact form makes
    // the answer independent of how the interval was subdivided.
    const std::array<double, 3> one_shot = stepServos(a, cmd, tau, tau);
    ASSERT_NEAR(one_shot[0], 1.0 - std::exp(-1.0), 1e-12);

    for (int k = 0; k < 100; ++k) a = stepServos(a, cmd, tau, tau / 100.0);
    ASSERT_NEAR(a[0], 1.0 - std::exp(-1.0), 1e-9);

    // No step, no motion — and no division by a zero tau.
    const std::array<double, 3> still = stepServos(a, cmd, tau, 0.0);
    ASSERT_NEAR(still[0], a[0], 1e-15);
    const std::array<double, 3> degenerate = stepServos(a, cmd, 0.0, 1.0 / 60.0);
    ASSERT_TRUE(std::isfinite(degenerate[0]));

    // Each leg keeps its own command.
    const std::array<double, 3> mixed =
        stepServos({0.0, 0.0, 0.0}, {1.0, -1.0, 0.0}, tau, tau);
    ASSERT_TRUE(mixed[0] > 0.0 && mixed[1] < 0.0);
    ASSERT_NEAR(mixed[2], 0.0, 1e-15);
}

// ---------------------------------------------------------------------------
// Layer 2: the loop, closed around the nonlinear plate
// ---------------------------------------------------------------------------

const ModelEntry& cascadeModel(const std::vector<ModelEntry>& models) {
    for (const auto& m : models)
        if (m.name.rfind("Ball-Balancer Cascade", 0) == 0) return m;
    std::fprintf(stderr, "FAIL: no cascade model in the library\n");
    std::exit(1);
}

double paramValue(const ModelEntry& e, const std::string& symbol) {
    for (const auto& p : e.params)
        if (p.symbol == symbol) return p.value;
    std::fprintf(stderr, "FAIL: no parameter '%s'\n", symbol.c_str());
    std::exit(1);
}

TableParams tableParamsOf(const ModelEntry& e) {
    TableParams tp;
    tp.R_ground  = paramValue(e, "Rg");
    tp.R_table   = paramValue(e, "Rt");
    tp.L1        = paramValue(e, "L1");
    tp.L2        = paramValue(e, "L2");
    tp.alpha_min = 10.0 * kDeg;
    tp.alpha_max = 80.0 * kDeg;
    return tp;
}

struct SimResult {
    double final_radius;   // [m] distance of the ball from its setpoint
    double peak_radius;    // [m]
    double settle_time;    // [s] after which it never again left 5 mm; -1 = never
    double final_tilt_deg; // [deg] the larger of |roll| and |pitch| at the end
    bool left_plate;
};

// The same four calls, in the same order, that PlateView::step makes:
// legCommand -> stepServos -> forward_kinematics -> stepBall.  Nothing here is
// linearised, and nothing here knows what Q and R were.
SimResult runClosedLoop(const ModelEntry& e,
                        const Eigen::MatrixXd& K,
                        const Eigen::Vector4d& ball0,
                        double x_sp,
                        double y_sp,
                        double duration) {
    const TableParams tp = tableParamsOf(e);
    const TableKinematics tk(tp);

    AutoBalanceDesign d;
    d.K = K;
    d.home_leg_rad = paramValue(e, "a0") * kDeg;
    d.servo_tau = paramValue(e, "\xcf\x84");
    d.alpha_min_rad = tp.alpha_min;
    d.alpha_max_rad = tp.alpha_max;

    const RollingBallDynamics dynamics(
        BallParams{kBallRadius, kBallMass, kBallFriction},
        PlateParams{tp.R_table, paramValue(e, "g")});

    std::array<double, 3> alpha = {d.home_leg_rad, d.home_leg_rad, d.home_leg_rad};
    Eigen::Vector4d ball = ball0;
    TablePose pose = tk.home_pose(d.home_leg_rad);

    const double dt = 1.0 / 60.0;
    const int steps = static_cast<int>(duration / dt);

    SimResult r{0.0, 0.0, -1.0, 0.0, false};
    for (int k = 0; k < steps; ++k) {
        const LegCommand c = legCommand(d, alpha, ball, x_sp, y_sp);
        alpha = stepServos(alpha, c.alpha_rad, d.servo_tau, dt);

        const FKResult fk = tk.forward_kinematics(alpha, pose);
        if (fk.converged) pose = fk.pose;

        ball = stepBall(dynamics, ball, pose, dt);

        const double radius = std::hypot(ball(0) - x_sp, ball(1) - y_sp);
        r.peak_radius = std::max(r.peak_radius, radius);
        if (radius > 0.005) r.settle_time = (k + 1) * dt;
        if (!ballOnPlate(ball, tp.R_table, kBallRadius)) {
            r.left_plate = true;
            r.final_radius = radius;
            return r;
        }
    }
    r.final_radius = std::hypot(ball(0) - x_sp, ball(1) - y_sp);
    r.final_tilt_deg = std::max(std::abs(pose.phi), std::abs(pose.theta)) / kDeg;
    if (r.settle_time >= duration - dt) r.settle_time = -1.0;
    return r;
}

Eigen::MatrixXd gainFor(const ModelEntry& e,
                        const Eigen::VectorXd& q,
                        const Eigen::VectorXd& r) {
    const LqrResult res = computeLQR(e.system,
                                     q.asDiagonal().toDenseMatrix(),
                                     r.asDiagonal().toDenseMatrix());
    ASSERT_TRUE(res.success);
    return res.K;
}

Eigen::MatrixXd defaultGain(const ModelEntry& e) {
    return gainFor(e, defaultLqrStateWeights(7), defaultLqrInputWeights(3));
}

// 60 mm out along x and 40 mm back along y, at rest: a quarter of the plate's
// radius, plainly visible in the 3D view and well outside anything a
// linearisation could excuse.
const Eigen::Vector4d kDisplaced(0.06, -0.04, 0.0, 0.0);

// The weights the application opens on.  If THIS does not balance the ball,
// nobody who switches the controller to LQR ever sees the product work.
void test_default_weights_balance_the_ball() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(e, defaultGain(e), kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// The nudge buttons in Plate Control, as a number: an impulse in velocity,
// rejected.  This is the disturbance the "visibly degrades" criterion is read
// against, so the good case has to be pinned first.
void test_default_weights_reject_a_nudge() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(
        e, defaultGain(e), Eigen::Vector4d(0.0, 0.0, 0.35, 0.25), 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.peak_radius > 0.05);      // it really was kicked
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// A ball at rest anywhere on a flat plate is an equilibrium, so a position
// setpoint is reachable with no feedforward at all.  That claim is arithmetic
// in `legCommand`; here it is the nonlinear plate that has to agree.
void test_setpoint_is_tracked() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const SimResult s = runClosedLoop(
        e, defaultGain(e), Eigen::Vector4d::Zero(), 0.08, 0.05, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time > 0.0);
    ASSERT_TRUE(s.settle_time < 5.0);
    ASSERT_TRUE(s.final_radius < 0.005);
}

// "Moving the poles toward the imaginary axis should visibly make the ball
// sluggish" — the ticket's own words, as a comparison the user can reproduce
// by dragging the ball-position weights down.
void test_looser_weights_are_visibly_sluggish() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd q_loose(7);
    q_loose << 1, 1, 1, 50, 50, 5, 5;
    const Eigen::MatrixXd K_loose = gainFor(e, q_loose, defaultLqrInputWeights(3));

    const SimResult loose = runClosedLoop(e, K_loose, kDisplaced, 0.0, 0.0, 20.0);
    const SimResult base = runClosedLoop(e, defaultGain(e), kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(loose.settle_time > 0.0);
    ASSERT_TRUE(base.settle_time > 0.0);
    // Not "a bit slower": four times slower, which is the difference between
    // watching it arrive and wondering whether it is moving.
    ASSERT_TRUE(loose.settle_time > 4.0 * base.settle_time);
}

// The third acceptance criterion, as a number.  Legs priced a hundred thousand
// times above the ball buys a controller that will barely move them: the plate
// stays flat to a hundredth of a degree, and the ball simply stays where it
// was put — and rolls straight off if it is nudged.
void test_poor_tuning_visibly_degrades() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    Eigen::VectorXd r_poor(3);
    r_poor << 1e5, 1e5, 1e5;
    const Eigen::MatrixXd K_poor =
        gainFor(e, defaultLqrStateWeights(7), r_poor);

    const SimResult held = runClosedLoop(e, K_poor, kDisplaced, 0.0, 0.0, 20.0);
    ASSERT_TRUE(held.settle_time < 0.0);          // never settles
    ASSERT_TRUE(held.final_radius > 0.05);        // and never really tried
    ASSERT_TRUE(held.final_tilt_deg < 0.1);

    const SimResult kicked = runClosedLoop(
        e, K_poor, Eigen::Vector4d(0.0, 0.0, 0.35, 0.25), 0.0, 0.0, 20.0);
    ASSERT_TRUE(kicked.left_plate);
}

// Why the defaults above are not unit weights, kept as a fact rather than a
// memory.  The rolling model's Coulomb resistance is a dead band: below
// asin(c_rr) of tilt the ball cannot be started at all, and a state feedback
// with no integral term parks inside it.  A future reader will meet this as
// "the ball stops 37 mm off centre and the plate just sits there tilted", and
// it is a property of the plant, not a defect in the loop.
void test_unit_weights_park_in_the_friction_dead_band() {
    const auto models = getBuiltinModels();
    const auto& e = cascadeModel(models);

    const Eigen::MatrixXd K_unit =
        gainFor(e, Eigen::VectorXd::Ones(7), Eigen::VectorXd::Ones(3));
    const SimResult s = runClosedLoop(e, K_unit, kDisplaced, 0.0, 0.0, 20.0);

    ASSERT_TRUE(!s.left_plate);
    ASSERT_TRUE(s.settle_time < 0.0);        // stalls, well short of centre
    ASSERT_TRUE(s.final_radius > 0.02);
    // Stalled against the tilt that rolling resistance exactly cancels — to
    // within a tenth of it, because sign(v) leaves the ball creeping rather
    // than stopping dead.  An order of magnitude below the 8 deg the working
    // tuning uses, which is the whole difference between the two runs.
    ASSERT_NEAR(s.final_tilt_deg, std::asin(kBallFriction) / kDeg, 0.06);
}

}  // namespace

int main() {
    test_state_vector_layout();
    test_state_deviation_follows_the_design_home();
    test_shape_gate();
    test_mechanism_identity_survives_the_float_round_trip();
    test_zero_error_commands_the_home_pose();
    test_command_is_home_minus_k_error();
    test_setpoint_moves_the_target();
    test_command_clamps_to_servo_travel();
    test_wrong_shape_commands_the_home_pose();
    test_servo_lag_is_first_order();
    test_default_weights_balance_the_ball();
    test_default_weights_reject_a_nudge();
    test_setpoint_is_tracked();
    test_looser_weights_are_visibly_sluggish();
    test_poor_tuning_visibly_degrades();
    test_unit_weights_park_in_the_friction_dead_band();
    std::printf("test_auto_balance: all passed\n");
    return 0;
}
