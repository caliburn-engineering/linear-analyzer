// tests/test_controller_builders.cpp
#include "analysis/controller_builders.h"
#include "analysis/frequency_response.h"
#include "analysis/system_connect.h"
#include "test_helpers.h"

#include <complex>
#include <cstdio>
#include <vector>

using namespace caliburn;
using cd = std::complex<double>;

namespace {

// Evaluate off the poles: a PI block has a pole at s = 0, so the DC check used
// by test_system_connect.cpp is invalid here.  This set spans below, at and
// above the derivative-filter corner (tau_f = 0.1 -> 10 rad/s).
const cd kS[3] = {cd(0.0, 0.1), cd(0.0, 1.0), cd(0.0, 10.0)};
constexpr double kTol = 1e-10;  // matches the algebraic-identity tests

void assertChannel(const LinearSystem& sys, int i, int j, cd expected, cd s) {
    const cd actual = evalTransferFunction(sys, i, j, s);
    ASSERT_NEAR(actual.real(), expected.real(), kTol);
    ASSERT_NEAR(actual.imag(), expected.imag(), kTol);
}

cd pidTF(const PIDParams& p, cd s) {
    return cd(p.Kp, 0.0) + cd(p.Ki, 0.0) / s
         + cd(p.Kd, 0.0) * s / (cd(p.tau_f, 0.0) * s + 1.0);
}

cd sectionTF(double alpha, double T, cd s) {
    return (s + 1.0 / T) / (s + 1.0 / (alpha * T));
}

cd leadLagTF(const LeadLagParams& p, cd s) {
    cd v(p.Kc, 0.0);
    switch (p.mode) {
        case CompensatorMode::Lead:
            v *= sectionTF(p.alpha_lead, p.T_lead, s); break;
        case CompensatorMode::Lag:
            v *= sectionTF(p.alpha_lag, p.T_lag, s); break;
        case CompensatorMode::LeadLag:
            v *= sectionTF(p.alpha_lead, p.T_lead, s)
               * sectionTF(p.alpha_lag, p.T_lag, s); break;
    }
    return v;
}

// Hand-built shape fixtures, not getBuiltinModels().  Preset numerics belong to
// test_model_library.cpp; coupling these tests to them means a physical-
// parameter tweak breaks an unrelated suite.  Shape is what the assembler
// cares about.
LinearSystem makePlant(int p, int m) {
    LinearSystem g;
    g.A = Eigen::MatrixXd::Constant(1, 1, -1.0);
    g.B = Eigen::MatrixXd::Ones(1, m);
    g.C = Eigen::MatrixXd::Ones(p, 1);
    g.D = Eigen::MatrixXd::Zero(p, m);
    return g;
}

// --------------------------------------------------------------- SISO builders

void test_pid_transfer_functions() {
    PIDParams p_only{2.5f, 0.0f, 0.0f, 0.1f};
    PIDParams pi{1.0f, 3.0f, 0.0f, 0.1f};
    PIDParams pd{1.0f, 0.0f, 0.25f, 0.1f};
    PIDParams pid{1.5f, 3.0f, 0.25f, 0.1f};

    for (const auto& p : {p_only, pi, pd, pid}) {
        const LinearSystem c = buildPID(p);
        ASSERT_EQ(c.outputs(), 1);
        ASSERT_EQ(c.inputs(), 1);
        for (const cd& s : kS) assertChannel(c, 0, 0, pidTF(p, s), s);
    }
}

void test_pid_minimal_realization() {
    // The rule the transfer-function oracle cannot see.
    ASSERT_EQ(buildPID(PIDParams{2.5f, 0.0f, 0.0f, 0.1f}).states(), 0);
    ASSERT_EQ(buildPID(PIDParams{1.0f, 3.0f, 0.0f, 0.1f}).states(), 1);
    ASSERT_EQ(buildPID(PIDParams{1.0f, 0.0f, 0.25f, 0.1f}).states(), 1);
    ASSERT_EQ(buildPID(PIDParams{1.5f, 3.0f, 0.25f, 0.1f}).states(), 2);
}

void test_leadlag_transfer_functions() {
    LeadLagParams lead{CompensatorMode::Lead, 2.0f, 0.1f, 1.0f, 10.0f, 1.0f};
    LeadLagParams lag{CompensatorMode::Lag, 2.0f, 0.1f, 1.0f, 10.0f, 0.5f};
    LeadLagParams both{CompensatorMode::LeadLag, 1.5f, 0.2f, 0.3f, 5.0f, 2.0f};

    for (const auto& p : {lead, lag, both}) {
        const LinearSystem c = buildLeadLag(p);
        for (const cd& s : kS) assertChannel(c, 0, 0, leadLagTF(p, s), s);
    }

    ASSERT_EQ(buildLeadLag(lead).states(), 1);
    ASSERT_EQ(buildLeadLag(lag).states(), 1);
    ASSERT_EQ(buildLeadLag(both).states(), 2);
}

void test_leadlag_degenerate_alpha() {
    // alpha == 1 zeroes the output coefficient -> pure gain, no state.
    LeadLagParams unity{CompensatorMode::Lead, 3.0f, 1.0f, 1.0f, 10.0f, 1.0f};
    ASSERT_EQ(buildLeadLag(unity).states(), 0);
    for (const cd& s : kS) assertChannel(buildLeadLag(unity), 0, 0, cd(3.0, 0.0), s);

    LeadLagParams half{CompensatorMode::LeadLag, 1.0f, 1.0f, 1.0f, 5.0f, 2.0f};
    ASSERT_EQ(buildLeadLag(half).states(), 1);

    LeadLagParams none{CompensatorMode::LeadLag, 1.0f, 1.0f, 1.0f, 1.0f, 2.0f};
    ASSERT_EQ(buildLeadLag(none).states(), 0);
}

// ------------------------------------------------------------------ assembler

void test_assembler_shape_and_placement() {
    const LinearSystem plant = makePlant(2, 2);
    PIDParams a{1.0f, 2.0f, 0.0f, 0.1f};   // 1 state
    PIDParams b{3.0f, 0.0f, 0.5f, 0.1f};   // 1 state
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{1, 1, b, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), plant.inputs());   // m
    ASSERT_EQ(c.inputs(), plant.outputs());   // p
    ASSERT_EQ(c.states(), 2);

    // Loop k occupies channel (in_k, out_k) — the controller is transposed
    // relative to the plant.
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 1, 1, pidTF(b, s), s);
        assertChannel(c, 0, 1, cd(0.0, 0.0), s);
        assertChannel(c, 1, 0, cd(0.0, 0.0), s);
    }
}

void test_assembler_unpaired_channel_is_inert() {
    // SIMO: 1 input, 2 outputs, one loop y0 -> u0.  Channel (0,1) of the 1x2
    // controller must be exactly zero at every frequency — this is the "inert
    // column" the whole loop-pairing model rests on, and nothing else tests it.
    const LinearSystem plant = makePlant(2, 1);
    PIDParams a{2.0f, 1.0f, 0.0f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), 1);
    ASSERT_EQ(c.inputs(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 0, 1, cd(0.0, 0.0), s);
    }
}

void test_assembler_miso_fan_in_sums() {
    // Two loops into one input — the primary SIMO control structure, e.g. the
    // Inverted Pendulum's u = PID_angle(e1) + PID_pos(e0).
    const LinearSystem plant = makePlant(2, 1);
    PIDParams a{1.0f, 2.0f, 0.0f, 0.1f};
    PIDParams b{4.0f, 0.0f, 0.5f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{1, 0, b, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), 1);
    ASSERT_EQ(c.inputs(), 2);
    ASSERT_EQ(c.states(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 0, 1, pidTF(b, s), s);
    }
}

void test_assembler_duplicates_sum() {
    // Unreachable through the pairing grid, but correct for a loop list built
    // programmatically — and it pins the += on D.
    const LinearSystem plant = makePlant(1, 1);
    PIDParams a{1.5f, 2.0f, 0.25f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{0, 0, a, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 4);
    for (const cd& s : kS) assertChannel(c, 0, 0, 2.0 * pidTF(a, s), s);
}

void test_assembler_empty_and_out_of_range() {
    const LinearSystem plant = makePlant(2, 1);

    const LinearSystem empty = buildLoopController({}, LoopKind::PID, plant);
    ASSERT_EQ(empty.outputs(), 1);
    ASSERT_EQ(empty.inputs(), 2);
    ASSERT_EQ(empty.states(), 0);
    ASSERT_NEAR(empty.D.cwiseAbs().sum(), 0.0, 1e-15);

    PIDParams a{2.0f, 1.0f, 0.0f, 0.1f};
    PIDParams bad{9.0f, 9.0f, 0.0f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{5, 0, bad, {}},
                               Loop{0, 3, bad, {}}};
    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 1);           // only the valid loop contributed
    for (const cd& s : kS) assertChannel(c, 0, 0, pidTF(a, s), s);
}

void test_assembler_leadlag_kind() {
    const LinearSystem plant = makePlant(2, 2);
    LeadLagParams l0{CompensatorMode::Lead, 2.0f, 0.1f, 1.0f, 10.0f, 1.0f};
    LeadLagParams l1{CompensatorMode::Lag, 1.0f, 0.1f, 1.0f, 20.0f, 0.5f};
    std::vector<Loop> loops = {Loop{0, 0, {}, l0}, Loop{1, 1, {}, l1}};

    const LinearSystem c = buildLoopController(loops, LoopKind::LeadLag, plant);
    ASSERT_EQ(c.states(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, leadLagTF(l0, s), s);
        assertChannel(c, 1, 1, leadLagTF(l1, s), s);
    }
}

// ------------------------------------------------------------ integration

void test_series_connect_is_square() {
    // The map's "closed loop stays p x p" claim, end to end, on the SIMO shape
    // where it is least obvious.  Also exercises the n = 0 path: the P-only
    // loop contributes no states.
    const LinearSystem plant = makePlant(2, 1);
    std::vector<Loop> loops = {Loop{0, 0, PIDParams{2.0f, 0.0f, 0.0f, 0.1f}, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 0);

    const LinearSystem open = seriesConnect(c, plant);
    ASSERT_EQ(open.outputs(), 2);
    ASSERT_EQ(open.inputs(), 2);
    ASSERT_EQ(open.outputs(), open.inputs());

    const LinearSystem closed = feedbackConnect(open);
    ASSERT_EQ(closed.outputs(), 2);
    ASSERT_EQ(closed.inputs(), 2);
}

}  // anonymous namespace

int main() {
    test_pid_transfer_functions();
    test_pid_minimal_realization();
    test_leadlag_transfer_functions();
    test_leadlag_degenerate_alpha();
    test_assembler_shape_and_placement();
    test_assembler_unpaired_channel_is_inert();
    test_assembler_miso_fan_in_sums();
    test_assembler_duplicates_sum();
    test_assembler_empty_and_out_of_range();
    test_assembler_leadlag_kind();
    test_series_connect_is_square();
    std::printf("All controller_builders tests passed.\n");
    return 0;
}
