// tests/test_loop_diagnostics.cpp
#include "analysis/loop_diagnostics.h"
#include "test_helpers.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace caliburn;

namespace {

// The Ball-Balancer shape: G = k/s^2 * [[0,1],[1,0]].  B(2,1) and B(3,0) mean
// input 0 drives output 1's chain, so the plant is exactly anti-diagonal.
LinearSystem makeAntiDiagonalIntegrator(double k) {
    LinearSystem g;
    g.A = Eigen::MatrixXd::Zero(4, 4);
    g.A(0, 2) = 1.0;
    g.A(1, 3) = 1.0;
    g.B = Eigen::MatrixXd::Zero(4, 2);
    g.B(2, 1) = k;
    g.B(3, 0) = k;
    g.C = Eigen::MatrixXd::Zero(2, 4);
    g.C(0, 0) = 1.0;
    g.C(1, 1) = 1.0;
    g.D = Eigen::MatrixXd::Zero(2, 2);
    return g;
}

// Coupled, genuinely complex on the jw axis, invertible.
LinearSystem makeCoupled() {
    LinearSystem g;
    g.A = (Eigen::MatrixXd(2, 2) << -1.0, 0.0, 0.0, -2.0).finished();
    g.B = (Eigen::MatrixXd(2, 2) << 1.0, 1.0, 1.0, 3.0).finished();
    g.C = Eigen::MatrixXd::Identity(2, 2);
    g.D = Eigen::MatrixXd::Zero(2, 2);
    return g;
}

LinearSystem makeSimo() {  // 1 input, 2 outputs
    LinearSystem g;
    g.A = (Eigen::MatrixXd(2, 2) << -1.0, 0.0, 0.0, -5.0).finished();
    g.B = (Eigen::MatrixXd(2, 1) << 1.0, 4.0).finished();
    g.C = Eigen::MatrixXd::Identity(2, 2);
    g.D = Eigen::MatrixXd::Zero(2, 1);
    return g;
}

const std::vector<float> kUnitScales2 = {1.0f, 1.0f};

// An integrator plant is fine, provably: by scaling invariance the RGA is
// [[0,1],[1,0]] at EVERY w > 0, with cond(G) = 1 throughout.  The infinite DC
// gain that made the spec's -C*A^-1*B heuristic undefined cancels exactly here.
void test_antidiagonal_rga_is_exact_at_every_omega() {
    const LinearSystem g = makeAntiDiagonalIntegrator(5.0 * 9.81 / 7.0);
    const std::vector<Loop> identity = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};

    for (double f : {1e-4, 1e-2, 1.0, 1e2, 1e4}) {
        const auto d = computeLoopDiagnostics(g, identity, kUnitScales2,
                                              f, 0.01, 100.0, 60);
        ASSERT_TRUE(d.has_rga);
        ASSERT_TRUE(!d.has_share);
        ASSERT_NEAR(d.lambda(0, 0), 0.0, 1e-9);
        ASSERT_NEAR(d.lambda(0, 1), 1.0, 1e-9);
        ASSERT_NEAR(d.lambda(1, 0), 1.0, 1e-9);
        ASSERT_NEAR(d.lambda(1, 1), 0.0, 1e-9);
        // 4.00 for the naive diagonal pairing — the loudest possible signal.
        ASSERT_NEAR(d.rga_number, 4.0, 1e-9);
    }

    const std::vector<Loop> swapped = {Loop{0, 1, {}, {}}, Loop{1, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, swapped, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(d.rga_number, 0.0, 1e-9);
}

// The transpose/adjoint trap.  Row and column sums of the true RGA are exactly
// 1 (S&P A.4.2); the adjoint form does not preserve that for complex G.  This
// is the only test that can fail if .adjoint() creeps in.
void test_rga_row_and_column_sums_are_one() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> loops = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    const auto d = computeLoopDiagnostics(g, loops, kUnitScales2,
                                          0.5, 0.01, 100.0, 60);
    ASSERT_TRUE(d.has_rga);
    for (int r = 0; r < 2; ++r) ASSERT_NEAR(d.lambda.row(r).sum(), 1.0, 1e-9);
    for (int c = 0; c < 2; ++c) ASSERT_NEAR(d.lambda.col(c).sum(), 1.0, 1e-9);
}

// A 1x1 sub-matrix gives L = [1] for every plant at every w.  Trigger on the
// loop list, not on plant.inputs() == 1.
void test_single_loop_on_square_plant_falls_back_to_share() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, one, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(d.has_share);
}

// MISO fan-out: distinct outputs 2, distinct inputs 1 -> non-square.
// Two loops, so it passes a naive loops.size() >= 2 check.
void test_miso_fan_out_falls_back_to_share() {
    const LinearSystem g = makeSimo();
    const std::vector<Loop> two = {Loop{0, 0, {}, {}}, Loop{1, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, two, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(d.has_share);
    ASSERT_EQ(d.share_in, 0);
}

void test_share_sums_to_one_and_respects_scales() {
    const LinearSystem g = makeSimo();
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};

    const auto flat = computeLoopDiagnostics(g, one, kUnitScales2,
                                             1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(flat.share[0] + flat.share[1], 1.0, 1e-12);
    ASSERT_TRUE(flat.scales_at_default);

    const std::vector<float> scaled = {1.0f, 10.0f};
    const auto s = computeLoopDiagnostics(g, one, scaled, 1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(s.share[0] + s.share[1], 1.0, 1e-12);
    ASSERT_TRUE(!s.scales_at_default);
    // Scaling output 1 down by 10 must raise output 0's share.
    ASSERT_TRUE(s.share[0] > flat.share[0]);
    // The dB column is unscaled, so it does not move.
    ASSERT_NEAR(s.mag_db[1], flat.mag_db[1], 1e-12);
}

void test_structurally_dead_channel() {
    const LinearSystem g = makeAntiDiagonalIntegrator(1.0);
    // The paired identity channels are dead at every w; the off-diagonal ones
    // are not.
    ASSERT_TRUE(isStructurallyDeadChannel(g, 0, 0, 0.01, 100.0, 60));
    ASSERT_TRUE(isStructurallyDeadChannel(g, 1, 1, 0.01, 100.0, 60));
    ASSERT_TRUE(!isStructurallyDeadChannel(g, 0, 1, 0.01, 100.0, 60));
    ASSERT_TRUE(!isStructurallyDeadChannel(g, 1, 0, 0.01, 100.0, 60));

    const std::vector<Loop> identity = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    const auto d = computeLoopDiagnostics(g, identity, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_EQ((int)d.loop_dead[0], 1);
    ASSERT_EQ((int)d.loop_dead[1], 1);
}

void test_empty_loop_list_shows_nothing() {
    const LinearSystem g = makeCoupled();
    const auto d = computeLoopDiagnostics(g, {}, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(!d.has_share);
}

void test_rga_number_sweep_matches_full_computation() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> loops = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    for (double f : {0.05, 0.5, 5.0, 50.0}) {
        const auto d = computeLoopDiagnostics(g, loops, kUnitScales2,
                                              f, 0.01, 100.0, 60);
        ASSERT_NEAR(rgaNumberAt(g, loops, f), d.rga_number, 1e-12);
    }
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};
    ASSERT_TRUE(std::isnan(rgaNumberAt(g, one, 1.0)));
}

}  // anonymous namespace

int main() {
    test_antidiagonal_rga_is_exact_at_every_omega();
    test_rga_row_and_column_sums_are_one();
    test_single_loop_on_square_plant_falls_back_to_share();
    test_miso_fan_out_falls_back_to_share();
    test_share_sums_to_one_and_respects_scales();
    test_structurally_dead_channel();
    test_empty_loop_list_shows_nothing();
    test_rga_number_sweep_matches_full_computation();
    std::printf("All loop_diagnostics tests passed.\n");
    return 0;
}
