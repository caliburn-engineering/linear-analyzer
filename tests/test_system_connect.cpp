// tests/test_system_connect.cpp
#include "test_helpers.h"
#include "analysis/system_connect.h"
#include "analysis/frequency_response.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <cstdio>

caliburn::LinearSystem make_g1() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

caliburn::LinearSystem make_g2() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -2).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 2).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);
    return sys;
}

void test_series_dimensions() {
    auto g1 = make_g1();
    auto g2 = make_g2();
    auto L = caliburn::seriesConnect(g1, g2);
    ASSERT_EQ(L.states(), 2);
    ASSERT_EQ(L.inputs(), 1);
    ASSERT_EQ(L.outputs(), 1);
}

void test_series_transfer_function() {
    auto g1 = make_g1();
    auto g2 = make_g2();
    auto L = caliburn::seriesConnect(g1, g2);

    auto G_dc = caliburn::evalTransferFunction(L, 0, 0, {0, 0});
    ASSERT_NEAR(G_dc.real(), 1.0, 1e-10);
    ASSERT_NEAR(G_dc.imag(), 0.0, 1e-10);

    auto G_j = caliburn::evalTransferFunction(L, 0, 0, {0, 1});
    ASSERT_NEAR(G_j.real(), 0.2, 1e-6);
    ASSERT_NEAR(G_j.imag(), -0.6, 1e-6);
}

void test_series_poles() {
    auto g1 = make_g1();
    auto g2 = make_g2();
    auto L = caliburn::seriesConnect(g1, g2);

    Eigen::EigenSolver<Eigen::MatrixXd> es(L.A, false);
    auto eigs = es.eigenvalues();
    std::vector<double> reals = {eigs(0).real(), eigs(1).real()};
    std::sort(reals.begin(), reals.end());
    ASSERT_NEAR(reals[0], -2.0, 1e-10);
    ASSERT_NEAR(reals[1], -1.0, 1e-10);
}

void test_feedback_closed_loop_poles() {
    auto g1 = make_g1();
    auto g2 = make_g2();
    auto L = caliburn::seriesConnect(g1, g2);
    auto T = caliburn::feedbackConnect(L);

    ASSERT_EQ(T.states(), 2);
    ASSERT_EQ(T.inputs(), 1);
    ASSERT_EQ(T.outputs(), 1);

    Eigen::EigenSolver<Eigen::MatrixXd> es(T.A, false);
    auto eigs = es.eigenvalues();

    ASSERT_NEAR(eigs(0).real(), -1.5, 1e-6);
    ASSERT_NEAR(eigs(1).real(), -1.5, 1e-6);
    double expected_imag = std::sqrt(7.0) / 2.0;
    ASSERT_NEAR(std::abs(eigs(0).imag()), expected_imag, 1e-4);
}

void test_feedback_dc_gain() {
    auto g1 = make_g1();
    auto g2 = make_g2();
    auto L = caliburn::seriesConnect(g1, g2);
    auto T = caliburn::feedbackConnect(L);

    auto G_dc = caliburn::evalTransferFunction(T, 0, 0, {0, 0});
    ASSERT_NEAR(G_dc.real(), 0.5, 1e-10);
}

void test_state_feedback() {
    caliburn::LinearSystem plant;
    plant.A = (Eigen::MatrixXd(2, 2) << 0, 1, -2, -3).finished();
    plant.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    plant.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    plant.D = Eigen::MatrixXd::Zero(1, 1);

    Eigen::MatrixXd K(1, 2);
    K << 1, 0;

    auto cl = caliburn::stateFeedbackClose(plant, K);
    ASSERT_EQ(cl.states(), 2);

    Eigen::EigenSolver<Eigen::MatrixXd> es(cl.A, false);
    auto eigs = es.eigenvalues();

    ASSERT_NEAR(eigs(0).real(), -1.5, 1e-6);
    ASSERT_NEAR(eigs(1).real(), -1.5, 1e-6);
    ASSERT_NEAR(std::abs(eigs(0).imag()), std::sqrt(3.0) / 2.0, 1e-4);
}

int main() {
    test_series_dimensions();
    test_series_transfer_function();
    test_series_poles();
    test_feedback_closed_loop_poles();
    test_feedback_dc_gain();
    test_state_feedback();
    std::printf("All system_connect tests passed.\n");
    return 0;
}
