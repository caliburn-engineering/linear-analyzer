// tests/test_pole_zero.cpp
#include "test_helpers.h"
#include "analysis/pole_zero.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

// Helper: find the pole closest to a given complex value
std::complex<double> find_nearest(
    const std::vector<std::complex<double>>& vec,
    std::complex<double> target) {
    double best = 1e30;
    std::complex<double> result = vec[0];
    for (const auto& v : vec) {
        double d = std::abs(v - target);
        if (d < best) { best = d; result = v; }
    }
    return result;
}

void test_first_order_pole() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.poles.size(), 1);
    ASSERT_NEAR(result.poles[0].real(), -1.0, 1e-10);
    ASSERT_NEAR(result.poles[0].imag(), 0.0, 1e-10);
    ASSERT_TRUE(result.is_stable);
    ASSERT_TRUE(result.zeros.empty());
}

void test_second_order_complex_poles() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1.4).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.poles.size(), 2);

    auto p1 = find_nearest(result.poles, {-0.7, 0.7141});
    ASSERT_NEAR(p1.real(), -0.7, 1e-4);
    ASSERT_NEAR(std::abs(p1.imag()), 0.7141, 1e-3);
    ASSERT_TRUE(result.is_stable);
}

void test_unstable_system() {
    caliburn::LinearSystem sys;
    sys.A = Eigen::MatrixXd::Zero(4, 4);
    sys.A(0, 1) = 1;
    sys.A(1, 2) = -0.981;
    sys.A(2, 3) = 1;
    sys.A(3, 2) = 21.582;
    sys.B = (Eigen::MatrixXd(4, 1) << 0, 1, 0, -2).finished();
    sys.C = Eigen::MatrixXd::Zero(2, 4);
    sys.C(0, 0) = 1; sys.C(1, 2) = 1;
    sys.D = Eigen::MatrixXd::Zero(2, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_TRUE(!result.is_stable);

    auto p = find_nearest(result.poles, {4.645, 0});
    ASSERT_NEAR(p.real(), 4.645, 0.01);
}

void test_transmission_zeros() {
    // G(s) = (s+3)/((s+1)(s+2)), zero at s = -3
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -2, -3).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 3, 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::computePoleZero(sys, 0, 0);

    auto p1 = find_nearest(result.poles, {-1, 0});
    auto p2 = find_nearest(result.poles, {-2, 0});
    ASSERT_NEAR(p1.real(), -1.0, 1e-6);
    ASSERT_NEAR(p2.real(), -2.0, 1e-6);

    ASSERT_EQ((int)result.zeros.size(), 1);
    ASSERT_NEAR(result.zeros[0].real(), -3.0, 1e-6);
    ASSERT_NEAR(result.zeros[0].imag(), 0.0, 1e-6);
}

void test_zeros_with_nonzero_D() {
    // G(s) = (s+1)/(s+2), zero at s = -1
    // State-space: A=[-2], B=[1], C=[-1], D=[1]
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -2).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.D = (Eigen::MatrixXd(1, 1) << 1).finished();

    auto result = caliburn::computePoleZero(sys, 0, 0);
    ASSERT_EQ((int)result.zeros.size(), 1);
    ASSERT_NEAR(result.zeros[0].real(), -1.0, 1e-6);
}

void test_root_locus_endpoints() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(1, 1) << -1).finished();
    sys.B = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.C = (Eigen::MatrixXd(1, 1) << 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto locus = caliburn::computeRootLocus(sys, 0, 0, 0.0, 50.0, 100);
    ASSERT_EQ((int)locus.size(), 100);

    ASSERT_NEAR(locus.front().poles[0].real(), -1.0, 1e-6);
    ASSERT_NEAR(locus.back().poles[0].real(), -51.0, 1e-6);
}

void test_state_feedback_locus() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    Eigen::MatrixXd K(1, 2);
    K << 1, 0;

    auto locus = caliburn::computeStateFeedbackLocus(sys, K, 0.0, 3.0, 50);
    ASSERT_EQ((int)locus.size(), 50);

    auto p0 = find_nearest(locus.front().poles, {-0.5, 0.866});
    ASSERT_NEAR(p0.real(), -0.5, 1e-3);
    ASSERT_NEAR(std::abs(p0.imag()), 0.866, 0.01);

    double expected_imag = std::sqrt(15.0) / 2.0;
    auto p3 = find_nearest(locus.back().poles, {-0.5, expected_imag});
    ASSERT_NEAR(p3.real(), -0.5, 1e-3);
    ASSERT_NEAR(std::abs(p3.imag()), expected_imag, 0.02);
}

int main() {
    test_first_order_pole();
    test_second_order_complex_poles();
    test_unstable_system();
    test_transmission_zeros();
    test_zeros_with_nonzero_D();
    test_root_locus_endpoints();
    test_state_feedback_locus();
    std::printf("All pole_zero tests passed.\n");
    return 0;
}
