// tests/test_system_properties.cpp
#include "test_helpers.h"
#include "analysis/system_properties.h"
#include <cstdio>

void test_controllable_system() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1.4).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::checkControllability(sys);
    ASSERT_TRUE(result.pass);
    ASSERT_EQ(result.rank, 2);
    ASSERT_EQ(result.required_rank, 2);
    ASSERT_EQ(result.matrix.rows(), 2);
    ASSERT_EQ(result.matrix.cols(), 2);
}

void test_uncontrollable_system() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 1, 0, 0, 2).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 1, 0).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 1).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::checkControllability(sys);
    ASSERT_TRUE(!result.pass);
    ASSERT_EQ(result.rank, 1);
    ASSERT_EQ(result.required_rank, 2);
}

void test_observable_system() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1.4).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::checkObservability(sys);
    ASSERT_TRUE(result.pass);
    ASSERT_EQ(result.rank, 2);
}

void test_unobservable_system() {
    caliburn::LinearSystem sys;
    sys.A = (Eigen::MatrixXd(2, 2) << 1, 0, 0, 2).finished();
    sys.B = (Eigen::MatrixXd(2, 1) << 1, 1).finished();
    sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
    sys.D = Eigen::MatrixXd::Zero(1, 1);

    auto result = caliburn::checkObservability(sys);
    ASSERT_TRUE(!result.pass);
    ASSERT_EQ(result.rank, 1);
}

void test_mimo_controllability() {
    caliburn::LinearSystem sys;
    constexpr double k = 5.0 / 7.0, g = 9.81;
    sys.A = Eigen::MatrixXd::Zero(4, 4);
    sys.A(0, 2) = 1; sys.A(1, 3) = 1;
    sys.B = Eigen::MatrixXd::Zero(4, 2);
    sys.B(2, 1) = k * g; sys.B(3, 0) = k * g;
    sys.C = Eigen::MatrixXd::Zero(2, 4);
    sys.C(0, 0) = 1; sys.C(1, 1) = 1;
    sys.D = Eigen::MatrixXd::Zero(2, 2);

    auto ctrb = caliburn::checkControllability(sys);
    ASSERT_TRUE(ctrb.pass);
    ASSERT_EQ(ctrb.rank, 4);
    ASSERT_EQ(ctrb.matrix.rows(), 4);
    ASSERT_EQ(ctrb.matrix.cols(), 8);
}

int main() {
    test_controllable_system();
    test_uncontrollable_system();
    test_observable_system();
    test_unobservable_system();
    test_mimo_controllability();
    std::printf("All system_properties tests passed.\n");
    return 0;
}
