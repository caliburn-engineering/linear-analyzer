// src/analysis/system_connect.cpp
#include "system_connect.h"
#include <Eigen/QR>
#include <cassert>

namespace caliburn {

LinearSystem seriesConnect(const LinearSystem& sys1,
                           const LinearSystem& sys2) {
    assert(sys1.outputs() == sys2.inputs());
    int n1 = sys1.states(), n2 = sys2.states();
    int n = n1 + n2;
    int m = sys1.inputs();
    int p = sys2.outputs();

    LinearSystem result;
    result.A = Eigen::MatrixXd::Zero(n, n);
    result.A.topLeftCorner(n1, n1) = sys1.A;
    result.A.bottomLeftCorner(n2, n1) = sys2.B * sys1.C;
    result.A.bottomRightCorner(n2, n2) = sys2.A;

    result.B = Eigen::MatrixXd::Zero(n, m);
    result.B.topRows(n1) = sys1.B;
    result.B.bottomRows(n2) = sys2.B * sys1.D;

    result.C = Eigen::MatrixXd::Zero(p, n);
    result.C.leftCols(n1) = sys2.D * sys1.C;
    result.C.rightCols(n2) = sys2.C;

    result.D = sys2.D * sys1.D;

    return result;
}

LinearSystem feedbackConnect(const LinearSystem& open_loop) {
    int p = open_loop.outputs();
    int m = open_loop.inputs();
    assert(p == m && "Feedback requires square system");

    const auto& A = open_loop.A;
    const auto& B = open_loop.B;
    const auto& C = open_loop.C;
    const auto& D = open_loop.D;

    Eigen::MatrixXd IpD = Eigen::MatrixXd::Identity(p, p) + D;
    Eigen::MatrixXd IpD_inv =
        IpD.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(p, p));

    LinearSystem result;
    result.A = A - B * IpD_inv * C;
    result.B = B * IpD_inv;
    result.C = IpD_inv * C;
    result.D = IpD_inv * D;

    return result;
}

LinearSystem stateFeedbackClose(const LinearSystem& plant,
                                const Eigen::MatrixXd& K) {
    LinearSystem result;
    result.A = plant.A - plant.B * K;
    result.B = plant.B;
    result.C = plant.C - plant.D * K;
    result.D = plant.D;
    return result;
}

}  // namespace caliburn
