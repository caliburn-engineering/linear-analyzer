// src/analysis/system_properties.cpp
#include "system_properties.h"
#include <Eigen/QR>

namespace caliburn {

PropertyResult checkControllability(const LinearSystem& sys) {
    int n = sys.states();
    int m = sys.inputs();

    Eigen::MatrixXd ctrb(n, n * m);
    Eigen::MatrixXd Ak_B = sys.B;
    for (int k = 0; k < n; ++k) {
        ctrb.block(0, k * m, n, m) = Ak_B;
        Ak_B = sys.A * Ak_B;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(ctrb);
    int rank = qr.rank();

    return {ctrb, rank, n, rank == n};
}

PropertyResult checkObservability(const LinearSystem& sys) {
    int n = sys.states();
    int p = sys.outputs();

    Eigen::MatrixXd obsv(n * p, n);
    Eigen::MatrixXd C_Ak = sys.C;
    for (int k = 0; k < n; ++k) {
        obsv.block(k * p, 0, p, n) = C_Ak;
        C_Ak = C_Ak * sys.A;
    }

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(obsv);
    int rank = qr.rank();

    return {obsv, rank, n, rank == n};
}

}  // namespace caliburn
