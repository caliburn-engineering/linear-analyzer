// src/analysis/lqr.cpp
#include "lqr.h"

namespace caliburn {

LqrResult computeLQR(const LinearSystem&, const Eigen::MatrixXd&,
                     const Eigen::MatrixXd&) {
    LqrResult r;
    r.error = "not implemented";
    return r;
}

}  // namespace caliburn
