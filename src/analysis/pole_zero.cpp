// src/analysis/pole_zero.cpp
#include "pole_zero.h"
#include <Eigen/Eigenvalues>
#include <cmath>
#include <limits>

namespace caliburn {

namespace {
void matchPoles(std::vector<std::complex<double>>& poles,
                const std::vector<std::complex<double>>& prev) {
    int n = static_cast<int>(poles.size());
    std::vector<bool> used(n, false);
    std::vector<std::complex<double>> matched(n);
    for (int i = 0; i < n; ++i) {
        double best_dist = std::numeric_limits<double>::max();
        int best_j = 0;
        for (int j = 0; j < n; ++j) {
            if (used[j]) continue;
            double dist = std::abs(poles[j] - prev[i]);
            if (dist < best_dist) {
                best_dist = dist;
                best_j = j;
            }
        }
        matched[i] = poles[best_j];
        used[best_j] = true;
    }
    poles = matched;
}
}  // anonymous namespace

PoleZeroResult computePoleZero(
    const LinearSystem& sys, int output_i, int input_j) {
    PoleZeroResult result;
    int n = sys.states();

    // Poles = eigenvalues of A
    Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A, false);
    result.poles.resize(n);
    result.is_stable = true;
    for (int i = 0; i < n; ++i) {
        result.poles[i] = es.eigenvalues()(i);
        if (result.poles[i].real() >= -1e-10) {
            result.is_stable = false;
        }
    }

    // Transmission zeros
    double D_ij = sys.D(output_i, input_j);

    if (std::abs(D_ij) > 1e-14) {
        Eigen::MatrixXd A_z = sys.A -
            sys.B.col(input_j) * sys.C.row(output_i) / D_ij;
        Eigen::EigenSolver<Eigen::MatrixXd> zes(A_z, false);
        result.zeros.resize(n);
        for (int i = 0; i < n; ++i) {
            result.zeros[i] = zes.eigenvalues()(i);
        }
    } else {
        Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n + 1, n + 1);
        P.topLeftCorner(n, n) = sys.A;
        P.topRightCorner(n, 1) = sys.B.col(input_j);
        P.bottomLeftCorner(1, n) = sys.C.row(output_i);
        P(n, n) = D_ij;

        Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n + 1, n + 1);
        Q.topLeftCorner(n, n) = Eigen::MatrixXd::Identity(n, n);

        Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges(P, Q);
        auto alphas = ges.alphas();
        auto betas = ges.betas();

        for (int k = 0; k < n + 1; ++k) {
            if (std::abs(betas(k)) > 1e-10) {
                result.zeros.push_back(alphas(k) / betas(k));
            }
        }
    }

    return result;
}

std::vector<RootLocusPoint> computeRootLocus(
    const LinearSystem& sys, int output_i, int input_j,
    double k_min, double k_max, int num_points) {
    int n = sys.states();
    double D_ij = sys.D(output_i, input_j);
    std::vector<RootLocusPoint> result;
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        double K = k_min + (k_max - k_min) * step /
                   std::max(num_points - 1, 1);

        double denom = 1.0 + K * D_ij;
        double K_eff = (std::abs(denom) > 1e-12) ? K / denom : K;

        Eigen::MatrixXd A_cl =
            sys.A - K_eff * sys.B.col(input_j) * sys.C.row(output_i);
        Eigen::EigenSolver<Eigen::MatrixXd> es(A_cl, false);

        std::vector<std::complex<double>> poles(n);
        for (int i = 0; i < n; ++i) {
            poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) {
            matchPoles(poles, result.back().poles);
        }

        result.push_back({K, poles});
    }

    return result;
}

std::vector<RootLocusPoint> computeStateFeedbackLocus(
    const LinearSystem& sys,
    const Eigen::MatrixXd& K,
    double alpha_min, double alpha_max, int num_points) {
    int n = sys.states();
    std::vector<RootLocusPoint> result;
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        double alpha = alpha_min + (alpha_max - alpha_min) * step /
                       std::max(num_points - 1, 1);

        Eigen::MatrixXd A_cl = sys.A - alpha * sys.B * K;
        Eigen::EigenSolver<Eigen::MatrixXd> es(A_cl, false);

        std::vector<std::complex<double>> poles(n);
        for (int i = 0; i < n; ++i) {
            poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) {
            matchPoles(poles, result.back().poles);
        }

        result.push_back({alpha, poles});
    }

    return result;
}

}  // namespace caliburn
