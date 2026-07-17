// src/analysis/pole_zero.h
#pragma once

#include "../linear_system.h"
#include <Eigen/Core>
#include <complex>
#include <vector>

namespace caliburn {

struct PoleZeroResult {
    std::vector<std::complex<double>> poles;
    std::vector<std::complex<double>> zeros;
    bool is_stable;
};

struct RootLocusPoint {
    double gain;
    std::vector<std::complex<double>> poles;
};

PoleZeroResult computePoleZero(
    const LinearSystem& sys, int output_i, int input_j);

std::vector<RootLocusPoint> computeRootLocus(
    const LinearSystem& sys, int output_i, int input_j,
    double k_min, double k_max, int num_points);

std::vector<RootLocusPoint> computeStateFeedbackLocus(
    const LinearSystem& sys,
    const Eigen::MatrixXd& K,
    double alpha_min, double alpha_max, int num_points);

}  // namespace caliburn
