// src/analysis/loop_diagnostics.cpp
#include "loop_diagnostics.h"
#include "frequency_response.h"

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace caliburn {

namespace {
constexpr double kPi = 3.14159265358979323846;

void pairedSets(const std::vector<Loop>& loops,
                std::vector<int>& outs, std::vector<int>& ins) {
    for (const auto& l : loops) {
        if (std::find(outs.begin(), outs.end(), l.out) == outs.end())
            outs.push_back(l.out);
        if (std::find(ins.begin(), ins.end(), l.in) == ins.end())
            ins.push_back(l.in);
    }
    std::sort(outs.begin(), outs.end());
    std::sort(ins.begin(), ins.end());
}

// The pairing permutation Pi, on the paired sub-matrix's own indices.
Eigen::MatrixXd pairingMatrix(const std::vector<Loop>& loops,
                              const std::vector<int>& outs,
                              const std::vector<int>& ins) {
    const int q = static_cast<int>(outs.size());
    Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(q, q);
    for (const auto& l : loops) {
        const int r =
            static_cast<int>(std::find(outs.begin(), outs.end(), l.out) - outs.begin());
        const int c =
            static_cast<int>(std::find(ins.begin(), ins.end(), l.in) - ins.begin());
        Pi(r, c) = 1.0;
    }
    return Pi;
}

Eigen::MatrixXd rgaOf(const Eigen::MatrixXcd& sub) {
    // L = G o (G^-1)^T on the COMPLEX matrix.  .transpose(), NEVER .adjoint()
    // — the two give wildly different answers and no real-matrix test can tell
    // them apart (S&P footnote 5, p. 83).
    return (sub.array() * sub.inverse().transpose().array()).real();
}
}  // anonymous namespace

Eigen::MatrixXcd evalPlantMatrix(const LinearSystem& sys, double freq_hz) {
    const std::complex<double> s(0.0, 2.0 * kPi * freq_hz);
    Eigen::MatrixXcd G(sys.outputs(), sys.inputs());
    for (int i = 0; i < sys.outputs(); ++i)
        for (int j = 0; j < sys.inputs(); ++j)
            G(i, j) = evalTransferFunction(sys, i, j, s);
    return G;
}

bool isStructurallyDeadChannel(const LinearSystem& sys, int output_i, int input_j,
                               double freq_min_hz, double freq_max_hz,
                               int num_points, double tol) {
    const double log_min = std::log10(freq_min_hz);
    const double log_max = std::log10(freq_max_hz);
    const int n = std::max(num_points, 2);
    for (int k = 0; k < n; ++k) {
        const double f =
            std::pow(10.0, log_min + (log_max - log_min) * k / (n - 1));
        const std::complex<double> s(0.0, 2.0 * kPi * f);
        if (std::abs(evalTransferFunction(sys, output_i, input_j, s)) > tol)
            return false;
    }
    return true;
}

LoopDiagnostics computeLoopDiagnostics(const LinearSystem& plant,
                                       const std::vector<Loop>& loops,
                                       const std::vector<float>& output_scales,
                                       double freq_hz,
                                       double freq_min_hz, double freq_max_hz,
                                       int num_points) {
    LoopDiagnostics d;
    const int p = plant.outputs();

    for (float v : output_scales)
        if (std::abs(v - 1.0f) > 1e-6f) { d.scales_at_default = false; break; }

    if (loops.empty()) {
        // Nothing paired: there is no input to take a column from.
        d.headline = "No loops \xe2\x80\x94 controller is None";
        return d;
    }

    const Eigen::MatrixXcd G = evalPlantMatrix(plant, freq_hz);

    d.loop_lambda.assign(loops.size(), std::nan(""));
    d.loop_dead.assign(loops.size(), 0);
    for (std::size_t k = 0; k < loops.size(); ++k) {
        d.loop_dead[k] = isStructurallyDeadChannel(
            plant, loops[k].out, loops[k].in,
            freq_min_hz, freq_max_hz, num_points) ? 1 : 0;
    }

    std::vector<int> outs, ins;
    pairedSets(loops, outs, ins);

    // The RGA exists only on a SQUARE paired sub-matrix of size >= 2.  A 1x1
    // sub-matrix yields L = [1] for every plant at every w, so a green 1.00
    // would be actively misleading.  Note this is a property of the LOOP LIST,
    // not the plant shape: it also excludes a multi-input plant with one loop,
    // and the MISO fan-out whose sub-matrix is non-square — the case that
    // passes a naive `loops.size() >= 2` check.  See issues #4, #8.
    if (outs.size() == ins.size() && outs.size() >= 2) {
        const int q = static_cast<int>(outs.size());
        Eigen::MatrixXcd sub(q, q);
        for (int r = 0; r < q; ++r)
            for (int c = 0; c < q; ++c) sub(r, c) = G(outs[r], ins[c]);

        d.has_rga = true;
        d.lambda = rgaOf(sub);
        d.sub_out = outs;
        d.sub_in = ins;

        // Severity comes from the RGA number, NEVER from |lambda_ii|: for
        // G = [[1,2],[1,1]], |lambda_11| = 1.0000 looks perfect while
        // lambda_11 = -1 and the RGA number is 8.
        d.rga_number = (d.lambda - pairingMatrix(loops, outs, ins))
                           .cwiseAbs().sum();

        for (std::size_t k = 0; k < loops.size(); ++k) {
            const int r = static_cast<int>(
                std::find(outs.begin(), outs.end(), loops[k].out) - outs.begin());
            const int c = static_cast<int>(
                std::find(ins.begin(), ins.end(), loops[k].in) - ins.begin());
            d.loop_lambda[k] = d.lambda(r, c);
        }

        char buf[192];
        std::snprintf(buf, sizeof(buf),
                      "RGA @ %.3g Hz \xe2\x80\x94 RGA number %.2f",
                      freq_hz, d.rga_number);
        d.headline = buf;
        return d;
    }

    // A share over a single output is 1.0 by construction and says nothing, so
    // a 1-output plant gets no diagnostic at all — one cell, the |g| readout
    // only.  This is the glossary's "on a 1x1 plant it is one cell with no
    // diagnostic at all"; the same reasoning that makes a 1x1 RGA vacuous.
    if (p < 2) {
        char buf[192];
        std::snprintf(buf, sizeof(buf),
                      "Single output \xe2\x80\x94 no pairing diagnostic applies");
        d.headline = buf;
        return d;
    }

    // Channel Share.  Not a pairing measure; never coloured or thresholded as
    // one.  lambda_i = |g_i|^2 / sum|g_k|^2 against the user's output scales.
    d.has_share = true;
    d.share_in = loops.front().in;
    d.share.assign(p, 0.0);
    d.mag_db.assign(p, 0.0);
    double total = 0.0;
    for (int i = 0; i < p; ++i) {
        const double scale =
            (i < static_cast<int>(output_scales.size()) && output_scales[i] > 1e-12f)
                ? output_scales[i]
                : 1.0;
        const double gi = std::abs(G(i, d.share_in)) / scale;
        d.share[i] = gi * gi;
        total += d.share[i];
        // The dB column is load-bearing, not decoration: the share is
        // normalised, so it cannot distinguish "drives both outputs strongly"
        // from "drives neither" — both read 50/50.  S&P's C3 has two halves,
        // and for a p x 1 matrix the column sum is exactly 1 by construction,
        // so the input-deletion half is vacuous.  Only the output half
        // survives, and it needs the magnitudes.
        d.mag_db[i] =
            20.0 * std::log10(std::max(std::abs(G(i, d.share_in)), 1e-300));
    }
    for (int i = 0; i < p; ++i)
        d.share[i] = (total > 0.0) ? d.share[i] / total : 0.0;

    // One computation throughout; the MESSAGE adapts to which of the three
    // cases put us here.
    const char* why =
        (plant.inputs() == 1)       ? "single-input plant"
      : (outs.size() != ins.size()) ? "pairing is non-square"
                                    : "one loop of several inputs \xe2\x80\x94 RGA needs two";
    char buf[192];
    std::snprintf(buf, sizeof(buf),
                  "Channel Share of u%d @ %.3g Hz \xe2\x80\x94 %s; scales %s",
                  d.share_in, freq_hz, why,
                  d.scales_at_default ? "at default (all 1)" : "user-set");
    d.headline = buf;
    return d;
}

double rgaNumberAt(const LinearSystem& plant, const std::vector<Loop>& loops,
                   double freq_hz) {
    // Cheap path: the sweep strip only ever wants the number.
    std::vector<int> outs, ins;
    pairedSets(loops, outs, ins);
    if (outs.size() != ins.size() || outs.size() < 2) return std::nan("");

    const Eigen::MatrixXcd G = evalPlantMatrix(plant, freq_hz);
    const int q = static_cast<int>(outs.size());
    Eigen::MatrixXcd sub(q, q);
    for (int r = 0; r < q; ++r)
        for (int c = 0; c < q; ++c) sub(r, c) = G(outs[r], ins[c]);

    return (rgaOf(sub) - pairingMatrix(loops, outs, ins)).cwiseAbs().sum();
}

}  // namespace caliburn
