// src/analysis/model_library.cpp
#include "model_library.h"
#include <cstdio>
#include <sstream>

namespace caliburn {

std::optional<Eigen::MatrixXd> parseMatrix(const std::string& text) {
    std::vector<std::vector<double>> rows;
    std::istringstream row_stream(text);
    std::string row_str;

    while (std::getline(row_stream, row_str, ';')) {
        std::istringstream col_stream(row_str);
        std::vector<double> row;
        double val;
        while (col_stream >> val) {
            row.push_back(val);
        }
        // Check for non-numeric content remaining
        if (col_stream.fail() && !col_stream.eof()) {
            return std::nullopt;
        }
        if (!row.empty()) {
            rows.push_back(row);
        }
    }

    if (rows.empty()) return std::nullopt;

    size_t cols = rows[0].size();
    for (const auto& r : rows) {
        if (r.size() != cols) return std::nullopt;
    }

    Eigen::MatrixXd mat(rows.size(), cols);
    for (size_t i = 0; i < rows.size(); ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mat(static_cast<int>(i), static_cast<int>(j)) = rows[i][j];
        }
    }
    return mat;
}

std::string matrixToString(const Eigen::MatrixXd& mat) {
    std::string text;
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) text += "; ";
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) text += " ";
            char num[32];
            std::snprintf(num, sizeof(num), "%g", mat(i, j));
            text += num;
        }
    }
    return text;
}

std::vector<ModelEntry> getBuiltinModels() {
    std::vector<ModelEntry> models;

    // 1. Simple Second-Order (ωn=1, ζ=0.7)
    {
        LinearSystem sys;
        sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -1, -1.4).finished();
        sys.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
        sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
        sys.D = Eigen::MatrixXd::Zero(1, 1);
        models.push_back({"Second-Order",
                           "Mass-spring-damper (\xcf\x89n=1, \xce\xb6=0.7)", sys});
    }

    // 2. Ball-Balancer (4 states, 2 inputs, 2 outputs)
    {
        constexpr double k = 5.0 / 7.0;
        constexpr double g = 9.81;
        LinearSystem sys;
        sys.A = Eigen::MatrixXd::Zero(4, 4);
        sys.A(0, 2) = 1.0;
        sys.A(1, 3) = 1.0;
        sys.B = Eigen::MatrixXd::Zero(4, 2);
        sys.B(2, 1) = k * g;
        sys.B(3, 0) = k * g;
        sys.C = Eigen::MatrixXd::Zero(2, 4);
        sys.C(0, 0) = 1.0;
        sys.C(1, 1) = 1.0;
        sys.D = Eigen::MatrixXd::Zero(2, 2);
        models.push_back({"Ball-Balancer",
                           "4 states, 2 inputs, 2 outputs", sys});
    }

    // 3. Inverted Pendulum on Cart (M=1, m=0.1, l=0.5)
    {
        constexpr double M = 1.0, m = 0.1, l = 0.5, g = 9.81;
        LinearSystem sys;
        sys.A = Eigen::MatrixXd::Zero(4, 4);
        sys.A(0, 1) = 1.0;
        sys.A(1, 2) = -m * g / M;
        sys.A(2, 3) = 1.0;
        sys.A(3, 2) = (M + m) * g / (M * l);
        sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / M, 0, -1.0 / (M * l))
                    .finished();
        sys.C = Eigen::MatrixXd::Zero(2, 4);
        sys.C(0, 0) = 1.0;
        sys.C(1, 2) = 1.0;
        sys.D = Eigen::MatrixXd::Zero(2, 1);
        models.push_back({"Inverted Pendulum",
                           "Cart-pendulum (unstable)", sys});
    }

    // 4. Quarter-Car Suspension
    {
        constexpr double ms = 300, mu = 50;
        constexpr double ks = 20000, kt = 200000, bs = 1000;
        LinearSystem sys;
        sys.A = Eigen::MatrixXd::Zero(4, 4);
        sys.A(0, 1) = 1.0;
        sys.A(1, 0) = -ks / ms;
        sys.A(1, 1) = -bs / ms;
        sys.A(1, 2) = ks / ms;
        sys.A(1, 3) = bs / ms;
        sys.A(2, 3) = 1.0;
        sys.A(3, 0) = ks / mu;
        sys.A(3, 1) = bs / mu;
        sys.A(3, 2) = -(ks + kt) / mu;
        sys.A(3, 3) = -bs / mu;
        sys.B = (Eigen::MatrixXd(4, 1) << 0, 0, 0, kt / mu).finished();
        sys.C = Eigen::MatrixXd::Zero(2, 4);
        sys.C(0, 0) = 1.0;
        sys.C(1, 2) = 1.0;
        sys.D = Eigen::MatrixXd::Zero(2, 1);
        models.push_back({"Quarter-Car",
                           "Suspension (sprung + unsprung)", sys});
    }

    // 5. Double Mass-Spring-Damper
    {
        constexpr double m1 = 1, m2 = 1;
        constexpr double k1 = 1, k2 = 1, b1 = 0.1, b2 = 0.1;
        LinearSystem sys;
        sys.A = Eigen::MatrixXd::Zero(4, 4);
        sys.A(0, 1) = 1.0;
        sys.A(1, 0) = -(k1 + k2) / m1;
        sys.A(1, 1) = -(b1 + b2) / m1;
        sys.A(1, 2) = k2 / m1;
        sys.A(1, 3) = b2 / m1;
        sys.A(2, 3) = 1.0;
        sys.A(3, 0) = k2 / m2;
        sys.A(3, 1) = b2 / m2;
        sys.A(3, 2) = -k2 / m2;
        sys.A(3, 3) = -b2 / m2;
        sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / m1, 0, 0).finished();
        sys.C = Eigen::MatrixXd::Zero(2, 4);
        sys.C(0, 0) = 1.0;
        sys.C(1, 2) = 1.0;
        sys.D = Eigen::MatrixXd::Zero(2, 1);
        models.push_back({"Double Mass-Spring",
                           "Two coupled oscillators", sys});
    }

    return models;
}

}  // namespace caliburn
