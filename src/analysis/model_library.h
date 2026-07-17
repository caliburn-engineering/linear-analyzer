// src/analysis/model_library.h
#pragma once

#include "../linear_system.h"
#include <Eigen/Core>
#include <optional>
#include <string>
#include <vector>

namespace caliburn {

struct ModelEntry {
    std::string name;
    std::string description;
    LinearSystem system;
};

// Built-in preset models.
std::vector<ModelEntry> getBuiltinModels();

// Parse a matrix from a MATLAB-style string: "0 1; -2 -3"
// Rows separated by ';', columns by whitespace.
// Returns nullopt on parse failure.
std::optional<Eigen::MatrixXd> parseMatrix(const std::string& text);

// Format a matrix as MATLAB-style string: "0 1; -2 -3"
std::string matrixToString(const Eigen::MatrixXd& mat);

}  // namespace caliburn
