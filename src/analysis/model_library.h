// src/analysis/model_library.h
#pragma once

#include "linear_system.h"
#include <Eigen/Core>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace caliburn {

struct PhysicalParam {
    std::string name;      // Display name ("Mass")
    std::string symbol;    // Short symbol ("m")
    std::string unit;      // Unit string ("kg")
    float value;           // Default/current value
    float min_val;         // Slider minimum
    float max_val;         // Slider maximum
    bool logarithmic = false;
};

struct DerivationStep {
    std::string title;     // "Equations of Motion (ODE)"
    std::string content;   // Multi-line explanation with equations
};

struct ModelEntry {
    std::string name;
    std::string description;
    LinearSystem system;

    // Physical modeling (optional — empty means no physical params)
    std::vector<PhysicalParam> params;
    std::vector<DerivationStep> derivation;
    std::function<LinearSystem(const std::vector<PhysicalParam>&)> builder;
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
