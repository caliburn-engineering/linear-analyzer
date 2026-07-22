// src/app_state.h
#pragma once

#include "linear_system.h"
#include "analysis/frequency_response.h"
#include "analysis/model_library.h"
#include "analysis/pole_zero.h"
#include "analysis/system_connect.h"
#include "analysis/system_properties.h"
#include "analysis/time_response.h"
#include "imgui.h"

#include <cstring>

namespace caliburn {

constexpr int NUM_SYSTEMS = 4;

inline const char* system_names[NUM_SYSTEMS] = {
    "Plant", "Controller", "Open-Loop", "Closed-Loop"};

inline ImVec4 system_colors[NUM_SYSTEMS] = {
    {0.22f, 0.74f, 0.97f, 1.0f},  // Plant — blue  #38bdf8
    {0.66f, 0.33f, 0.97f, 1.0f},  // Controller — purple  #a855f7
    {0.20f, 0.83f, 0.60f, 1.0f},  // Open-Loop — green  #34d399
    {0.98f, 0.45f, 0.09f, 1.0f},  // Closed-Loop — orange  #f97316
};

inline ImU32 system_colors_u32[NUM_SYSTEMS] = {
    IM_COL32(56, 189, 248, 255),
    IM_COL32(168, 85, 247, 255),
    IM_COL32(52, 211, 153, 255),
    IM_COL32(249, 115, 22, 255),
};

enum class ControllerType { None, StateSpace, GainMatrix };
enum class InputType { Step, Impulse, Ramp };
enum class PZMode { PoleZero, UnityFB, StateFB };

struct AppState {
    // --- Model ---
    int preset_index = 0;
    LinearSystem plant;
    ControllerType ctrl_type = ControllerType::None;
    LinearSystem ctrl_ss;
    Eigen::MatrixXd ctrl_K;

    // --- Derived systems ---
    LinearSystem systems[NUM_SYSTEMS];
    bool system_valid[NUM_SYSTEMS] = {true, false, false, false};

    // --- Channel selection ---
    int output_i = 0;
    int input_j = 0;
    bool show_all_channels = false;

    // --- Analysis results ---
    FrequencyResponse bode[NUM_SYSTEMS];
    PoleZeroResult pole_zero[NUM_SYSTEMS];
    TimeResponse time_resp[NUM_SYSTEMS];
    PropertyResult controllability{};
    PropertyResult observability{};
    PropertyResult cl_controllability{};
    PropertyResult cl_observability{};
    std::vector<RootLocusPoint> root_locus;

    // --- Physical parameters (mutable copy of current preset's params) ---
    std::vector<PhysicalParam> current_params;

    // All-channel results (populated when show_all_channels is true)
    std::vector<FrequencyResponse> bode_grid[NUM_SYSTEMS];

    // --- Recompute flag ---
    bool needs_recompute = true;

    // --- Panel visibility ---
    bool show_pole_zero = true;
    bool show_bode = true;
    bool show_nyquist = true;
    bool show_time_response = true;

    // --- Trace visibility ---
    bool trace_visible[NUM_SYSTEMS] = {true, false, false, false};

    // --- Bode settings ---
    float freq_min_hz = 0.01f;
    float freq_max_hz = 100.0f;
    int num_freq_points = 500;

    // --- Time response settings ---
    InputType input_type = InputType::Step;
    float amplitude = 1.0f;
    float slope = 1.0f;
    float duration = 10.0f;
    float dt_sim = 0.01f;
    int time_input_j = 0;

    // --- Pole-zero / root locus settings ---
    PZMode pz_mode = PZMode::PoleZero;
    float rl_k_min = 0.0f;
    float rl_k_max = 100.0f;
    int rl_num_points = 200;
    float rl_current_k = 1.0f;
    float rl_alpha_min = 0.0f;
    float rl_alpha_max = 2.0f;
    float rl_current_alpha = 1.0f;

    // --- Transfer function parameterization (1st/2nd order SISO) ---
    bool plant_is_1st_order = false;
    bool plant_is_2nd_order = false;
    float tf_K = 1.0f;
    float tf_tau = 1.0f;
    float tf_wn = 1.0f;
    float tf_zeta = 0.7f;

    // --- Text fields for matrix entry ---
    char A_text[1024] = "";
    char B_text[1024] = "";
    char C_text[1024] = "";
    char D_text[1024] = "";
    char Ac_text[1024] = "";
    char Bc_text[1024] = "";
    char Cc_text[1024] = "";
    char Dc_text[1024] = "";
    char K_text[1024] = "";

    // --- Cursors (for synchronized vertical lines + value readout) ---
    double time_cursor_x = -1.0;
    double bode_cursor_x = -1.0;
    double time_x_min = 0.0;
    double time_x_max = 10.0;
};

// Helper: write matrix to text buffer for ImGui display
inline void matrixToTextBuf(const Eigen::MatrixXd& mat, char* buf,
                            size_t buf_size) {
    std::string s = matrixToString(mat);
    std::strncpy(buf, s.c_str(), buf_size - 1);
    buf[buf_size - 1] = '\0';
}

}  // namespace caliburn
