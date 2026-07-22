// src/panels/model_panel.cpp
#include "model_panel.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

namespace caliburn {

void extractTFFromSS(AppState& state) {
    const auto& plant = state.plant;
    state.plant_is_1st_order = false;
    state.plant_is_2nd_order = false;

    if (plant.states() == 1 && plant.inputs() == 1 && plant.outputs() == 1) {
        double a = -plant.A(0, 0);  // a = 1/tau
        if (a > 1e-10) {
            state.plant_is_1st_order = true;
            state.tf_tau = static_cast<float>(1.0 / a);
            state.tf_K = static_cast<float>(plant.C(0, 0) * plant.B(0, 0) / a);
        }
    } else if (plant.states() == 2 && plant.inputs() == 1 && plant.outputs() == 1) {
        double neg_wn2 = plant.A(1, 0);  // = -wn^2
        if (neg_wn2 < -1e-10) {
            double wn = std::sqrt(-neg_wn2);
            double zeta = -plant.A(1, 1) / (2.0 * wn);
            double K = plant.C(0, 0) / (wn * wn);
            if (wn > 0 && zeta >= 0) {
                state.plant_is_2nd_order = true;
                state.tf_wn = static_cast<float>(wn);
                state.tf_zeta = static_cast<float>(zeta);
                state.tf_K = static_cast<float>(K);
            }
        }
    }
}

void applyTFToSS(AppState& state) {
    if (state.plant_is_1st_order) {
        double tau = state.tf_tau;
        double K = state.tf_K;
        state.plant.A = (Eigen::MatrixXd(1, 1) << -1.0 / tau).finished();
        state.plant.B = (Eigen::MatrixXd(1, 1) << 1.0 / tau).finished();
        state.plant.C = (Eigen::MatrixXd(1, 1) << K).finished();
        state.plant.D = Eigen::MatrixXd::Zero(1, 1);
    } else if (state.plant_is_2nd_order) {
        double wn = state.tf_wn;
        double zeta = state.tf_zeta;
        double K = state.tf_K;
        state.plant.A = (Eigen::MatrixXd(2, 2) << 0, 1, -wn * wn, -2.0 * zeta * wn).finished();
        state.plant.B = (Eigen::MatrixXd(2, 1) << 0, 1).finished();
        state.plant.C = (Eigen::MatrixXd(1, 2) << K * wn * wn, 0).finished();
        state.plant.D = Eigen::MatrixXd::Zero(1, 1);
    }
    matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
    matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
    matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
    matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
}

namespace {

bool drawMatrixSliders(const char* name, Eigen::MatrixXd& mat,
                       char* text_buf, size_t buf_size) {
    bool changed = false;
    int rows = static_cast<int>(mat.rows());
    int cols = static_cast<int>(mat.cols());

    char table_id[64];
    std::snprintf(table_id, sizeof(table_id), "##sliders_%s", name);

    if (ImGui::BeginTable(table_id, cols)) {
        for (int i = 0; i < rows; ++i) {
            ImGui::TableNextRow();
            for (int j = 0; j < cols; ++j) {
                ImGui::TableNextColumn();
                char slider_id[64];
                std::snprintf(slider_id, sizeof(slider_id), "##%s_%d_%d", name, i, j);
                float val = static_cast<float>(mat(i, j));
                float absv = std::abs(val);
                float range = std::max(10.0f, absv * 3.0f);
                ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::SliderFloat(slider_id, &val, -range, range, "%.3f")) {
                    mat(i, j) = static_cast<double>(val);
                    matrixToTextBuf(mat, text_buf, buf_size);
                    changed = true;
                }
            }
        }
        ImGui::EndTable();
    }
    return changed;
}

}  // anonymous namespace

void drawModelPanel(AppState& state, const std::vector<ModelEntry>& presets) {
    ImGui::Begin("Model Configuration");

    // --- Plant section ---
    ImGui::SeparatorText("Plant G(s)");

    // Preset dropdown
    auto getter = [](void* data, int idx) -> const char* {
        auto* v = static_cast<const std::vector<ModelEntry>*>(data);
        return (*v)[idx].name.c_str();
    };
    if (ImGui::Combo("Preset", &state.preset_index, getter,
                     const_cast<void*>(
                         static_cast<const void*>(&presets)),
                     static_cast<int>(presets.size()))) {
        const auto& model = presets[state.preset_index];
        state.plant = model.system;
        state.current_params = model.params;
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        state.output_i = 0;
        state.input_j = 0;
        extractTFFromSS(state);
        state.needs_recompute = true;
    }

    if (!presets.empty()) {
        ImGui::TextDisabled("%s", presets[state.preset_index].description.c_str());
    }

    // --- Physical parameters ---
    if (!state.current_params.empty()) {
        const auto& preset = presets[state.preset_index];
        if (ImGui::CollapsingHeader("Physical Parameters",
                                     ImGuiTreeNodeFlags_DefaultOpen)) {
            bool param_changed = false;
            for (size_t i = 0; i < state.current_params.size(); ++i) {
                auto& pp = state.current_params[i];
                char label[64];
                std::snprintf(label, sizeof(label), "%s [%s]",
                              pp.name.c_str(), pp.unit.c_str());
                ImGui::TextUnformatted(label);
                ImGui::SameLine();
                char slider_id[32];
                std::snprintf(slider_id, sizeof(slider_id), "##pp_%d",
                              static_cast<int>(i));
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGuiSliderFlags flags = pp.logarithmic
                    ? ImGuiSliderFlags_Logarithmic : 0;
                param_changed |= ImGui::SliderFloat(
                    slider_id, &pp.value, pp.min_val, pp.max_val,
                    "%.4g", flags);
            }
            if (param_changed && preset.builder) {
                state.plant = preset.builder(state.current_params);
                matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
                matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
                matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
                matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
                extractTFFromSS(state);
                state.needs_recompute = true;
            }
        }
    }

    // --- Derivation steps ---
    if (!presets.empty() && !presets[state.preset_index].derivation.empty()) {
        if (ImGui::CollapsingHeader("Derivation")) {
            for (const auto& step : presets[state.preset_index].derivation) {
                if (ImGui::TreeNode(step.title.c_str())) {
                    ImGui::TextWrapped("%s", step.content.c_str());
                    ImGui::TreePop();
                }
            }
        }
    }

    // Matrix text fields
    ImGui::InputText("A", state.A_text, sizeof(state.A_text));
    ImGui::InputText("B", state.B_text, sizeof(state.B_text));
    ImGui::InputText("C", state.C_text, sizeof(state.C_text));
    ImGui::InputText("D", state.D_text, sizeof(state.D_text));

    if (ImGui::Button("Apply Plant")) {
        auto A = parseMatrix(state.A_text);
        auto B = parseMatrix(state.B_text);
        auto C = parseMatrix(state.C_text);
        auto D = parseMatrix(state.D_text);
        if (A && B && C && D) {
            state.plant = {*A, *B, *C, *D};
            state.output_i = 0;
            state.input_j = 0;
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    }

    // System info
    ImGui::Text("n=%d  m=%d  p=%d",
                state.plant.states(), state.plant.inputs(),
                state.plant.outputs());
    ImGui::SameLine();

    // Stability indicator helper lambda
    auto showStability = [](const PoleZeroResult& pz) {
        bool has_positive = false;
        bool has_zero = false;
        for (const auto& p : pz.poles) {
            if (p.real() > 1e-10) has_positive = true;
            else if (std::abs(p.real()) <= 1e-10) has_zero = true;
        }
        if (has_positive) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1), "UNSTABLE");
        } else if (has_zero) {
            ImGui::TextColored(ImVec4(1, 0.8f, 0.2f, 1), "MARGINAL");
        } else {
            ImGui::TextColored(ImVec4(0.3f, 1, 0.3f, 1), "STABLE");
        }
    };

    if (state.system_valid[0]) {
        showStability(state.pole_zero[0]);
    }
    if (state.system_valid[3]) {
        ImGui::SameLine();
        ImGui::Text("CL:");
        ImGui::SameLine();
        showStability(state.pole_zero[3]);
    }

    // --- Transfer function display (1st/2nd order SISO only) ---
    if (state.plant_is_1st_order) {
        ImGui::SeparatorText("Transfer Function");
        ImGui::TextDisabled("G(s) = K / (\xcf\x84s + 1)");
        bool tf_changed = false;
        ImGui::TextUnformatted("K"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_K", &state.tf_K, -10.0f, 10.0f);
        ImGui::TextUnformatted("\xcf\x84"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_tau", &state.tf_tau, 0.01f, 100.0f,
                                          "%.3f", ImGuiSliderFlags_Logarithmic);
        if (tf_changed) {
            applyTFToSS(state);
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    } else if (state.plant_is_2nd_order) {
        ImGui::SeparatorText("Transfer Function");
        ImGui::TextDisabled("G(s) = K\xc2\xb7\xcf\x89n\xc2\xb2 / (s\xc2\xb2 + 2\xce\xb6\xcf\x89n\xc2\xb7s + \xcf\x89n\xc2\xb2)");
        bool tf_changed = false;
        ImGui::TextUnformatted("K"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_K", &state.tf_K, -10.0f, 10.0f);
        ImGui::TextUnformatted("\xcf\x89n"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_wn", &state.tf_wn, 0.01f, 100.0f,
                                          "%.3f", ImGuiSliderFlags_Logarithmic);
        ImGui::TextUnformatted("\xce\xb6"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        tf_changed |= ImGui::SliderFloat("##tf_zeta", &state.tf_zeta, 0.0f, 5.0f);
        if (tf_changed) {
            applyTFToSS(state);
            extractTFFromSS(state);
            state.needs_recompute = true;
        }
    }

    // --- State-space matrix sliders ---
    ImGui::SeparatorText("State-Space Matrices");
    bool slider_changed = false;
    ImGui::TextUnformatted("A:");
    slider_changed |= drawMatrixSliders("A", state.plant.A, state.A_text, sizeof(state.A_text));
    ImGui::TextUnformatted("B:");
    slider_changed |= drawMatrixSliders("B", state.plant.B, state.B_text, sizeof(state.B_text));
    ImGui::TextUnformatted("C:");
    slider_changed |= drawMatrixSliders("C", state.plant.C, state.C_text, sizeof(state.C_text));
    ImGui::TextUnformatted("D:");
    slider_changed |= drawMatrixSliders("D", state.plant.D, state.D_text, sizeof(state.D_text));
    if (slider_changed) {
        extractTFFromSS(state);
        state.needs_recompute = true;
    }

    // --- Controller section ---
    ImGui::SeparatorText("Controller");
    const char* ctrl_types[] = {"None", "State-Space", "Gain Matrix K"};
    int ctrl_idx = static_cast<int>(state.ctrl_type);
    if (ImGui::Combo("Type", &ctrl_idx, ctrl_types, 3)) {
        state.ctrl_type = static_cast<ControllerType>(ctrl_idx);
        state.needs_recompute = true;
        if (state.ctrl_type != ControllerType::None) {
            state.trace_visible[3] = true;
            if (state.ctrl_type == ControllerType::StateSpace) {
                state.trace_visible[1] = true;
            }
        } else {
            state.trace_visible[1] = false;
            state.trace_visible[3] = false;
        }
    }

    if (state.ctrl_type == ControllerType::StateSpace) {
        ImGui::InputText("Ac", state.Ac_text, sizeof(state.Ac_text));
        ImGui::InputText("Bc", state.Bc_text, sizeof(state.Bc_text));
        ImGui::InputText("Cc", state.Cc_text, sizeof(state.Cc_text));
        ImGui::InputText("Dc", state.Dc_text, sizeof(state.Dc_text));
        if (ImGui::Button("Apply Controller")) {
            auto Ac = parseMatrix(state.Ac_text);
            auto Bc = parseMatrix(state.Bc_text);
            auto Cc = parseMatrix(state.Cc_text);
            auto Dc = parseMatrix(state.Dc_text);
            if (Ac && Bc && Cc && Dc &&
                Ac->rows() == Ac->cols() &&
                Bc->rows() == Ac->rows() &&
                Cc->cols() == Ac->cols() &&
                Dc->rows() == Cc->rows() &&
                Dc->cols() == Bc->cols()) {
                state.ctrl_ss = {*Ac, *Bc, *Cc, *Dc};
                state.needs_recompute = true;
            }
        }
        if (state.ctrl_ss.A.size() > 0) {
            bool ctrl_slider_changed = false;
            ImGui::TextUnformatted("Ac:");
            ctrl_slider_changed |= drawMatrixSliders("Ac", state.ctrl_ss.A, state.Ac_text, sizeof(state.Ac_text));
            ImGui::TextUnformatted("Bc:");
            ctrl_slider_changed |= drawMatrixSliders("Bc", state.ctrl_ss.B, state.Bc_text, sizeof(state.Bc_text));
            ImGui::TextUnformatted("Cc:");
            ctrl_slider_changed |= drawMatrixSliders("Cc", state.ctrl_ss.C, state.Cc_text, sizeof(state.Cc_text));
            ImGui::TextUnformatted("Dc:");
            ctrl_slider_changed |= drawMatrixSliders("Dc", state.ctrl_ss.D, state.Dc_text, sizeof(state.Dc_text));
            if (ctrl_slider_changed) state.needs_recompute = true;
        }
    } else if (state.ctrl_type == ControllerType::GainMatrix) {
        ImGui::TextDisabled("K must be %d x %d (inputs x states)",
                            state.plant.inputs(), state.plant.states());
        ImGui::InputText("K", state.K_text, sizeof(state.K_text));
        if (ImGui::Button("Apply K")) {
            auto K = parseMatrix(state.K_text);
            if (K && K->rows() == state.plant.inputs() &&
                K->cols() == state.plant.states()) {
                state.ctrl_K = *K;
                state.needs_recompute = true;
            } else if (K) {
                ImGui::OpenPopup("K_dim_err");
            }
        }
        if (ImGui::BeginPopup("K_dim_err")) {
            ImGui::TextColored(ImVec4(1, 0.3f, 0.3f, 1),
                "K must be %d x %d, got %d x %d",
                state.plant.inputs(), state.plant.states(),
                (int)state.ctrl_K.rows(), (int)state.ctrl_K.cols());
            ImGui::EndPopup();
        }
        if (state.ctrl_K.size() > 0) {
            ImGui::TextUnformatted("K:");
            if (drawMatrixSliders("K", state.ctrl_K, state.K_text, sizeof(state.K_text)))
                state.needs_recompute = true;
        }
    }

    // --- Channel selector ---
    ImGui::SeparatorText("Channel");
    int p = state.plant.outputs();
    int m = state.plant.inputs();
    if (ImGui::SliderInt("Output i", &state.output_i, 0, std::max(p - 1, 0))) {
        state.needs_recompute = true;
    }
    if (ImGui::SliderInt("Input j", &state.input_j, 0, std::max(m - 1, 0))) {
        state.needs_recompute = true;
    }
    if (p <= 3 && m <= 3) {
        if (ImGui::Checkbox("Show All Channels", &state.show_all_channels)) {
            state.needs_recompute = true;
        }
    }

    // --- Trace toggles (skip system 2: open-loop combined) ---
    ImGui::SeparatorText("Traces");
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (s == 2) continue;
        if (!state.system_valid[s]) continue;
        ImGui::PushStyleColor(ImGuiCol_CheckMark, system_colors[s]);
        ImGui::Checkbox(system_names[s], &state.trace_visible[s]);
        ImGui::PopStyleColor();
        // Find next visible system for SameLine
        bool has_next = false;
        for (int ns = s + 1; ns < NUM_SYSTEMS; ++ns) {
            if (ns == 2) continue;
            if (state.system_valid[ns]) { has_next = true; break; }
        }
        if (has_next) ImGui::SameLine();
    }

    // --- Reset ---
    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 1.0f));
    if (ImGui::Button("Reset All", ImVec2(-1, 0))) {
        state = AppState{};
        state.plant = presets[0].system;
        state.current_params = presets[0].params;
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        extractTFFromSS(state);
        state.needs_recompute = true;
    }
    ImGui::PopStyleColor();

    ImGui::End();
}

}  // namespace caliburn
