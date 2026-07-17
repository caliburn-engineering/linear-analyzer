// src/panels/model_panel.cpp
#include "model_panel.h"
#include <cstdio>

namespace caliburn {

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
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        state.output_i = 0;
        state.input_j = 0;
        state.needs_recompute = true;
    }

    if (!presets.empty()) {
        ImGui::TextDisabled("%s", presets[state.preset_index].description.c_str());
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

    // Plant stability
    if (state.system_valid[0]) {
        showStability(state.pole_zero[0]);
    }

    // Closed-loop stability (when controller is active)
    if (state.system_valid[3]) {
        ImGui::SameLine();
        ImGui::Text("CL:");
        ImGui::SameLine();
        showStability(state.pole_zero[3]);
    }

    // --- Controller section ---
    ImGui::SeparatorText("Controller");
    const char* ctrl_types[] = {"None", "State-Space", "Gain Matrix K"};
    int ctrl_idx = static_cast<int>(state.ctrl_type);
    if (ImGui::Combo("Type", &ctrl_idx, ctrl_types, 3)) {
        state.ctrl_type = static_cast<ControllerType>(ctrl_idx);
        state.needs_recompute = true;
        // Update trace visibility
        if (state.ctrl_type != ControllerType::None) {
            state.trace_visible[3] = true;  // show closed-loop
            if (state.ctrl_type == ControllerType::StateSpace) {
                state.trace_visible[1] = true;
                state.trace_visible[2] = true;
            }
        } else {
            state.trace_visible[1] = false;
            state.trace_visible[2] = false;
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
            if (Ac && Bc && Cc && Dc) {
                state.ctrl_ss = {*Ac, *Bc, *Cc, *Dc};
                state.needs_recompute = true;
            }
        }
    } else if (state.ctrl_type == ControllerType::GainMatrix) {
        ImGui::InputText("K", state.K_text, sizeof(state.K_text));
        if (ImGui::Button("Apply K")) {
            auto K = parseMatrix(state.K_text);
            if (K) {
                state.ctrl_K = *K;
                state.needs_recompute = true;
            }
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

    // --- Trace toggles ---
    ImGui::SeparatorText("Traces");
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (!state.system_valid[s]) continue;
        ImGui::PushStyleColor(ImGuiCol_CheckMark, system_colors[s]);
        ImGui::Checkbox(system_names[s], &state.trace_visible[s]);
        ImGui::PopStyleColor();
        if (s < NUM_SYSTEMS - 1 && state.system_valid[s + 1])
            ImGui::SameLine();
    }

    // --- Reset ---
    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 1.0f));
    if (ImGui::Button("Reset All", ImVec2(-1, 0))) {
        state = AppState{};
        state.plant = presets[0].system;
        matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
        matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
        matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
        matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));
        state.needs_recompute = true;
    }
    ImGui::PopStyleColor();

    ImGui::End();
}

}  // namespace caliburn
