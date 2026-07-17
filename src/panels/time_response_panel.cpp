// src/panels/time_response_panel.cpp
#include "time_response_panel.h"
#include "implot.h"
#include <cmath>
#include <algorithm>

namespace caliburn {

void drawTimeResponsePanel(AppState& state) {
    if (!state.show_time_response) return;
    ImGui::Begin("Time Response", &state.show_time_response);

    int itype = static_cast<int>(state.input_type);
    bool changed = false;
    changed |= ImGui::RadioButton("Step", &itype, 0); ImGui::SameLine();
    changed |= ImGui::RadioButton("Impulse", &itype, 1); ImGui::SameLine();
    changed |= ImGui::RadioButton("Ramp", &itype, 2);
    if (changed) {
        state.input_type = static_cast<InputType>(itype);
        state.needs_recompute = true;
    }

    if (state.input_type == InputType::Ramp) {
        if (ImGui::SliderFloat("Slope", &state.slope, 0.1f, 10.0f))
            state.needs_recompute = true;
    } else {
        if (ImGui::SliderFloat("Amplitude", &state.amplitude, 0.1f, 10.0f))
            state.needs_recompute = true;
    }
    if (ImGui::SliderFloat("Duration [s]", &state.duration, 1.0f, 30.0f))
        state.needs_recompute = true;
    if (ImGui::SliderInt("Input ch", &state.time_input_j, 0,
                         std::max(state.plant.inputs() - 1, 0)))
        state.needs_recompute = true;

    float avail_h = ImGui::GetContentRegionAvail().y;
    float plot_h = avail_h / 2.0f;

    // --- Output plot ---
    if (ImPlot::BeginPlot("Output", ImVec2(-1, plot_h))) {
        ImPlot::SetupAxes("Time [s]", "y(t)");
        ImPlot::SetupAxisLinks(ImAxis_X1, &state.time_x_min, &state.time_x_max);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 0, ImPlotCond_Once);
        ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

        for (int s : {0, 3}) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& resp = state.time_resp[s];
            if (resp.points.empty()) continue;

            int p = static_cast<int>(resp.points[0].output.size());
            int n_pts = static_cast<int>(resp.points.size());

            for (int oi = 0; oi < p; ++oi) {
                std::vector<double> t(n_pts), y(n_pts);
                for (int k = 0; k < n_pts; ++k) {
                    t[k] = resp.points[k].time;
                    y[k] = resp.points[k].output(oi);
                }

                ImPlot::SetNextLineStyle(system_colors[s], 2.0f);
                const char* tag = (s == 0) ? "P" : "T";
                char lbl[8];
                std::snprintf(lbl, sizeof(lbl), "y%d%s", oi, tag);
                ImPlot::PlotLine(lbl, t.data(), y.data(), n_pts);
            }
        }

        if (ImPlot::IsPlotHovered()) {
            double mx = ImPlot::GetPlotMousePos().x;
            ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(255, 255, 255, 80));
            ImPlot::PlotInfLines("##cursor", &mx, 1);
            ImPlot::PopStyleColor();
            state.time_cursor_x = mx;
        }

        ImPlot::EndPlot();
    }

    // --- Input plot ---
    if (ImPlot::BeginPlot("Input", ImVec2(-1, plot_h))) {
        ImPlot::SetupAxes("Time [s]", "u(t)");
        ImPlot::SetupAxisLinks(ImAxis_X1, &state.time_x_min, &state.time_x_max);
        ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 0, ImPlotCond_Once);
        ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

        if (state.system_valid[0] && !state.time_resp[0].points.empty()) {
            const auto& resp = state.time_resp[0];
            int n_pts = static_cast<int>(resp.points.size());
            int m = static_cast<int>(resp.points[0].input.size());

            for (int ui = 0; ui < m; ++ui) {
                std::vector<double> t(n_pts), u(n_pts);
                for (int k = 0; k < n_pts; ++k) {
                    t[k] = resp.points[k].time;
                    u[k] = resp.points[k].input(ui);
                }

                ImPlot::SetNextLineStyle(system_colors[0], 2.0f);
                char lbl[8];
                std::snprintf(lbl, sizeof(lbl), "u%d", ui);
                ImPlot::PlotLine(lbl, t.data(), u.data(), n_pts);
            }
        }

        if (state.time_cursor_x >= 0) {
            ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(255, 255, 255, 80));
            ImPlot::PlotInfLines("##cursor", &state.time_cursor_x, 1);
            ImPlot::PopStyleColor();
            state.time_cursor_x = -1.0;
        }

        ImPlot::EndPlot();
    }

    ImGui::End();
}

}  // namespace caliburn
