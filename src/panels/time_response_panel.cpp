// src/panels/time_response_panel.cpp
#include "time_response_panel.h"
#include "panel_utils.h"
#include "implot.h"
#include <cmath>
#include <algorithm>

namespace caliburn {

void drawTimeResponsePanel(AppState& state) {
    if (!state.show_time_response) return;
    ImGui::Begin("Time Response", &state.show_time_response);
    drawHelpMarker(
        "Time Response\n\n"
        "Output y(t) and input u(t) for step, impulse, or ramp.\n\n"
        "Step response characteristics:\n"
        "  Rise time: 10% to 90% of final value.\n"
        "    Faster with higher bandwidth.\n"
        "  Settling time: Time to stay within 2%% of final value.\n"
        "    Related to dominant pole real part.\n"
        "  Overshoot: Peak above steady-state (%%).\n"
        "    Related to damping ratio (zeta).\n"
        "    Approx exp(-pi*zeta/sqrt(1-zeta^2))*100%%.\n"
        "  Steady-state error: Final value offset.\n"
        "    Depends on system type (# of integrators).\n\n"
        "Compare Plant vs Closed-Loop to see controller effect.");
    ImGui::SameLine();

    int itype = static_cast<int>(state.input_type);
    bool changed = false;
    changed |= ImGui::RadioButton("Impulse", &itype, 1); ImGui::SameLine();
    changed |= ImGui::RadioButton("Step", &itype, 0); ImGui::SameLine();
    changed |= ImGui::RadioButton("Ramp", &itype, 2);
    if (changed) {
        state.input_type = static_cast<InputType>(itype);
        state.needs_recompute = true;
    }

    if (ImGui::BeginTable("##time_ctrls", 3)) {
        ImGui::TableNextColumn();
        if (state.input_type == InputType::Ramp) {
            ImGui::TextUnformatted("Slope"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##slope", &state.slope, 0.1f, 10.0f))
                state.needs_recompute = true;
        } else {
            ImGui::TextUnformatted("Amplitude"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##amp", &state.amplitude, 0.1f, 10.0f))
                state.needs_recompute = true;
        }
        ImGui::TableNextColumn();
        ImGui::TextUnformatted("Duration [s]"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::SliderFloat("##dur", &state.duration, 1.0f, 30.0f))
            state.needs_recompute = true;
        ImGui::TableNextColumn();
        // Rule C: the plant is excited on input_j, the closed loop on ref_r.
        // The channel selector lives in the model panel; a parallel
        // time-domain index pair was never deliberate (issue #9).
        ImGui::TextDisabled("plant u%d  /  closed-loop r%d",
                            state.input_j, state.ref_r);
        ImGui::EndTable();
    }

    float avail_h = ImGui::GetContentRegionAvail().y;
    float plot_h = avail_h / 2.0f;

    // Precompute output y-limits
    double out_lo = 1e30, out_hi = -1e30;
    bool has_out = false;
    for (int s : {0, 3}) {
        if (!state.trace_visible[s] || !state.system_valid[s]) continue;
        const auto& resp = state.time_resp[s];
        for (const auto& pt : resp.points) {
            for (int i = 0; i < pt.output.size(); ++i) {
                out_lo = std::min(out_lo, pt.output(i));
                out_hi = std::max(out_hi, pt.output(i));
                has_out = true;
            }
        }
    }
    if (!has_out) { out_lo = -1; out_hi = 1; }
    { double m = std::max((out_hi - out_lo) * 0.1, 0.05); out_lo -= m; out_hi += m; }

    // Precompute input y-limits
    double in_lo = 1e30, in_hi = -1e30;
    bool has_in = false;
    if (state.system_valid[0] && !state.time_resp[0].points.empty()) {
        for (const auto& pt : state.time_resp[0].points) {
            for (int i = 0; i < pt.input.size(); ++i) {
                in_lo = std::min(in_lo, pt.input(i));
                in_hi = std::max(in_hi, pt.input(i));
                has_in = true;
            }
        }
    }
    if (!has_in) { in_lo = -1; in_hi = 1; }
    { double m = std::max((in_hi - in_lo) * 0.1, 0.05); in_lo -= m; in_hi += m; }

    // --- Output plot ---
    if (ImPlot::BeginPlot("Output", ImVec2(-1, plot_h))) {
        ImPlot::SetupAxes("Time [s]", "y(t)");
        ImPlot::SetupAxisLimits(ImAxis_X1, 0, state.duration,
                                ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, out_lo, out_hi, ImPlotCond_Always);
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

                const char* tag = (s == 0) ? "P" : "T";
                char lbl[8];
                std::snprintf(lbl, sizeof(lbl), "y%d%s", oi, tag);
                ImPlot::PlotLine(lbl, t.data(), y.data(), n_pts,
                                 ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                            ImPlotProp_LineWeight, 2.0f));
            }
        }

        if (ImPlot::IsPlotHovered()) {
            double mx = ImPlot::GetPlotMousePos().x;
            state.time_cursor_x = mx;
        }
        if (state.time_cursor_x >= 0) {
            ImPlot::PlotInfLines("##cursor", &state.time_cursor_x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 1, 1, 0.31f)));
        }

        ImPlot::EndPlot();
    }

    // --- Input plot ---
    if (ImPlot::BeginPlot("Input", ImVec2(-1, plot_h))) {
        ImPlot::SetupAxes("Time [s]", "u(t)");
        ImPlot::SetupAxisLimits(ImAxis_X1, 0, state.duration,
                                ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, in_lo, in_hi, ImPlotCond_Always);
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

                char lbl[8];
                std::snprintf(lbl, sizeof(lbl), "u%d", ui);
                ImPlot::PlotLine(lbl, t.data(), u.data(), n_pts,
                                 ImPlotSpec(ImPlotProp_LineColor, system_colors[0],
                                            ImPlotProp_LineWeight, 2.0f));
            }
        }

        if (ImPlot::IsPlotHovered()) {
            double mx = ImPlot::GetPlotMousePos().x;
            state.time_cursor_x = mx;
        }
        if (state.time_cursor_x >= 0) {
            ImPlot::PlotInfLines("##cursor", &state.time_cursor_x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 1, 1, 0.31f)));
        }

        ImPlot::EndPlot();
    }

    // Cursor readout tooltip
    if (state.time_cursor_x >= 0) {
        double tx = state.time_cursor_x;
        ImGui::BeginTooltip();
        ImGui::Text("t = %.3f s", tx);
        for (int s : {0, 3}) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& resp = state.time_resp[s];
            if (resp.points.empty()) continue;
            int n = static_cast<int>(resp.points.size());
            // Find nearest sample
            int best = 0;
            double best_d = 1e30;
            for (int k = 0; k < n; ++k) {
                double d = std::abs(resp.points[k].time - tx);
                if (d < best_d) { best_d = d; best = k; }
            }
            const auto& pt = resp.points[best];
            for (int oi = 0; oi < pt.output.size(); ++oi) {
                const char* tag = (s == 0) ? "P" : "T";
                ImGui::TextColored(system_colors[s], "y%d%s = %.4f", oi, tag, pt.output(oi));
            }
        }
        // Input values
        if (state.system_valid[0] && !state.time_resp[0].points.empty()) {
            const auto& resp = state.time_resp[0];
            int n = static_cast<int>(resp.points.size());
            int best = 0;
            double best_d = 1e30;
            for (int k = 0; k < n; ++k) {
                double d = std::abs(resp.points[k].time - tx);
                if (d < best_d) { best_d = d; best = k; }
            }
            const auto& pt = resp.points[best];
            for (int ui = 0; ui < pt.input.size(); ++ui) {
                ImGui::Text("u%d = %.4f", ui, pt.input(ui));
            }
        }
        ImGui::EndTooltip();
        state.time_cursor_x = -1.0;
    }

    ImGui::End();
}

}  // namespace caliburn
