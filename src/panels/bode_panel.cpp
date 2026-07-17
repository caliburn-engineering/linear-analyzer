// src/panels/bode_panel.cpp
#include "bode_panel.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

namespace {
void plotBodeSingle(const char* suffix, AppState& state) {
    float avail_h = ImGui::GetContentRegionAvail().y;
    float plot_h = (avail_h - 40) / 2.0f;

    // --- Magnitude ---
    char mag_id[64];
    std::snprintf(mag_id, sizeof(mag_id), "Magnitude [dB]%s", suffix);
    if (ImPlot::BeginPlot(mag_id, ImVec2(-1, plot_h))) {
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxes("Frequency [Hz]", "|G| [dB]");
        ImPlot::SetupAxisLimits(ImAxis_X1, state.freq_min_hz, state.freq_max_hz,
                                ImPlotCond_Always);

        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& bode = state.bode[s];
            if (bode.points.empty()) continue;
            int n = static_cast<int>(bode.points.size());
            std::vector<double> f(n), m(n);
            for (int k = 0; k < n; ++k) {
                f[k] = bode.points[k].freq_hz;
                m[k] = bode.points[k].magnitude_db;
            }
            ImPlot::PlotLine(system_names[s], f.data(), m.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                        ImPlotProp_LineWeight, 2.0f));
        }

        double zero_db = 0.0;
        ImPlot::PlotInfLines("##0dB", &zero_db, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        if (state.system_valid[2] && state.trace_visible[2] &&
            state.bode[2].phase_crossover_hz > 0) {
            double x = state.bode[2].phase_crossover_hz;
            double gm = state.bode[2].gain_margin_db;
            ImPlot::PlotInfLines("##gm_line", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(1, 0.39f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
            ImPlot::Annotation(x, 0.0, ImVec4(1, 0.4f, 0.4f, 1),
                               ImVec2(5, -10), true,
                               "GM=%.1f dB", gm);
        }

        ImPlot::EndPlot();
    }

    // --- Phase ---
    char phase_id[64];
    std::snprintf(phase_id, sizeof(phase_id), "Phase [deg]%s", suffix);
    if (ImPlot::BeginPlot(phase_id, ImVec2(-1, plot_h))) {
        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
        ImPlot::SetupAxes("Frequency [Hz]", "\xe2\x88\xa0G [\xc2\xb0]");
        ImPlot::SetupAxisLimits(ImAxis_X1, state.freq_min_hz, state.freq_max_hz,
                                ImPlotCond_Always);

        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& bode = state.bode[s];
            if (bode.points.empty()) continue;
            int n = static_cast<int>(bode.points.size());
            std::vector<double> f(n), p(n);
            for (int k = 0; k < n; ++k) {
                f[k] = bode.points[k].freq_hz;
                p[k] = bode.points[k].phase_deg;
            }
            ImPlot::PlotLine(system_names[s], f.data(), p.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                        ImPlotProp_LineWeight, 2.0f));
        }

        double neg180 = -180.0;
        ImPlot::PlotInfLines("##-180", &neg180, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        if (state.system_valid[2] && state.trace_visible[2] &&
            state.bode[2].gain_crossover_hz > 0) {
            double x = state.bode[2].gain_crossover_hz;
            double pm = state.bode[2].phase_margin_deg;
            ImPlot::PlotInfLines("##pm_line", &x, 1,
                                 ImPlotSpec(ImPlotProp_LineColor,
                                            ImVec4(0.39f, 0.78f, 0.39f, 0.47f),
                                            ImPlotProp_LineWeight, 1.0f));
            ImPlot::Annotation(x, -180.0, ImVec4(0.4f, 0.8f, 0.4f, 1),
                               ImVec2(5, -10), true,
                               "PM=%.1f\xc2\xb0", pm);
        }

        ImPlot::EndPlot();
    }
}
}  // anonymous namespace

void drawBodePanel(AppState& state) {
    if (!state.show_bode) return;
    ImGui::Begin("Bode Plot", &state.show_bode);

    bool freq_changed = false;
    freq_changed |= ImGui::SliderFloat(
        "f min [Hz]", &state.freq_min_hz, 0.001f, 1.0f, "%.3f",
        ImGuiSliderFlags_Logarithmic);
    ImGui::SameLine();
    freq_changed |= ImGui::SliderFloat(
        "f max [Hz]", &state.freq_max_hz, 1.0f, 10000.0f, "%.0f",
        ImGuiSliderFlags_Logarithmic);
    if (freq_changed) state.needs_recompute = true;

    if (state.show_all_channels) {
        int p = state.plant.outputs();
        int m = state.plant.inputs();
        if (ImGui::BeginTable("##bode_grid", m, ImGuiTableFlags_Borders)) {
            for (int i = 0; i < p; ++i) {
                ImGui::TableNextRow();
                for (int j = 0; j < m; ++j) {
                    ImGui::TableNextColumn();
                    ImGui::Text("(%d,%d)", i, j);
                    for (int s = 0; s < NUM_SYSTEMS; ++s) {
                        if (!state.trace_visible[s] || !state.system_valid[s]) continue;
                        int idx = i * m + j;
                        if (idx < (int)state.bode_grid[s].size()) {
                            const auto& bode = state.bode_grid[s][idx];
                            if (bode.points.empty()) continue;
                            int nn = static_cast<int>(bode.points.size());
                            std::vector<double> ff(nn), mm(nn);
                            for (int k = 0; k < nn; ++k) {
                                ff[k] = bode.points[k].freq_hz;
                                mm[k] = bode.points[k].magnitude_db;
                            }
                            char pid[64];
                            std::snprintf(pid, sizeof(pid), "##bode_%d_%d_%d", i, j, s);
                            if (ImPlot::BeginPlot(pid, ImVec2(-1, 120))) {
                                ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                                ImPlot::PlotLine(system_names[s], ff.data(), mm.data(), nn,
                                                 ImPlotSpec(ImPlotProp_LineColor,
                                                            system_colors[s]));
                                ImPlot::EndPlot();
                            }
                        }
                    }
                }
            }
            ImGui::EndTable();
        }
    } else {
        plotBodeSingle("##single", state);
    }

    if (state.system_valid[2] && state.trace_visible[2]) {
        const auto& ol = state.bode[2];
        ImGui::Separator();
        ImGui::Text("Open-Loop Margins:");
        ImGui::SameLine();
        if (ol.phase_crossover_hz > 0) {
            ImGui::Text("GM = %.1f dB @ %.3f Hz",
                        ol.gain_margin_db, ol.phase_crossover_hz);
        } else {
            ImGui::Text("GM = \xe2\x88\x9e");
        }
        ImGui::SameLine();
        ImGui::Text("|");
        ImGui::SameLine();
        if (ol.gain_crossover_hz > 0) {
            ImGui::Text("PM = %.1f\xc2\xb0 @ %.3f Hz",
                        ol.phase_margin_deg, ol.gain_crossover_hz);
        } else {
            ImGui::Text("PM = \xe2\x88\x9e");
        }
    }

    ImGui::End();
}

}  // namespace caliburn
