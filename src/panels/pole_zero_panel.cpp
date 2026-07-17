// src/panels/pole_zero_panel.cpp
#include "pole_zero_panel.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

void drawPoleZeroPanel(AppState& state) {
    if (!state.show_pole_zero) return;
    ImGui::Begin("Pole-Zero / Root Locus", &state.show_pole_zero);

    // Mode selector
    int mode = static_cast<int>(state.pz_mode);
    ImGui::RadioButton("Pole-Zero", &mode, 0); ImGui::SameLine();
    ImGui::RadioButton("Unity FB", &mode, 1); ImGui::SameLine();
    ImGui::RadioButton("State FB", &mode, 2);
    if (static_cast<PZMode>(mode) != state.pz_mode) {
        state.pz_mode = static_cast<PZMode>(mode);
        state.needs_recompute = true;
    }

    // Root locus controls
    if (state.pz_mode == PZMode::UnityFB) {
        ImGui::SliderFloat("K", &state.rl_current_k,
                           state.rl_k_min, state.rl_k_max, "%.2f",
                           ImGuiSliderFlags_Logarithmic);
        if (ImGui::SliderFloat("K min", &state.rl_k_min, 0.0f, 1.0f)) {
            state.needs_recompute = true;
        }
        ImGui::SameLine();
        if (ImGui::SliderFloat("K max", &state.rl_k_max, 1.0f, 1000.0f)) {
            state.needs_recompute = true;
        }
    } else if (state.pz_mode == PZMode::StateFB) {
        ImGui::SliderFloat("\xce\xb1", &state.rl_current_alpha,
                           state.rl_alpha_min, state.rl_alpha_max);
        if (ImGui::SliderFloat("\xce\xb1 max", &state.rl_alpha_max, 0.1f, 10.0f)) {
            state.needs_recompute = true;
        }
    }

    // Plot
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (ImPlot::BeginPlot("##pz_plot", avail, ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("\xcf\x83 (Real)", "j\xcf\x89 (Imag)");

        // Stable region shading (Re < 0)
        {
            auto limits = ImPlot::GetPlotLimits();
            double shade_x[4] = {limits.X.Min, 0, 0, limits.X.Min};
            double shade_y[4] = {limits.Y.Min, limits.Y.Min, limits.Y.Max, limits.Y.Max};
            ImPlot::PushStyleColor(ImPlotCol_Fill, IM_COL32(50, 200, 50, 20));
            ImPlot::PlotShaded("##stable", shade_x, shade_y, 4);
            ImPlot::PopStyleColor();
        }

        // Stability boundary (imaginary axis)
        double zero = 0.0;
        ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(255, 255, 255, 40));
        ImPlot::PlotInfLines("##jw_axis", &zero, 1);
        ImPlot::PopStyleColor();

        // Unit circle
        {
            constexpr int N = 128;
            double cx[N + 1], cy[N + 1];
            for (int i = 0; i <= N; ++i) {
                double th = 2.0 * M_PI * i / N;
                cx[i] = std::cos(th);
                cy[i] = std::sin(th);
            }
            ImPlot::PushStyleColor(ImPlotCol_Line, IM_COL32(255, 255, 255, 30));
            ImPlot::PlotLine("##unit_circle", cx, cy, N + 1);
            ImPlot::PopStyleColor();
        }

        if (state.pz_mode == PZMode::PoleZero) {
            for (int s = 0; s < NUM_SYSTEMS; ++s) {
                if (!state.trace_visible[s] || !state.system_valid[s]) continue;
                const auto& pz = state.pole_zero[s];

                std::vector<double> pre, pim;
                for (const auto& p : pz.poles) {
                    pre.push_back(p.real());
                    pim.push_back(p.imag());
                }
                if (!pre.empty()) {
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Cross, 8,
                                               system_colors[s], 2.0f);
                    char label_p[32];
                    std::snprintf(label_p, sizeof(label_p), "P(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_p, pre.data(), pim.data(),
                                       static_cast<int>(pre.size()));
                }

                std::vector<double> zre, zim;
                for (const auto& z : pz.zeros) {
                    zre.push_back(z.real());
                    zim.push_back(z.imag());
                }
                if (!zre.empty()) {
                    ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 8,
                                               system_colors[s], 2.0f,
                                               ImVec4(0, 0, 0, 0));
                    char label_z[32];
                    std::snprintf(label_z, sizeof(label_z), "Z(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_z, zre.data(), zim.data(),
                                       static_cast<int>(zre.size()));
                }
            }
        } else {
            if (!state.root_locus.empty()) {
                int n_poles = static_cast<int>(state.root_locus[0].poles.size());
                int n_steps = static_cast<int>(state.root_locus.size());

                for (int p = 0; p < n_poles; ++p) {
                    std::vector<double> re(n_steps), im(n_steps);
                    for (int k = 0; k < n_steps; ++k) {
                        re[k] = state.root_locus[k].poles[p].real();
                        im[k] = state.root_locus[k].poles[p].imag();
                    }
                    ImPlot::SetNextLineStyle(system_colors[0], 1.5f);
                    char lbl[16];
                    std::snprintf(lbl, sizeof(lbl), "##rl%d", p);
                    ImPlot::PlotLine(lbl, re.data(), im.data(), n_steps);
                }

                double current_gain = (state.pz_mode == PZMode::UnityFB)
                                          ? state.rl_current_k
                                          : state.rl_current_alpha;
                int best = 0;
                double best_d = 1e30;
                for (int k = 0; k < n_steps; ++k) {
                    double d = std::abs(state.root_locus[k].gain - current_gain);
                    if (d < best_d) { best_d = d; best = k; }
                }
                std::vector<double> cre, cim;
                for (const auto& pole : state.root_locus[best].poles) {
                    cre.push_back(pole.real());
                    cim.push_back(pole.imag());
                }
                ImPlot::SetNextMarkerStyle(ImPlotMarker_Diamond, 10,
                                           system_colors[3], 2.5f);
                ImPlot::PlotScatter("##current", cre.data(), cim.data(),
                                   static_cast<int>(cre.size()));
            }

            const auto& pz = state.pole_zero[0];
            std::vector<double> pre, pim;
            for (const auto& p : pz.poles) {
                pre.push_back(p.real());
                pim.push_back(p.imag());
            }
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Cross, 6,
                                       ImVec4(0.5f, 0.5f, 0.5f, 0.7f), 1.5f);
            ImPlot::PlotScatter("OL Poles", pre.data(), pim.data(),
                               static_cast<int>(pre.size()));
        }

        ImPlot::EndPlot();
    }

    ImGui::End();
}

}  // namespace caliburn
