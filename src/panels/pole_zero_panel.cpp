// src/panels/pole_zero_panel.cpp
#include "pole_zero_panel.h"
#include "panel_utils.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

void drawPoleZeroPanel(AppState& state) {
    if (!state.show_pole_zero) return;
    ImGui::Begin("Pole-Zero / Root Locus", &state.show_pole_zero);
    drawHelpMarker(
        "Pole-Zero / Root Locus Plot\n\n"
        "Poles (x) are roots of det(sI - A) = 0.\n"
        "Zeros (o) are transmission zeros.\n\n"
        "Stability: All poles must be in the left half-plane (Re < 0).\n"
        "Green shading = stable region.\n\n"
        "Pole location:\n"
        "  Real part: decay/growth rate\n"
        "  Imag part: oscillation frequency\n"
        "  Distance from origin: speed of response\n"
        "  Angle from neg. real axis: damping ratio\n\n"
        "Root Locus: Poles move as gain K varies.\n"
        "Start at OL poles (K=0), end at OL zeros (K->inf).\n"
        "Unstable when a locus crosses into the right half-plane.");
    ImGui::SameLine();

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
    if (state.pz_mode == PZMode::PlantLocus) {
        if (ImGui::BeginTable("##rl_ctrls", 3)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##K", &state.rl_current_k,
                               state.rl_k_min, state.rl_k_max, "%.2f",
                               ImGuiSliderFlags_Logarithmic);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K min"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmin", &state.rl_k_min, 0.0f, 1.0f))
                state.needs_recompute = true;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmax", &state.rl_k_max, 1.0f, 1000.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    } else if (state.pz_mode == PZMode::StateFB) {
        if (ImGui::BeginTable("##sfb_ctrls", 2)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##alpha", &state.rl_current_alpha,
                               state.rl_alpha_min, state.rl_alpha_max);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1 max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##alphamax", &state.rl_alpha_max, 0.1f, 10.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    }

    // Compute axis limits from visible data (mode-dependent)
    double pz_xmin = 0, pz_xmax = 0, pz_ymin = 0, pz_ymax = 0;
    bool pz_has_data = false;
    auto pz_extend = [&](double re, double im) {
        if (!pz_has_data) { pz_xmin = pz_xmax = re; pz_ymin = pz_ymax = im; pz_has_data = true; }
        else { pz_xmin = std::min(pz_xmin, re); pz_xmax = std::max(pz_xmax, re);
               pz_ymin = std::min(pz_ymin, im); pz_ymax = std::max(pz_ymax, im); }
    };
    if (state.pz_mode == PZMode::PoleZero) {
        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            for (const auto& p : state.pole_zero[s].poles) pz_extend(p.real(), p.imag());
            for (const auto& z : state.pole_zero[s].zeros) pz_extend(z.real(), z.imag());
        }
    } else {
        for (const auto& step : state.root_locus)
            for (const auto& p : step.poles) pz_extend(p.real(), p.imag());
        if (state.system_valid[0])
            for (const auto& p : state.pole_zero[0].poles) pz_extend(p.real(), p.imag());
    }
    if (!pz_has_data) { pz_xmin = -2; pz_xmax = 2; pz_ymin = -2; pz_ymax = 2; }
    {
        double mx = std::max((pz_xmax - pz_xmin) * 0.15, 0.5);
        double my = std::max((pz_ymax - pz_ymin) * 0.15, 0.5);
        pz_xmin -= mx; pz_xmax += mx; pz_ymin -= my; pz_ymax += my;
        pz_xmin = std::min(pz_xmin, -1.0);
        pz_xmax = std::max(pz_xmax, 1.0);
    }

    // Plot
    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (ImPlot::BeginPlot("##pz_plot", avail, ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("\xcf\x83 (Real)", "j\xcf\x89 (Imag)",
                          ImPlotAxisFlags_NoGridLines,
                          ImPlotAxisFlags_NoGridLines);
        ImPlot::SetupAxisLimits(ImAxis_X1, pz_xmin, pz_xmax, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, pz_ymin, pz_ymax, ImPlotCond_Always);

        // Stable region shading (Re < 0)
        {
            auto limits = ImPlot::GetPlotLimits();
            double shade_x[4] = {limits.X.Min, 0, 0, limits.X.Min};
            double shade_y[4] = {limits.Y.Min, limits.Y.Min, limits.Y.Max, limits.Y.Max};
            ImPlot::PlotShaded("##stable", shade_x, shade_y, 4, 0.0,
                               ImPlotSpec(ImPlotProp_FillColor,
                                          ImVec4(0.20f, 0.78f, 0.20f, 0.08f)));
        }

        // Stability boundary (imaginary axis)
        double zero = 0.0;
        ImPlot::PlotInfLines("##jw_axis", &zero, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f)));

        // Real axis
        ImPlot::PlotInfLines("##sigma_axis", &zero, 1,
                             ImPlotSpec(ImPlotProp_LineColor,
                                        ImVec4(1, 1, 1, 0.16f),
                                        ImPlotProp_Flags,
                                        (int)ImPlotInfLinesFlags_Horizontal));

        // S-plane grid (wn circles + zeta lines)
        drawSPlaneGrid();

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
                    char label_p[32];
                    std::snprintf(label_p, sizeof(label_p), "P(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_p, pre.data(), pim.data(),
                                       static_cast<int>(pre.size()),
                                       ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Cross,
                                                  ImPlotProp_MarkerSize, 8.0f,
                                                  ImPlotProp_MarkerLineColor, system_colors[s],
                                                  ImPlotProp_LineWeight, 2.0f));
                }

                std::vector<double> zre, zim;
                for (const auto& z : pz.zeros) {
                    zre.push_back(z.real());
                    zim.push_back(z.imag());
                }
                if (!zre.empty()) {
                    char label_z[32];
                    std::snprintf(label_z, sizeof(label_z), "Z(%s)", system_names[s]);
                    ImPlot::PlotScatter(label_z, zre.data(), zim.data(),
                                       static_cast<int>(zre.size()),
                                       ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Circle,
                                                  ImPlotProp_MarkerSize, 8.0f,
                                                  ImPlotProp_MarkerFillColor, ImVec4(0, 0, 0, 0),
                                                  ImPlotProp_MarkerLineColor, system_colors[s],
                                                  ImPlotProp_LineWeight, 2.0f));
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
                    char lbl[16];
                    std::snprintf(lbl, sizeof(lbl), "##rl%d", p);
                    ImPlot::PlotLine(lbl, re.data(), im.data(), n_steps,
                                     ImPlotSpec(ImPlotProp_LineColor, system_colors[0],
                                                ImPlotProp_LineWeight, 1.5f));
                }

                double current_gain = (state.pz_mode == PZMode::PlantLocus)
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
                ImPlot::PlotScatter("##current", cre.data(), cim.data(),
                                   static_cast<int>(cre.size()),
                                   ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Diamond,
                                              ImPlotProp_MarkerSize, 10.0f,
                                              ImPlotProp_MarkerLineColor, system_colors[3],
                                              ImPlotProp_LineWeight, 2.5f));
            }

            const auto& pz = state.pole_zero[0];
            std::vector<double> pre, pim;
            for (const auto& p : pz.poles) {
                pre.push_back(p.real());
                pim.push_back(p.imag());
            }
            ImPlot::PlotScatter("OL Poles", pre.data(), pim.data(),
                               static_cast<int>(pre.size()),
                               ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Cross,
                                          ImPlotProp_MarkerSize, 6.0f,
                                          ImPlotProp_MarkerLineColor, ImVec4(0.5f, 0.5f, 0.5f, 0.7f),
                                          ImPlotProp_LineWeight, 1.5f));
        }

        ImPlot::EndPlot();
    }

    // A suppressed trace must be distinguishable from breakage.
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (s == 2 || state.channel_reason[s].empty()) continue;
        ImGui::TextDisabled("%s: %s", system_names[s],
                            state.channel_reason[s].c_str());
    }

    ImGui::End();
}

}  // namespace caliburn
