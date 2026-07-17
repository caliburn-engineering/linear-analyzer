// src/panels/nyquist_panel.cpp
#include "nyquist_panel.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

void drawNyquistPanel(AppState& state) {
    if (!state.show_nyquist) return;
    ImGui::Begin("Nyquist Plot", &state.show_nyquist);

    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (ImPlot::BeginPlot("##nyquist", avail, ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("Re{G(j\xcf\x89)}", "Im{G(j\xcf\x89)}");

        // Critical point (-1, 0)
        {
            double neg1 = -1.0, zero = 0.0;
            ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 8,
                                       ImVec4(1, 0.2f, 0.2f, 1), 2.0f,
                                       ImVec4(0, 0, 0, 0));
            ImPlot::PlotScatter("(-1,0)", &neg1, &zero, 1);
        }

        for (int s = 0; s < NUM_SYSTEMS; ++s) {
            if (!state.trace_visible[s] || !state.system_valid[s]) continue;
            const auto& bode = state.bode[s];
            if (bode.nyquist.empty()) continue;

            int n = static_cast<int>(bode.nyquist.size());
            std::vector<double> re(n), im(n);
            for (int k = 0; k < n; ++k) {
                re[k] = bode.nyquist[k].real();
                im[k] = bode.nyquist[k].imag();
            }

            ImPlot::SetNextLineStyle(system_colors[s], 2.0f);
            char lbl[32];
            std::snprintf(lbl, sizeof(lbl), "%s+", system_names[s]);
            ImPlot::PlotLine(lbl, re.data(), im.data(), n);

            std::vector<double> im_neg(n);
            for (int k = 0; k < n; ++k) im_neg[k] = -im[k];
            ImVec4 col_dim = system_colors[s];
            col_dim.w = 0.5f;
            ImPlot::SetNextLineStyle(col_dim, 1.5f);
            std::snprintf(lbl, sizeof(lbl), "%s-", system_names[s]);
            ImPlot::PlotLine(lbl, re.data(), im_neg.data(), n);

            for (int k = 0; k < n; k += std::max(n / 10, 1)) {
                ImPlot::SetNextMarkerStyle(ImPlotMarker_Right, 4,
                                           system_colors[s], 1.0f);
                ImPlot::PlotScatter("##arr", &re[k], &im[k], 1);
            }
        }

        // Frequency tooltip on hover
        if (ImPlot::IsPlotHovered()) {
            ImPlotPoint mouse = ImPlot::GetPlotMousePos();
            for (int s = 0; s < NUM_SYSTEMS; ++s) {
                if (!state.trace_visible[s] || !state.system_valid[s]) continue;
                const auto& bode = state.bode[s];
                if (bode.nyquist.empty()) continue;
                double best_d = 1e30;
                int best_k = 0;
                int n = static_cast<int>(bode.nyquist.size());
                for (int k = 0; k < n; ++k) {
                    double dx = bode.nyquist[k].real() - mouse.x;
                    double dy = bode.nyquist[k].imag() - mouse.y;
                    double d = dx * dx + dy * dy;
                    if (d < best_d) { best_d = d; best_k = k; }
                }
                if (best_d < 0.5) {
                    ImPlot::Annotation(
                        bode.nyquist[best_k].real(),
                        bode.nyquist[best_k].imag(),
                        system_colors[s], ImVec2(10, -10), true,
                        "f=%.3f Hz", bode.points[best_k].freq_hz);
                }
                break;
            }
        }

        ImPlot::EndPlot();
    }

    ImGui::End();
}

}  // namespace caliburn
