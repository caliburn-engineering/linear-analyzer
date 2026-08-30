// src/panels/nyquist_panel.cpp
#include "nyquist_panel.h"
#include "panel_utils.h"
#include "implot.h"
#include <cmath>

namespace caliburn {

void drawNyquistPanel(AppState& state) {
    if (!state.show_nyquist) return;
    ImGui::Begin("Nyquist Plot", &state.show_nyquist);
    drawHelpMarker(
        "Nyquist Plot\n\n"
        "Plots G(jw) in the complex plane as w varies 0 to inf.\n\n"
        "Nyquist Stability Criterion:\n"
        "CL unstable poles = OL unstable poles + clockwise\n"
        "encirclements of the critical point (-1, 0).\n\n"
        "Features:\n"
        "  (-1, 0): Red circle. Passing through = marginal stability.\n"
        "  Solid curve: Positive frequencies (w > 0).\n"
        "  Faded curve: Mirror for negative frequencies.\n"
        "  Arrows: Direction of increasing frequency.\n"
        "  Hover: Shows frequency at that curve point.\n\n"
        "Distance from (-1,0) relates to gain/phase margins.");

    // Compute axis limits from Nyquist data + critical point
    double nq_xmin = -1.0, nq_xmax = 0.0, nq_ymin = 0.0, nq_ymax = 0.0;
    bool nq_has = false;
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (!state.trace_visible[s] || !state.system_valid[s]) continue;
        for (const auto& c : state.bode[s].nyquist) {
            double re = c.real(), im = c.imag();
            if (!nq_has) { nq_xmin = nq_xmax = re; nq_ymin = std::min(im, -im); nq_ymax = std::max(im, -im); nq_has = true; }
            else { nq_xmin = std::min(nq_xmin, re); nq_xmax = std::max(nq_xmax, re);
                   nq_ymin = std::min({nq_ymin, im, -im}); nq_ymax = std::max({nq_ymax, im, -im}); }
        }
    }
    nq_xmin = std::min(nq_xmin, -1.0);  // always include critical point
    if (!nq_has) { nq_xmin = -2; nq_xmax = 2; nq_ymin = -2; nq_ymax = 2; }
    {
        double mx = std::max((nq_xmax - nq_xmin) * 0.1, 0.2);
        double my = std::max((nq_ymax - nq_ymin) * 0.1, 0.2);
        nq_xmin -= mx; nq_xmax += mx; nq_ymin -= my; nq_ymax += my;
        nq_xmin = std::min(nq_xmin, -1.0);
        nq_xmax = std::max(nq_xmax, 1.0);
    }

    ImVec2 avail = ImGui::GetContentRegionAvail();
    if (ImPlot::BeginPlot("##nyquist", avail, ImPlotFlags_Equal)) {
        ImPlot::SetupAxes("Re{G(j\xcf\x89)}", "Im{G(j\xcf\x89)}",
                          ImPlotAxisFlags_NoGridLines,
                          ImPlotAxisFlags_NoGridLines);
        ImPlot::SetupAxisLimits(ImAxis_X1, nq_xmin, nq_xmax, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, nq_ymin, nq_ymax, ImPlotCond_Always);

        // Polar grid (concentric circles + radial lines)
        drawPolarGrid();

        // Critical point (-1, 0)
        {
            double neg1 = -1.0, zero = 0.0;
            ImPlot::PlotScatter("(-1,0)", &neg1, &zero, 1,
                                ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Circle,
                                           ImPlotProp_MarkerSize, 8.0f,
                                           ImPlotProp_MarkerFillColor, ImVec4(0, 0, 0, 0),
                                           ImPlotProp_MarkerLineColor, ImVec4(1, 0.2f, 0.2f, 1),
                                           ImPlotProp_LineWeight, 2.0f));
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

            char lbl[32];
            std::snprintf(lbl, sizeof(lbl), "%s+", system_names[s]);
            ImPlot::PlotLine(lbl, re.data(), im.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, system_colors[s],
                                        ImPlotProp_LineWeight, 2.0f));

            std::vector<double> im_neg(n);
            for (int k = 0; k < n; ++k) im_neg[k] = -im[k];
            ImVec4 col_dim = system_colors[s];
            col_dim.w = 0.5f;
            std::snprintf(lbl, sizeof(lbl), "%s-", system_names[s]);
            ImPlot::PlotLine(lbl, re.data(), im_neg.data(), n,
                             ImPlotSpec(ImPlotProp_LineColor, col_dim,
                                        ImPlotProp_LineWeight, 1.5f));

            // Direction arrows along curve
            int step = std::max(n / 10, 1);
            ImDrawList* dl = ImPlot::GetPlotDrawList();
            ImU32 arrow_col = system_colors_u32[s];
            for (int k = step / 2; k < n - 1; k += step) {
                // Use points a few indices apart for a stable direction
                int ka = std::max(k - 2, 0);
                int kb = std::min(k + 2, n - 1);
                ImVec2 pa = ImPlot::PlotToPixels(ImPlotPoint(re[ka], im[ka]));
                ImVec2 pb = ImPlot::PlotToPixels(ImPlotPoint(re[kb], im[kb]));
                float pdx = pb.x - pa.x;
                float pdy = pb.y - pa.y;
                float plen = std::sqrt(pdx * pdx + pdy * pdy);
                if (plen < 4.0f) continue;

                // Unit direction and perpendicular in pixel space
                float ux = pdx / plen, uy = pdy / plen;
                ImVec2 tip = ImPlot::PlotToPixels(ImPlotPoint(re[k], im[k]));
                float sz = 7.0f;
                ImVec2 a1(tip.x - ux * sz * 2.5f + uy * sz,
                          tip.y - uy * sz * 2.5f - ux * sz);
                ImVec2 a2(tip.x - ux * sz * 2.5f - uy * sz,
                          tip.y - uy * sz * 2.5f + ux * sz);
                dl->AddTriangleFilled(tip, a1, a2, arrow_col);
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
