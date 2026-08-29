// src/panels/loop_ui_prototype.cpp
//
// PROTOTYPE — THROWAWAY. Wayfinder ticket #6. See loop_ui_prototype.h.
//
// Three variants of the controller section, switchable from the floating bar:
//   A — Ledger      : one table row per loop, gains inline as columns.
//   B — Master/detail: compact loop list + full-width gain block for one loop.
//   C — Pairing grid : the p x m diagnostic matrix IS the pairing editor.
//
// The diagnostic numbers are real (evalTransferFunction on the live plant), so
// Ball-Balancer's dead identity pairing shows up as it actually would.

#include "loop_ui_prototype.h"
#include "../analysis/frequency_response.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace caliburn::proto {
namespace {

// ---------------------------------------------------------------- proto state

enum class Kind { None, PID, LeadLag };
enum class LLMode { Lead, Lag, LeadLag };

struct PID { float Kp = 1.0f, Ki = 0.0f, Kd = 0.0f, tau_f = 0.01f; };
struct LL  { LLMode mode = LLMode::Lead; float Kc = 1.0f, alpha_lead = 0.1f,
             T_lead = 1.0f, alpha_lag = 10.0f, T_lag = 1.0f; };

struct Loop { int out = 0, in = 0; PID pid; LL ll; };

struct Proto {
    int variant = 0;                 // 0=A 1=B 2=C
    Kind kind = Kind::PID;
    std::vector<Loop> loops;
    int selected = 0;                // variants B and C
    float omega_hz = 1.0f;
    bool sweep = false;              // variant C: show RGA-number-vs-omega strip
    std::vector<float> scales;       // per output, default 1
    int ki_scale = 2;                // which bench scale the variants use
    bool bench_open = false;
    float bench[4] = {0.5f, 0.5f, 0.5f, 0.5f};
    int seeded_p = -1, seeded_m = -1;
};
Proto g;

const char* kVariantNames[3] = {"Ledger", "Master / detail", "Pairing grid"};

// ------------------------------------------------------------------- reseeding

// #2: identity y[k] <- u[k], reseeded wholesale by one (p, m) dimension guard.
void reseedIfNeeded(const AppState& s) {
    int p = s.plant.outputs(), m = s.plant.inputs();
    if (p == g.seeded_p && m == g.seeded_m) return;
    g.seeded_p = p; g.seeded_m = m;
    g.loops.clear();
    for (int k = 0; k < std::min(p, m); ++k) { Loop l; l.out = k; l.in = k; g.loops.push_back(l); }
    g.scales.assign(static_cast<size_t>(p), 1.0f);
    g.selected = 0;
}

bool scalesAtDefault() {
    for (float v : g.scales) if (std::abs(v - 1.0f) > 1e-6f) return false;
    return true;
}

// ----------------------------------------------------------------- diagnostics

Eigen::MatrixXcd Gjw(const LinearSystem& sys, double f_hz) {
    std::complex<double> s(0.0, 2.0 * 3.14159265358979323846 * f_hz);
    Eigen::MatrixXcd G(sys.outputs(), sys.inputs());
    for (int i = 0; i < sys.outputs(); ++i)
        for (int j = 0; j < sys.inputs(); ++j)
            G(i, j) = evalTransferFunction(sys, i, j, s);
    return G;
}

struct Diag {
    bool has_rga = false;
    Eigen::MatrixXd lambda;              // sub-matrix RGA, real part
    std::vector<int> sub_out, sub_in;    // sub-matrix index maps
    double rga_number = 0.0;

    bool has_share = false;
    int share_in = 0;
    std::vector<double> share, mag_db;   // per plant output

    std::vector<double> loop_lambda;     // per loop, NaN when no RGA
    std::vector<bool>   loop_dead;       // per loop
    std::string headline;
};

// #8: a paired channel is structurally dead when |g| is below tol at EVERY omega.
bool deadChannel(const LinearSystem& sys, int i, int j) {
    for (int k = 0; k < 60; ++k) {
        double f = 0.01 * std::pow(10.0, 4.0 * k / 59.0);
        if (std::abs(evalTransferFunction(sys, i, j, {0.0, 2.0 * 3.14159265358979323846 * f})) > 1e-9)
            return false;
    }
    return true;
}

Diag computeDiag(const AppState& s) {
    Diag d;
    const auto& G = Gjw(s.plant, g.omega_hz);
    int p = s.plant.outputs();

    std::vector<int> outs, ins;
    for (const auto& l : g.loops) {
        if (std::find(outs.begin(), outs.end(), l.out) == outs.end()) outs.push_back(l.out);
        if (std::find(ins.begin(), ins.end(), l.in) == ins.end()) ins.push_back(l.in);
    }
    std::sort(outs.begin(), outs.end());
    std::sort(ins.begin(), ins.end());

    d.loop_lambda.assign(g.loops.size(), std::nan(""));
    d.loop_dead.resize(g.loops.size());
    for (size_t k = 0; k < g.loops.size(); ++k)
        d.loop_dead[k] = deadChannel(s.plant, g.loops[k].out, g.loops[k].in);

    // #4/#8: RGA exists only on a square paired sub-matrix of size >= 2.
    if (outs.size() == ins.size() && outs.size() >= 2) {
        int q = static_cast<int>(outs.size());
        Eigen::MatrixXcd sub(q, q);
        for (int r = 0; r < q; ++r)
            for (int c = 0; c < q; ++c) sub(r, c) = G(outs[r], ins[c]);
        Eigen::MatrixXcd L = sub.array() * sub.inverse().transpose().array();  // never adjoint()
        d.has_rga = true;
        d.lambda = L.real();
        d.sub_out = outs; d.sub_in = ins;

        Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(q, q);
        for (const auto& l : g.loops) {
            int r = static_cast<int>(std::find(outs.begin(), outs.end(), l.out) - outs.begin());
            int c = static_cast<int>(std::find(ins.begin(), ins.end(), l.in) - ins.begin());
            Pi(r, c) = 1.0;
        }
        d.rga_number = (d.lambda - Pi).cwiseAbs().sum();

        for (size_t k = 0; k < g.loops.size(); ++k) {
            int r = static_cast<int>(std::find(outs.begin(), outs.end(), g.loops[k].out) - outs.begin());
            int c = static_cast<int>(std::find(ins.begin(), ins.end(), g.loops[k].in) - ins.begin());
            d.loop_lambda[k] = d.lambda(r, c);
        }
        char buf[160];
        std::snprintf(buf, sizeof buf, "RGA @ %.3g Hz  -  RGA number %.2f", g.omega_hz, d.rga_number);
        d.headline = buf;
    } else if (!g.loops.empty()) {
        // #8: Channel Share. Not a pairing measure. Scale-dependent.
        d.has_share = true;
        d.share_in = g.loops.front().in;
        double tot = 0.0;
        d.share.assign(p, 0.0); d.mag_db.assign(p, 0.0);
        for (int i = 0; i < p; ++i) {
            double gi = std::abs(G(i, d.share_in)) / std::max(1e-12f, g.scales[i]);
            d.share[i] = gi * gi;
            tot += d.share[i];
            d.mag_db[i] = 20.0 * std::log10(std::max(std::abs(G(i, d.share_in)), 1e-300));
        }
        for (int i = 0; i < p; ++i) d.share[i] = tot > 0 ? d.share[i] / tot : 0.0;
        char buf[160];
        std::snprintf(buf, sizeof buf, "Channel Share of u%d @ %.3g Hz  -  scales %s",
                      d.share_in, g.omega_hz, scalesAtDefault() ? "at default (all 1)" : "user-set");
        d.headline = buf;
    } else {
        d.headline = "No loops - controller is None";
    }
    return d;
}

double rgaNumberAt(const AppState& s, double f_hz) {
    double saved = g.omega_hz;
    g.omega_hz = static_cast<float>(f_hz);
    Diag d = computeDiag(s);
    g.omega_hz = static_cast<float>(saved);
    return d.has_rga ? d.rga_number : std::nan("");
}

// ------------------------------------------------------------------- widgets

ImVec4 severityColor(double rga_number) {
    if (rga_number < 0.5) return ImVec4(0.30f, 1.00f, 0.30f, 1.0f);
    if (rga_number < 2.0) return ImVec4(1.00f, 0.80f, 0.20f, 1.0f);
    return ImVec4(1.00f, 0.35f, 0.35f, 1.0f);
}

// Candidate scales for Ki / Kd, which must reach 0 and still resolve 0.001.
bool gainSlider(const char* id, float* v, float hi, int mode) {
    ImGui::SetNextItemWidth(-FLT_MIN);
    switch (mode) {
        case 0:  // plain linear
            return ImGui::SliderFloat(id, v, 0.0f, hi, "%.4g");
        case 1:  // ImGui's own logarithmic flag, range including zero
            return ImGui::SliderFloat(id, v, 0.0f, hi, "%.4g", ImGuiSliderFlags_Logarithmic);
        case 2: {  // hand-rolled piecewise: linear 0-1 over first half, log 1-hi over second
            const float lin_hi = 1.0f;
            float t = (*v <= lin_hi) ? 0.5f * (*v / lin_hi)
                                     : 0.5f + 0.5f * std::log(*v / lin_hi) / std::log(hi / lin_hi);
            char fmt[32]; std::snprintf(fmt, sizeof fmt, "%.4g", *v);
            if (ImGui::SliderFloat(id, &t, 0.0f, 1.0f, fmt)) {
                *v = (t <= 0.5f) ? 2.0f * t * lin_hi
                                 : lin_hi * std::pow(hi / lin_hi, (t - 0.5f) * 2.0f);
                return true;
            }
            return false;
        }
        default: {  // drag, speed proportional to current magnitude
            float speed = std::max(0.0005f, std::abs(*v) * 0.02f);
            bool c = ImGui::DragFloat(id, v, speed, 0.0f, hi, "%.4g");
            *v = std::clamp(*v, 0.0f, hi);
            return c;
        }
    }
}

void pairCombos(AppState& s, Loop& l, float width) {
    int p = s.plant.outputs(), m = s.plant.inputs();
    std::vector<std::string> ys, us;
    for (int i = 0; i < p; ++i) ys.push_back("y" + std::to_string(i));
    for (int j = 0; j < m; ++j) us.push_back("u" + std::to_string(j));
    auto getter = [](void* data, int idx) -> const char* {
        return (*static_cast<std::vector<std::string>*>(data))[idx].c_str();
    };
    ImGui::SetNextItemWidth(width);
    ImGui::BeginDisabled(p == 1);
    ImGui::Combo("##out", &l.out, getter, &ys, p);
    ImGui::EndDisabled();
    ImGui::SameLine(0, 4);
    ImGui::TextUnformatted("<-");
    ImGui::SameLine(0, 4);
    ImGui::SetNextItemWidth(width);
    ImGui::BeginDisabled(m == 1);
    ImGui::Combo("##in", &l.in, getter, &us, m);
    ImGui::EndDisabled();
}

// Full-width labelled gain block (variants B and C).
void gainsBlock(Loop& l) {
    if (g.kind == Kind::PID) {
        ImGui::TextUnformatted("Kp"); ImGui::SameLine(34);
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::SliderFloat("##Kp", &l.pid.Kp, 0.01f, 100.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
        ImGui::TextUnformatted("Ki"); ImGui::SameLine(34); gainSlider("##Ki", &l.pid.Ki, 50.0f, g.ki_scale);
        ImGui::TextUnformatted("Kd"); ImGui::SameLine(34); gainSlider("##Kd", &l.pid.Kd, 20.0f, g.ki_scale);
        ImGui::TextUnformatted("tf"); ImGui::SameLine(34);
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::SliderFloat("##tf", &l.pid.tau_f, 0.001f, 1.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
        char f[160];
        std::snprintf(f, sizeof f, "C(s) = %.3g + %.3g/s + %.3g s/(%.3g s + 1)",
                      l.pid.Kp, l.pid.Ki, l.pid.Kd, l.pid.tau_f);
        ImGui::TextDisabled("%s", f);
        int n = (l.pid.Ki != 0.0f ? 1 : 0) + (l.pid.Kd != 0.0f ? 1 : 0);
        ImGui::TextDisabled("states: %d", n);
    } else {
        const char* modes[] = {"Lead", "Lag", "Lead-Lag"};
        int mi = static_cast<int>(l.ll.mode);
        ImGui::SetNextItemWidth(140);
        if (ImGui::Combo("Mode", &mi, modes, 3)) l.ll.mode = static_cast<LLMode>(mi);
        ImGui::TextUnformatted("Kc"); ImGui::SameLine(34);
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::SliderFloat("##Kc", &l.ll.Kc, 0.1f, 100.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
        if (l.ll.mode != LLMode::Lag) {
            ImGui::TextDisabled("lead section");
            ImGui::TextUnformatted("a"); ImGui::SameLine(34);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##al", &l.ll.alpha_lead, 0.01f, 0.99f, "%.4g", ImGuiSliderFlags_Logarithmic);
            ImGui::TextUnformatted("T"); ImGui::SameLine(34);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##Tl", &l.ll.T_lead, 0.001f, 100.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
        }
        if (l.ll.mode != LLMode::Lead) {
            ImGui::TextDisabled("lag section");
            ImGui::TextUnformatted("a"); ImGui::SameLine(34);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##ag", &l.ll.alpha_lag, 1.01f, 100.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
            ImGui::TextUnformatted("T"); ImGui::SameLine(34);
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##Tg", &l.ll.T_lag, 0.001f, 100.0f, "%.4g", ImGuiSliderFlags_Logarithmic);
        }
        ImGui::TextDisabled("states: %d", l.ll.mode == LLMode::LeadLag ? 2 : 1);
    }
}

void kindSelector(AppState& state) {
    const char* kinds[] = {"None", "PID", "Lead/Lag"};
    int ki = static_cast<int>(g.kind);
    ImGui::SetNextItemWidth(160);
    if (ImGui::Combo("Type", &ki, kinds, 3)) {
        g.kind = static_cast<Kind>(ki);
        if (g.kind == Kind::None) g.loops.clear();
        else if (g.loops.empty()) { g.seeded_p = -1; reseedIfNeeded(state); }
    }
    ImGui::SameLine();
    ImGui::TextDisabled("(%d loop%s)", (int)g.loops.size(), g.loops.size() == 1 ? "" : "s");
}

void omegaControl(float width) {
    ImGui::SetNextItemWidth(width);
    ImGui::SliderFloat("##omega", &g.omega_hz, 0.01f, 100.0f, "w = %.3g Hz",
                       ImGuiSliderFlags_Logarithmic);
}

void addRemoveButtons(AppState& s, int idx_for_remove) {
    if (ImGui::SmallButton("+ Add loop")) {
        Loop l;
        l.out = g.loops.empty() ? 0 : std::min(s.plant.outputs() - 1, g.loops.back().out + 1);
        l.in  = g.loops.empty() ? 0 : std::min(s.plant.inputs() - 1, g.loops.back().in + 1);
        g.loops.push_back(l);
    }
    if (idx_for_remove >= 0 && idx_for_remove < (int)g.loops.size()) {
        ImGui::SameLine();
        if (ImGui::SmallButton("Remove")) {
            g.loops.erase(g.loops.begin() + idx_for_remove);
            g.selected = std::max(0, g.selected - 1);
        }
        ImGui::SameLine();
        if (ImGui::SmallButton("Up") && idx_for_remove > 0) {
            std::swap(g.loops[idx_for_remove], g.loops[idx_for_remove - 1]);
            g.selected = idx_for_remove - 1;
        }
        ImGui::SameLine();
        if (ImGui::SmallButton("Down") && idx_for_remove + 1 < (int)g.loops.size()) {
            std::swap(g.loops[idx_for_remove], g.loops[idx_for_remove + 1]);
            g.selected = idx_for_remove + 1;
        }
    }
}

// Duplicate pairings are flagged, never blocked (#2).
int duplicateCount() {
    int dup = 0;
    for (size_t a = 0; a < g.loops.size(); ++a)
        for (size_t b = a + 1; b < g.loops.size(); ++b)
            if (g.loops[a].out == g.loops[b].out && g.loops[a].in == g.loops[b].in) ++dup;
    return dup;
}

void scalesRow(AppState& s) {
    int p = s.plant.outputs();
    char hdr[96];
    std::snprintf(hdr, sizeof hdr, "Output scales: %s###scales",
                  scalesAtDefault() ? "default (all 1)" : "user-set");
    if (ImGui::TreeNode(hdr)) {
        ImGui::TextDisabled("max allowed deviation per output - engineering intent, not derivable");
        for (int i = 0; i < p; ++i) {
            ImGui::PushID(i);
            char lbl[16]; std::snprintf(lbl, sizeof lbl, "y%d", i);
            ImGui::SetNextItemWidth(120);
            ImGui::DragFloat(lbl, &g.scales[i], 0.01f, 0.001f, 1000.0f, "%.3g");
            ImGui::PopID();
        }
        ImGui::TreePop();
    }
}

void deadWarnings(const Diag& d) {
    for (size_t k = 0; k < g.loops.size(); ++k) {
        if (!d.loop_dead[k]) continue;
        ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1),
                           "loop %d  y%d <- u%d is structurally dead: |g| below tolerance at every w",
                           (int)k, g.loops[k].out, g.loops[k].in);
    }
    int dup = duplicateCount();
    if (dup > 0)
        ImGui::TextColored(ImVec4(1, 0.8f, 0.2f, 1),
                           "%d duplicate pairing%s - gains sum, which is legal", dup, dup == 1 ? "" : "s");
}

// ---------------------------------------------------------------- variant A

void variantA(AppState& s, const Diag& d) {
    kindSelector(s);
    if (g.kind == Kind::None) return;

    addRemoveButtons(s, -1);
    ImGui::SameLine();
    omegaControl(-FLT_MIN);

    int ncol = (g.kind == Kind::PID) ? 8 : 8;
    ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                            ImGuiTableFlags_SizingStretchProp;
    if (ImGui::BeginTable("##loops", ncol, flags)) {
        ImGui::TableSetupColumn("#", ImGuiTableColumnFlags_WidthFixed, 18);
        ImGui::TableSetupColumn("pair", ImGuiTableColumnFlags_WidthFixed, 130);
        if (g.kind == Kind::PID) {
            ImGui::TableSetupColumn("Kp"); ImGui::TableSetupColumn("Ki");
            ImGui::TableSetupColumn("Kd"); ImGui::TableSetupColumn("tf");
        } else {
            ImGui::TableSetupColumn("mode", ImGuiTableColumnFlags_WidthFixed, 78);
            ImGui::TableSetupColumn("Kc"); ImGui::TableSetupColumn("a"); ImGui::TableSetupColumn("T");
        }
        ImGui::TableSetupColumn(d.has_rga ? "lambda" : "share", ImGuiTableColumnFlags_WidthFixed, 50);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 22);
        ImGui::TableHeadersRow();

        int to_erase = -1;
        for (size_t k = 0; k < g.loops.size(); ++k) {
            ImGui::PushID((int)k);
            auto& l = g.loops[k];
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::Text("%d", (int)k);
            ImGui::TableNextColumn(); pairCombos(s, l, 46);
            if (g.kind == Kind::PID) {
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##Kp", &l.pid.Kp, 0.05f, 0.01f, 100.0f, "%.3g");
                ImGui::TableNextColumn(); gainSlider("##Ki", &l.pid.Ki, 50.0f, g.ki_scale);
                ImGui::TableNextColumn(); gainSlider("##Kd", &l.pid.Kd, 20.0f, g.ki_scale);
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##tf", &l.pid.tau_f, 0.001f, 0.001f, 1.0f, "%.3g");
            } else {
                ImGui::TableNextColumn();
                const char* modes[] = {"Lead", "Lag", "L-L"};
                int mi = static_cast<int>(l.ll.mode); ImGui::SetNextItemWidth(-FLT_MIN);
                if (ImGui::Combo("##mode", &mi, modes, 3)) l.ll.mode = static_cast<LLMode>(mi);
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##Kc", &l.ll.Kc, 0.05f, 0.1f, 100.0f, "%.3g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##a", l.ll.mode == LLMode::Lag ? &l.ll.alpha_lag : &l.ll.alpha_lead,
                                 0.01f, 0.01f, 100.0f, "%.3g");
                ImGui::TableNextColumn(); ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::DragFloat("##T", l.ll.mode == LLMode::Lag ? &l.ll.T_lag : &l.ll.T_lead,
                                 0.05f, 0.001f, 100.0f, "%.3g");
            }
            ImGui::TableNextColumn();
            if (d.has_rga) {
                double lam = d.loop_lambda[k];
                ImGui::TextColored(severityColor(d.rga_number), "%.2f", lam);
            } else if (d.has_share) {
                ImGui::Text("%.0f%%", 100.0 * d.share[l.out]);
            } else ImGui::TextDisabled("-");
            ImGui::TableNextColumn();
            if (ImGui::SmallButton("x")) to_erase = (int)k;
            ImGui::PopID();
        }
        ImGui::EndTable();
        if (to_erase >= 0) g.loops.erase(g.loops.begin() + to_erase);
    }

    ImGui::TextColored(d.has_rga ? severityColor(d.rga_number) : ImVec4(0.7f, 0.7f, 0.7f, 1),
                       "%s", d.headline.c_str());
    if (d.has_share) {
        for (int i = 0; i < (int)d.share.size(); ++i)
            ImGui::TextDisabled("  y%d  %5.1f%%   |g| %7.1f dB", i, 100.0 * d.share[i], d.mag_db[i]);
    }
    deadWarnings(d);
    scalesRow(s);
}

// ---------------------------------------------------------------- variant B

void variantB(AppState& s, const Diag& d) {
    kindSelector(s);
    if (g.kind == Kind::None) return;

    ImGui::TextColored(d.has_rga ? severityColor(d.rga_number) : ImVec4(0.7f, 0.7f, 0.7f, 1),
                       "%s", d.headline.c_str());
    omegaControl(-FLT_MIN);
    deadWarnings(d);

    ImGui::BeginChild("##looplist", ImVec2(0, std::min(6, (int)g.loops.size() + 1) * 20.0f + 8),
                      ImGuiChildFlags_Borders);
    for (size_t k = 0; k < g.loops.size(); ++k) {
        const auto& l = g.loops[k];
        char row[192];
        if (g.kind == Kind::PID)
            std::snprintf(row, sizeof row, "y%d <- u%d   Kp %.3g  Ki %.3g  Kd %.3g",
                          l.out, l.in, l.pid.Kp, l.pid.Ki, l.pid.Kd);
        else
            std::snprintf(row, sizeof row, "y%d <- u%d   %s  Kc %.3g", l.out, l.in,
                          l.ll.mode == LLMode::Lead ? "Lead" : l.ll.mode == LLMode::Lag ? "Lag" : "Lead-Lag",
                          l.ll.Kc);
        ImGui::PushID((int)k);
        if (ImGui::Selectable(row, g.selected == (int)k)) g.selected = (int)k;
        if (d.loop_dead[k]) { ImGui::SameLine(); ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1), " dead"); }
        else if (d.has_rga) {
            ImGui::SameLine();
            ImGui::TextColored(severityColor(std::abs(d.loop_lambda[k] - 1.0)), " l=%.2f", d.loop_lambda[k]);
        }
        ImGui::PopID();
    }
    ImGui::EndChild();

    addRemoveButtons(s, g.selected);

    if (g.selected >= 0 && g.selected < (int)g.loops.size()) {
        auto& l = g.loops[g.selected];
        ImGui::SeparatorText("Selected loop");
        pairCombos(s, l, 60);
        gainsBlock(l);
    }
    scalesRow(s);
}

// ---------------------------------------------------------------- variant C

void variantC(AppState& s, const Diag& d) {
    kindSelector(s);
    if (g.kind == Kind::None) return;

    int p = s.plant.outputs(), m = s.plant.inputs();
    ImGui::TextDisabled("click a cell to pair / unpair - the diagnostic is the editor");

    if (ImGui::BeginTable("##grid", m + 2, ImGuiTableFlags_Borders | ImGuiTableFlags_SizingFixedFit)) {
        ImGui::TableSetupColumn("scale", ImGuiTableColumnFlags_WidthFixed, 70);
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 26);
        for (int j = 0; j < m; ++j) {
            char h[8]; std::snprintf(h, sizeof h, "u%d", j);
            ImGui::TableSetupColumn(h, ImGuiTableColumnFlags_WidthFixed, 72);
        }
        ImGui::TableHeadersRow();

        for (int i = 0; i < p; ++i) {
            ImGui::PushID(i);
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::DragFloat("##sc", &g.scales[i], 0.01f, 0.001f, 1000.0f, "%.3g");
            ImGui::TableNextColumn();
            ImGui::Text("y%d", i);
            for (int j = 0; j < m; ++j) {
                ImGui::TableNextColumn();
                ImGui::PushID(j);
                int found = -1;
                for (size_t k = 0; k < g.loops.size(); ++k)
                    if (g.loops[k].out == i && g.loops[k].in == j) { found = (int)k; break; }

                char cell[32];
                if (d.has_rga) {
                    int r = (int)(std::find(d.sub_out.begin(), d.sub_out.end(), i) - d.sub_out.begin());
                    int c = (int)(std::find(d.sub_in.begin(), d.sub_in.end(), j) - d.sub_in.begin());
                    if (r < (int)d.sub_out.size() && c < (int)d.sub_in.size())
                        std::snprintf(cell, sizeof cell, "%.2f", d.lambda(r, c));
                    else std::snprintf(cell, sizeof cell, "-");
                } else if (d.has_share && j == d.share_in) {
                    std::snprintf(cell, sizeof cell, "%.0f%%", 100.0 * d.share[i]);
                } else {
                    std::snprintf(cell, sizeof cell, "-");
                }

                bool paired = found >= 0;
                if (paired) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.20f, 0.45f, 0.75f, 1.0f));
                else        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.16f, 0.16f, 0.18f, 1.0f));
                if (ImGui::Button(cell, ImVec2(-FLT_MIN, 24))) {
                    if (paired) { g.loops.erase(g.loops.begin() + found); g.selected = 0; }
                    else { Loop l; l.out = i; l.in = j; g.loops.push_back(l); g.selected = (int)g.loops.size() - 1; }
                }
                ImGui::PopStyleColor();
                if (paired && deadChannel(s.plant, i, j)) {
                    ImGui::SameLine(0, 2);
                    ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1), "!");
                }
                ImGui::PopID();
            }
            ImGui::PopID();
        }
        ImGui::EndTable();
    }

    omegaControl(-FLT_MIN);
    ImGui::Checkbox("sweep w across the Bode grid", &g.sweep);
    if (g.sweep) {
        float vals[40]; int n = 0;
        for (int k = 0; k < 40; ++k) {
            double f = 0.01 * std::pow(10.0, 4.0 * k / 39.0);
            double v = rgaNumberAt(s, f);
            vals[n++] = std::isnan(v) ? 0.0f : static_cast<float>(v);
        }
        ImGui::PlotLines("##sweep", vals, n, 0, "RGA number vs w  (0.01 - 100 Hz)", 0.0f, FLT_MAX,
                         ImVec2(-FLT_MIN, 50));
    }
    ImGui::TextColored(d.has_rga ? severityColor(d.rga_number) : ImVec4(0.7f, 0.7f, 0.7f, 1),
                       "%s", d.headline.c_str());
    deadWarnings(d);

    if (g.selected >= 0 && g.selected < (int)g.loops.size()) {
        auto& l = g.loops[g.selected];
        char hdr[64]; std::snprintf(hdr, sizeof hdr, "Gains for y%d <- u%d", l.out, l.in);
        ImGui::SeparatorText(hdr);
        gainsBlock(l);
    } else {
        ImGui::TextDisabled("select a paired cell to edit its gains");
    }
}

// ------------------------------------------------------------------- bench

void scaleBench() {
    if (std::getenv("PROTO_BENCH")) ImGui::SetNextItemOpen(true, ImGuiCond_Once);
    if (!ImGui::CollapsingHeader("Ki / Kd slider scale bench")) return;
    ImGui::TextDisabled("spec asks for linear 0-1 then log 1-50; must reach exactly 0.");
    const char* names[4] = {"0  linear 0-50", "1  ImGuiSliderFlags_Logarithmic",
                            "2  piecewise lin 0-1 / log 1-50", "3  DragFloat, speed ~ value"};
    for (int i = 0; i < 4; ++i) {
        ImGui::PushID(i);
        ImGui::RadioButton("##use", &g.ki_scale, i);
        ImGui::SameLine();
        ImGui::TextUnformatted(names[i]);
        gainSlider("##bench", &g.bench[i], 50.0f, i);
        ImGui::TextDisabled("  value %.6g   %s", g.bench[i], g.bench[i] == 0.0f ? "(exactly zero)" : "");
        ImGui::PopID();
    }
    ImGui::TextDisabled("the radio picks which scale the variants above use for Ki / Kd");
}

}  // namespace

// ---------------------------------------------------------------------- entry

void drawControllerSection(AppState& state) {
    static bool env_read = false;
    if (!env_read) {
        env_read = true;
        if (const char* v = std::getenv("PROTO_VARIANT")) g.variant = std::atoi(v) % 3;
    }
    reseedIfNeeded(state);
    for (auto& l : g.loops) {
        l.out = std::min(l.out, std::max(0, state.plant.outputs() - 1));
        l.in  = std::min(l.in,  std::max(0, state.plant.inputs() - 1));
    }

    // PROTOTYPE: park the panel at the bottom so the controller section is on
    // screen for a screenshot, with the tail of the plant section above it.
    static int frames = 0;
    if (frames < 200) { ++frames; ImGui::SetScrollY(ImGui::GetScrollMaxY()); }

    char hdr[64];
    std::snprintf(hdr, sizeof hdr, "Controller  [PROTOTYPE %c - %s]",
                  'A' + g.variant, kVariantNames[g.variant]);
    ImGui::SeparatorText(hdr);

    Diag d = computeDiag(state);
    switch (g.variant) {
        case 0: variantA(state, d); break;
        case 1: variantB(state, d); break;
        default: variantC(state, d); break;
    }
    ImGui::Spacing();
    scaleBench();
}

void drawSwitcher() {
    ImGuiIO& io = ImGui::GetIO();
    if (!io.WantTextInput) {
        if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow))  g.variant = (g.variant + 2) % 3;
        if (ImGui::IsKeyPressed(ImGuiKey_RightArrow)) g.variant = (g.variant + 1) % 3;
    }
    ImGuiViewport* vp = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(ImVec2(vp->WorkPos.x + vp->WorkSize.x * 0.5f,
                                   vp->WorkPos.y + vp->WorkSize.y - 16),
                            ImGuiCond_Always, ImVec2(0.5f, 1.0f));
    ImGui::SetNextWindowBgAlpha(0.95f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.95f, 0.85f, 0.15f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.05f, 0.05f, 0.05f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.85f, 0.72f, 0.05f, 1.0f));
    ImGui::Begin("##proto_switcher", nullptr,
                 ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
                 ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoNav |
                 ImGuiWindowFlags_NoFocusOnAppearing);
    if (ImGui::Button("<")) g.variant = (g.variant + 2) % 3;
    ImGui::SameLine();
    ImGui::Text("PROTOTYPE  %c - %s", 'A' + g.variant, kVariantNames[g.variant]);
    ImGui::SameLine();
    if (ImGui::Button(">")) g.variant = (g.variant + 1) % 3;
    ImGui::SameLine();
    ImGui::TextUnformatted("  (arrow keys)");
    ImGui::End();
    ImGui::PopStyleColor(3);
}

}  // namespace caliburn::proto
