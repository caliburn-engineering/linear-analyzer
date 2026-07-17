// src/panels/properties_panel.cpp
#include "properties_panel.h"

namespace caliburn {

namespace {
void drawPropertySection(const char* label, const PropertyResult& prop) {
    ImVec4 col = prop.pass ? ImVec4(0.3f, 1, 0.3f, 1) : ImVec4(1, 0.3f, 0.3f, 1);
    ImGui::TextColored(col, "%s: %s (rank %d / %d)",
                       label, prop.pass ? "PASS" : "FAIL",
                       prop.rank, prop.required_rank);

    if (ImGui::TreeNode(label)) {
        for (int i = 0; i < prop.matrix.rows() && i < 20; ++i) {
            std::string row;
            for (int j = 0; j < prop.matrix.cols() && j < 20; ++j) {
                char buf[16];
                std::snprintf(buf, sizeof(buf), "%8.4f", prop.matrix(i, j));
                row += buf;
            }
            ImGui::TextUnformatted(row.c_str());
        }
        ImGui::TreePop();
    }
}
}  // anonymous namespace

void drawPropertiesPanel(const AppState& state) {
    ImGui::Begin("System Properties");

    ImGui::SeparatorText("Plant");
    drawPropertySection("Controllability", state.controllability);
    drawPropertySection("Observability", state.observability);

    if (state.system_valid[3]) {
        ImGui::SeparatorText("Closed-Loop");
        drawPropertySection("CL Controllability", state.cl_controllability);
        drawPropertySection("CL Observability", state.cl_observability);
    }

    ImGui::End();
}

}  // namespace caliburn
