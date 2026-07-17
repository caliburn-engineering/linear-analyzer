// src/visualizer.cpp
#include "app_state.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <cstdio>

#include "panels/model_panel.h"
#include "panels/properties_panel.h"
#include "panels/pole_zero_panel.h"
#include "panels/bode_panel.h"
#include "panels/nyquist_panel.h"
#include "panels/time_response_panel.h"

int main() {
    // --- Init GLFW ---
    glfwSetErrorCallback([](int err, const char* desc) {
        std::fprintf(stderr, "GLFW error %d: %s\n", err, desc);
    });
    if (!glfwInit()) return 1;

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    GLFWwindow* window = glfwCreateWindow(
        1600, 1000,
        "Linear System Analyzer \xe2\x80\x94 Caliburn", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::fprintf(stderr, "Failed to initialize GLAD\n");
        return 1;
    }

    // --- Init ImGui + ImPlot ---
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImPlot::CreateContext();

    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.IniFilename = "imgui.ini";

    // --- Load Unicode font ---
    ImFontConfig font_cfg;
    font_cfg.OversampleH = 2;
    font_cfg.OversampleV = 2;

    static const ImWchar glyph_ranges[] = {
        0x0020, 0x00FF,  // Basic Latin + Latin Supplement
        0x0370, 0x03FF,  // Greek and Coptic
        0x2070, 0x209F,  // Superscripts and Subscripts
        0x2200, 0x22FF,  // Mathematical Operators
        0x2190, 0x21FF,  // Arrows
        0,
    };

    ImFont* font = io.Fonts->AddFontFromFileTTF(
        "vendor/fonts/NotoSans-Regular.ttf", 16.0f, &font_cfg, glyph_ranges);
    if (!font) {
        std::fprintf(stderr, "Warning: NotoSans font not found, using default\n");
    }

    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = 6.0f;
    style.FrameRounding = 4.0f;
    style.GrabRounding = 4.0f;

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    // --- Init app state ---
    caliburn::AppState state;
    auto presets = caliburn::getBuiltinModels();
    state.plant = presets[0].system;
    caliburn::matrixToTextBuf(state.plant.A, state.A_text, sizeof(state.A_text));
    caliburn::matrixToTextBuf(state.plant.B, state.B_text, sizeof(state.B_text));
    caliburn::matrixToTextBuf(state.plant.C, state.C_text, sizeof(state.C_text));
    caliburn::matrixToTextBuf(state.plant.D, state.D_text, sizeof(state.D_text));

    glEnable(GL_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // --- Main loop ---
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- Full-viewport dockspace ---
        ImGuiViewport* vp = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(vp->WorkPos);
        ImGui::SetNextWindowSize(vp->WorkSize);
        ImGui::SetNextWindowViewport(vp->ID);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
        ImGui::Begin("DockSpace", nullptr,
            ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
            ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBackground);
        ImGui::PopStyleVar(3);
        ImGui::DockSpace(ImGui::GetID("MainDockSpace"), ImVec2(0, 0),
                         ImGuiDockNodeFlags_PassthruCentralNode);
        ImGui::End();

        // --- Recompute analysis if needed ---
        if (state.needs_recompute) {
            state.needs_recompute = false;

            // Plant is always system 0
            state.systems[0] = state.plant;
            state.system_valid[0] = true;

            // Build derived systems
            if (state.ctrl_type == caliburn::ControllerType::StateSpace) {
                state.systems[1] = state.ctrl_ss;
                state.system_valid[1] = true;
                state.systems[2] = caliburn::seriesConnect(state.ctrl_ss, state.plant);
                state.system_valid[2] = true;
                state.systems[3] = caliburn::feedbackConnect(state.systems[2]);
                state.system_valid[3] = true;
            } else if (state.ctrl_type == caliburn::ControllerType::GainMatrix) {
                state.system_valid[1] = false;
                state.system_valid[2] = false;
                state.systems[3] = caliburn::stateFeedbackClose(state.plant, state.ctrl_K);
                state.system_valid[3] = true;
            } else {
                state.system_valid[1] = false;
                state.system_valid[2] = false;
                state.system_valid[3] = false;
            }

            // Frequency response
            for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
                if (!state.system_valid[s]) continue;
                state.bode[s] = caliburn::computeBode(
                    state.systems[s], state.output_i, state.input_j,
                    state.freq_min_hz, state.freq_max_hz, state.num_freq_points);
                state.pole_zero[s] = caliburn::computePoleZero(
                    state.systems[s], state.output_i, state.input_j);
            }

            // All-channel Bode (when enabled)
            if (state.show_all_channels) {
                for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
                    if (!state.system_valid[s]) continue;
                    int p = state.systems[s].outputs();
                    int m = state.systems[s].inputs();
                    state.bode_grid[s].resize(p * m);
                    for (int i = 0; i < p; ++i) {
                        for (int j = 0; j < m; ++j) {
                            state.bode_grid[s][i * m + j] = caliburn::computeBode(
                                state.systems[s], i, j,
                                state.freq_min_hz, state.freq_max_hz,
                                state.num_freq_points);
                        }
                    }
                }
            }

            // Time response for Plant and Closed-Loop
            for (int s : {0, 3}) {
                if (!state.system_valid[s]) continue;
                switch (state.input_type) {
                    case caliburn::InputType::Step:
                        state.time_resp[s] = caliburn::computeStepResponse(
                            state.systems[s], state.time_input_j,
                            state.amplitude, state.duration, state.dt_sim);
                        break;
                    case caliburn::InputType::Impulse:
                        state.time_resp[s] = caliburn::computeImpulseResponse(
                            state.systems[s], state.time_input_j,
                            state.amplitude, state.duration, state.dt_sim);
                        break;
                    case caliburn::InputType::Ramp:
                        state.time_resp[s] = caliburn::computeRampResponse(
                            state.systems[s], state.time_input_j,
                            state.slope, state.duration, state.dt_sim);
                        break;
                }
            }

            // System properties
            state.controllability = caliburn::checkControllability(state.plant);
            state.observability = caliburn::checkObservability(state.plant);
            if (state.system_valid[3]) {
                state.cl_controllability =
                    caliburn::checkControllability(state.systems[3]);
                state.cl_observability =
                    caliburn::checkObservability(state.systems[3]);
            }

            // Root locus
            if (state.pz_mode == caliburn::PZMode::UnityFB) {
                state.root_locus = caliburn::computeRootLocus(
                    state.plant, state.output_i, state.input_j,
                    state.rl_k_min, state.rl_k_max, state.rl_num_points);
            } else if (state.pz_mode == caliburn::PZMode::StateFB &&
                       state.ctrl_type == caliburn::ControllerType::GainMatrix) {
                state.root_locus = caliburn::computeStateFeedbackLocus(
                    state.plant, state.ctrl_K,
                    state.rl_alpha_min, state.rl_alpha_max,
                    state.rl_num_points);
            }
        }

        // --- Panel toggle bar ---
        ImGui::Begin("##toggles", nullptr,
                     ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoScrollbar |
                     ImGuiWindowFlags_AlwaysAutoResize);
        auto toggleBtn = [](const char* label, bool& flag) {
            if (flag) {
                ImGui::PushStyleColor(ImGuiCol_Button,
                                      ImVec4(0.2f, 0.5f, 0.8f, 1.0f));
            }
            if (ImGui::Button(label)) flag = !flag;
            if (flag) ImGui::PopStyleColor();
            ImGui::SameLine();
        };
        toggleBtn("Pole-Zero", state.show_pole_zero);
        toggleBtn("Bode", state.show_bode);
        toggleBtn("Nyquist", state.show_nyquist);
        toggleBtn("Time Resp", state.show_time_response);
        ImGui::End();

        // --- Draw panels ---
        caliburn::drawModelPanel(state, presets);
        caliburn::drawPropertiesPanel(state);
        caliburn::drawPoleZeroPanel(state);
        caliburn::drawBodePanel(state);
        caliburn::drawNyquistPanel(state);
        caliburn::drawTimeResponsePanel(state);

        // --- Render ---
        ImGui::Render();
        int fb_w, fb_h;
        glfwGetFramebufferSize(window, &fb_w, &fb_h);
        glViewport(0, 0, fb_w, fb_h);
        glClearColor(0.08f, 0.08f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // --- Cleanup ---
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImPlot::DestroyContext();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
