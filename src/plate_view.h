// src/plate_view.h
#pragma once

#include "ball_sim.h"
#include "comparison_panel.h"
#include "plot_panel.h"
#include "renderer.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <memory>
#include <vector>

struct GLFWwindow;

namespace caliburn {

/// The ball-balancer half of the merged application: the 3-RRS table, the ball
/// rolling on it, the 3D scene, and the panels that drive and plot them.
///
/// It owns no GLFW window, no ImGui context and no main loop — the analyzer's
/// entry point owns all three, and this is what the merge amounted to.  Frame
/// order is fixed by the caller: attach() once before the ImGui backend is
/// initialised, initGL() once after the GL loader is up, then step() ->
/// drawPanels() -> drawScene() every frame.
class PlateView {
public:
    PlateView();
    ~PlateView();

    // Non-copyable, non-movable: the plot configs hold pointers to the series
    // members, so relocating the object would dangle them.
    PlateView(const PlateView&) = delete;
    PlateView& operator=(const PlateView&) = delete;

    /// Install the scroll handler.  Must run BEFORE ImGui_ImplGlfw_InitForOpenGL,
    /// so that ImGui's backend chains to it rather than replacing it — the
    /// other order silently kills scroll-wheel zoom in every ImPlot panel.
    void attach(GLFWwindow* window);

    /// Build GL resources.  Requires a current context with GLAD loaded.
    void initGL();

    /// Release GL resources.  Must run while the context is still current —
    /// leaving it to the destructor deletes buffers after glfwTerminate().
    void shutdownGL();

    /// Advance kinematics, ball and plot buffers by one frame.
    /// Call between glfwPollEvents() and ImGui::NewFrame().
    void step(GLFWwindow* window, float dt);

    /// Draw the "Plate Control", "Plate Plots" and "Model Comparison" panels.
    void drawPanels();

    /// Draw the 3D scene into the currently-set GL viewport.
    void drawScene(float aspect);

    /// Camera orbiting, for the scroll callback.
    OrbitCamera& camera() { return camera_; }

private:
    void drawControls();
    void drawMechanism();
    void resetAll();
    void resetBall();
    void setAllServos(float degrees);

    TableKinematics tk_;
    RollingBallDynamics ball_dynamics_;

    std::unique_ptr<LineRenderer> renderer_;

    // --- Servos ---
    float alpha_deg_[3] = {45.0f, 45.0f, 45.0f};
    bool link_servos_ = false;

    // --- Camera ---
    OrbitCamera camera_;
    bool dragging_ = false;
    double last_mx_ = 0.0, last_my_ = 0.0;

    // --- Display ---
    bool show_axes_ = true;
    bool show_grid_ = true;
    bool show_joints_ = true;

    // --- Animation ---
    bool animate_ = false;
    float anim_time_ = 0.0f;
    float anim_speed_ = 1.0f;
    float anim_amplitude_ = 5.0f;

    // --- Computed each frame ---
    TablePose pose_{};
    TablePose home_{};
    FKResult fk_result_{};
    double condition_num_ = 0.0;
    double manipulability_ = 0.0;
    float sim_time_ = 0.0f;

    // --- Ball ---
    Eigen::Vector4d ball_state_ = Eigen::Vector4d::Zero();
    bool ball_enabled_ = true;
    bool ball_on_plate_ = true;
    bool ball_auto_reset_ = true;
    float ball_nudge_ = 0.15f;  // [m/s]

    // --- Plots ---
    PlotState plot_state_;
    TimeSeries s_phi_, s_theta_;
    TimeSeries s_a0_, s_a1_, s_a2_;
    TimeSeries s_cond_, s_zc_;
    TimeSeries s_bx_, s_by_;
    std::vector<PlotConfig> plots_;

    ComparisonPanel comparison_panel_;
};

}  // namespace caliburn
