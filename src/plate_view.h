// src/plate_view.h
#pragma once

#include "auto_balance.h"
#include "ball_sim.h"
#include "comparison_panel.h"
#include "plot_panel.h"
#include "renderer.h"
#include "table_kinematics.h"

#include <Eigen/Core>
#include <array>
#include <memory>
#include <string>
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

    /// Hand over the design surface's current answer, once per frame.
    ///
    /// `offered` is the caller's claim — that LQR is the selected controller,
    /// that the solve succeeded, and that the plant it was solved against is
    /// the cascade.  Only the model panel can know any of that.  The plate
    /// adds the two checks it alone can make, the gain's shape and whether the
    /// mechanism in `d` is the one being simulated, and states its own reason
    /// when either fails.
    ///
    /// `d.servo_tau` is honoured whether or not the design is offered: the
    /// legs have first-order lag in manual driving too, and the model panel's
    /// tau slider is the one place that number lives.
    void setDesign(const AutoBalanceDesign& d, bool offered,
                   const std::string& reason);

    /// Camera orbiting, for the scroll callback.
    OrbitCamera& camera() { return camera_; }

private:
    void drawControls();
    void drawBalanceControls();
    void drawMechanism();
    void resetAll();
    void resetBall();

    /// Command all three legs to one angle.  A command, not a teleport: the
    /// Home/Low/High buttons ask, and the legs arrive one lag later like every
    /// other command.  `snapServos` is the reset path, which has no lag to
    /// respect because it is not driving anything.
    void commandAllServos(float degrees);
    void snapServos(float degrees);

    /// Where the legs actually are, in radians — the argument every kinematics
    /// call wants, assembled in one place instead of five.
    std::array<double, 3> legsRad() const;

    /// True when the design on hand can actually drive this plate.
    bool designUsable() const;

    /// True while the balance loop, and not the sliders, owns the leg command.
    /// The single authority: `step` acts on it and `drawControls` greys the
    /// manual controls on it, so the two cannot disagree about who is driving.
    bool loopDriving() const;

    TableKinematics tk_;
    RollingBallDynamics ball_dynamics_;

    std::unique_ptr<LineRenderer> renderer_;

    // --- Servos ---
    // Two arrays, since the loop was closed: `alpha_cmd_deg_` is what the
    // sliders, the animation or the controller ASK for, `alpha_deg_` is where
    // the legs actually are.  Before the split the slider *was* the leg angle,
    // and u = -Kx around that is an algebraic loop on the leg states — the
    // servo lag the plant model claims had to become real.
    //
    // While the loop is closed the controller writes the command array, so the
    // (disabled) sliders read out what it is doing.
    float alpha_cmd_deg_[3] = {45.0f, 45.0f, 45.0f};
    float alpha_deg_[3] = {45.0f, 45.0f, 45.0f};
    bool link_servos_ = false;

    // --- Balance loop ---
    AutoBalanceDesign design_{};
    bool design_offered_ = false;
    std::string design_reason_ = "select LQR as the controller type";
    bool balance_engaged_ = false;
    bool balance_saturated_ = false;
    float sp_x_mm_ = 0.0f;  // ball setpoint, plate frame
    float sp_y_mm_ = 0.0f;

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
