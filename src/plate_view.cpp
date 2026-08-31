// src/plate_view.cpp
#include "plate_view.h"

// No GL header here: this file makes no GL calls of its own, and renderer.h
// (via plate_view.h) is the one place that decides between GLAD and Emscripten's
// ES 3.0 headers.
#include <GLFW/glfw3.h>
#include "imgui.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace caliburn {
namespace {

constexpr double kDeg = M_PI / 180.0;
constexpr int kBufSize = 4000;

constexpr double kBallRadius = kPlateBall.radius;  // [m], used when drawing

namespace col {
    constexpr std::array<float,4> ground    = {0.5f, 0.5f, 0.5f, 1.0f};
    constexpr std::array<float,4> ground_f  = {0.3f, 0.3f, 0.3f, 0.12f};
    constexpr std::array<float,4> table_c   = {0.3f, 0.7f, 1.0f, 1.0f};
    constexpr std::array<float,4> table_f   = {0.2f, 0.5f, 0.8f, 0.12f};
    constexpr std::array<float,4> leg_L1    = {1.0f, 0.8f, 0.2f, 1.0f};
    constexpr std::array<float,4> leg_L2    = {0.2f, 1.0f, 0.4f, 1.0f};
    constexpr std::array<float,4> joint_sph = {1.0f, 0.3f, 0.3f, 1.0f};
    constexpr std::array<float,4> joint_rev = {0.3f, 1.0f, 0.3f, 1.0f};
    constexpr std::array<float,4> grid      = {0.25f, 0.25f, 0.25f, 0.5f};
    constexpr std::array<float,4> ball      = {0.98f, 0.45f, 0.09f, 1.0f};
    constexpr std::array<float,4> ball_f    = {0.98f, 0.45f, 0.09f, 0.45f};
    constexpr std::array<float,4> ball_off  = {0.9f, 0.2f, 0.2f, 1.0f};
}

TableParams defaultTableParams() {
    TableParams p;
    p.R_ground  = 0.300;
    p.R_table   = 0.300;
    p.L1        = 0.150;
    p.L2        = 0.150;
    p.alpha_min = 10.0 * kDeg;
    p.alpha_max = 80.0 * kDeg;
    return p;
}

OrbitCamera defaultCamera() {
    OrbitCamera cam;
    cam.azimuth = 45.0f;
    cam.elevation = 30.0f;
    cam.distance = 1.0f;
    cam.target = Eigen::Vector3f(0, 0, 0.12f);
    return cam;
}

// The plate's own gravity.  Fixed here and compared against the model's `g`
// rather than followed: the dynamics object is built once, and a design solved
// on the moon must be refused, not quietly run on Earth.
constexpr double kGravity = 9.81;

RollingBallDynamics ballDynamicsFor(const TableParams& table) {
    PlateParams plate{table.R_table, kGravity};
    return RollingBallDynamics(kPlateBall, plate);
}

// Is anybody touching anything, on this frame?
//
// A poll, deliberately, and not a "has the visitor arrived" fact: the latch is
// one line at the call site in `stepAttract`, where it can be read next to the
// flag it sets.  Attract mode's disturbances are for an empty room — a ball
// being kicked by nobody while the visitor is driving the plate themselves is
// a puzzle, not a demonstration — so any click, any scroll, any key ends them.
//
// Read off ImGui rather than GLFW because ImGui is the one place that already
// merges mouse, wheel, keyboard and touch into a single input state — and
// touch is how half of the visitors this ticket exists for will arrive.  The
// named-key range is the sweep ImGui itself documents for input mappers, and
// it covers the mouse buttons too; the wheel is checked separately because a
// wheel is an axis rather than a button and reads as motion, not a press.
//
// This is the one part of attract mode with no test behind it: it needs a live
// ImGui context and a real input event, and wrapping an abstraction around one
// bool to get a seam would cost more than it proves.  It is checked in a
// browser instead — click, then thirty seconds of idle with no further kick.
bool visitorIsTouchingSomething() {
    const ImGuiIO& io = ImGui::GetIO();
    if (io.MouseWheel != 0.0f || io.MouseWheelH != 0.0f) return true;
    for (int k = ImGuiKey_NamedKey_BEGIN; k < ImGuiKey_NamedKey_END; ++k)
        if (ImGui::IsKeyDown(static_cast<ImGuiKey>(k))) return true;
    return false;
}

// The scroll wheel is the one input this class cannot read by polling: GLFW
// only delivers it as an event.  One instance drives the app, so a file-scope
// pointer is enough and keeps the callback free of captures.
PlateView* g_plate_view = nullptr;

void scrollCallback(GLFWwindow*, double, double yoffset) {
    if (!g_plate_view || ImGui::GetIO().WantCaptureMouse) return;
    OrbitCamera& cam = g_plate_view->camera();
    cam.distance *= (1.0f - 0.1f * static_cast<float>(yoffset));
    cam.distance = std::clamp(cam.distance, 0.2f, 3.0f);
}

}  // namespace

PlateView::PlateView()
    // tk_ is declared first, so its params are live here: the ball's plate and
      // the drawn plate are one object, not two that have to be kept in step.
    : tk_(defaultTableParams()),
      ball_dynamics_(ballDynamicsFor(tk_.params())),
      plot_state_(kBufSize),
      s_phi_("\xcf\x86", kBufSize),                 // phi
      s_theta_("\xce\xb8", kBufSize),               // theta
      s_a0_("\xce\xb1\xe2\x82\x80", kBufSize),      // alpha_0
      s_a1_("\xce\xb1\xe2\x82\x81", kBufSize),      // alpha_1
      s_a2_("\xce\xb1\xe2\x82\x82", kBufSize),      // alpha_2
      s_cond_("\xce\xba", kBufSize),                // kappa
      s_zc_("z", kBufSize),
      s_bx_("x", kBufSize),
      s_by_("y", kBufSize) {
    camera_ = defaultCamera();

    home_ = tk_.home_pose();
    pose_ = home_;

    // Attract mode's opening frame.  The ball is placed rather than kicked:
    // at t = 0 there is nothing on screen yet for a kick to be a change from,
    // and a ball already off centre is moving on the first frame that draws.
    ball_state_ = attractStart(attract_);

    plots_ = {
        {"Ball Position [mm]", {&s_bx_, &s_by_}, 0},
        {"Table Angles [deg]", {&s_phi_, &s_theta_}, 1},
        {"Servo Angles [deg]", {&s_a0_, &s_a1_, &s_a2_}, 2},
        {"Condition Number",   {&s_cond_}, 3},
        {"Table Height [mm]",  {&s_zc_}, 4},
    };
}

void PlateView::attach(GLFWwindow* window) {
    g_plate_view = this;
    glfwSetScrollCallback(window, scrollCallback);
}

PlateView::~PlateView() {
    // The scroll callback outlives this object otherwise, and GLFW would call
    // it with a dangling pointer.
    if (g_plate_view == this) g_plate_view = nullptr;
}

void PlateView::initGL() {
    renderer_ = std::make_unique<LineRenderer>();
}

void PlateView::shutdownGL() {
    renderer_.reset();
}

void PlateView::commandAllServos(float degrees) {
    for (int i = 0; i < 3; ++i) alpha_cmd_deg_[i] = degrees;
}

void PlateView::snapServos(float degrees) {
    commandAllServos(degrees);
    for (int i = 0; i < 3; ++i) alpha_deg_[i] = degrees;
}

std::array<double, 3> PlateView::legsRad() const {
    return {alpha_deg_[0] * kDeg, alpha_deg_[1] * kDeg, alpha_deg_[2] * kDeg};
}

bool PlateView::designUsable() const {
    return design_offered_ && gainFitsCascade(design_) &&
           samePlant(design_.mechanism, design_.gravity, tk_.params(), kGravity);
}

bool PlateView::loopDriving() const {
    // The ball is part of the precondition, not a detail: with the simulation
    // off there is nothing to balance, and the loop would hold the plate at
    // whatever tilt the frozen ball state asks for, forever.
    return balance_engaged_ && ball_enabled_ && designUsable();
}

void PlateView::setDesign(const AutoBalanceDesign& d, bool offered,
                          const std::string& reason) {
    design_ = d;
    design_offered_ = offered;
    design_reason_ = reason;

    if (offered && !gainFitsCascade(d)) {
        design_reason_ = "the gain is not 3 x 7 - this is not the cascade plant";
    } else if (offered && !designUsable()) {
        // The physical sliders move the plant the gain is designed against;
        // the simulated plate keeps the geometry and gravity it was built
        // with.  Refusing is the honest answer — engaging would drive one
        // plate with a gain solved for another, and nothing on screen would
        // say so.
        design_reason_ = "plant geometry or gravity differs from the plate";
    }

    // Losing the design mid-run drops the loop rather than freezing the last
    // command: a stale gain is not a controller.  The one place this happens.
    if (balance_engaged_ && !designUsable()) balance_engaged_ = false;

    // And the one place it is engaged without being asked.  Attract mode
    // cannot open with `balance_engaged_` simply set true: on the first frame
    // the LQR solve has not run yet, so the drop above would clear the flag
    // and nothing would ever set it again.  Engaging on the first usable
    // design instead means the demo starts balancing the moment it CAN, which
    // is a frame later and is what "already stabilising at load" amounts to.
    if (attract_running_ && !balance_engaged_ && ball_enabled_ && designUsable())
        balance_engaged_ = true;
}

void PlateView::resetBall() {
    ball_state_.setZero();
    ball_on_plate_ = true;
}

void PlateView::resetAll() {
    snapServos(45.0f);
    animate_ = false;
    balance_engaged_ = false;
    balance_saturated_ = false;
    balance_clipped_ = false;
    sp_x_mm_ = 0.0f;
    sp_y_mm_ = 0.0f;
    anim_time_ = 0.0f;
    home_ = tk_.home_pose();
    camera_ = defaultCamera();
    plot_state_.clear();
    plot_state_.paused = false;
    for (auto& pc : plots_)
        for (auto* s : pc.series)
            s->data.clear();
    plot_state_.markers.clear();
    sim_time_ = 0.0f;
    resetBall();
}

void PlateView::step(GLFWwindow* window, float dt) {
    ImGuiIO& io = ImGui::GetIO();

    // --- Mouse orbit (the central dock node is empty, so a drag there is ours) ---
    if (!io.WantCaptureMouse) {
        double mx, my;
        glfwGetCursorPos(window, &mx, &my);
        if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
            if (dragging_) {
                camera_.azimuth -= static_cast<float>(mx - last_mx_) * 0.3f;
                camera_.elevation += static_cast<float>(my - last_my_) * 0.3f;
                camera_.elevation = std::clamp(camera_.elevation, 2.0f, 89.0f);
            }
            dragging_ = true;
        } else {
            dragging_ = false;
        }
        last_mx_ = mx;
        last_my_ = my;
    }

    if (plot_state_.paused) return;

    // --- Where the leg commands come from ---
    //
    // Three writers, one array, in strict precedence: the closed loop, then
    // the animation, then whatever the sliders last left there.  The
    // controller reads the legs where they ARE and the ball where it IS, and
    // the servos move afterwards — the same causal order the closed-loop test
    // runs, which is what makes that test evidence about this code.
    const std::array<double, 3> alpha_rad = legsRad();

    if (loopDriving()) {
        const LegCommand c = legCommand(tk_, design_, alpha_rad, ball_state_,
                                        sp_x_mm_ * 1e-3, sp_y_mm_ * 1e-3);
        balance_saturated_ = c.saturated;
        balance_clipped_ = c.clipped_to_workspace;
        for (int i = 0; i < 3; ++i)
            alpha_cmd_deg_[i] = static_cast<float>(c.alpha_rad[i] / kDeg);
    } else if (animate_) {
        balance_saturated_ = false;
        balance_clipped_ = false;
        anim_time_ += dt * anim_speed_;
        const float amp = anim_amplitude_;
        const float w = 2.0f * static_cast<float>(M_PI) * 0.5f * anim_time_;
        alpha_cmd_deg_[0] = 45.0f + amp * std::sin(w);
        alpha_cmd_deg_[1] = 45.0f + amp * std::sin(w + 2.0f * static_cast<float>(M_PI) / 3.0f);
        alpha_cmd_deg_[2] = 45.0f + amp * std::sin(w + 4.0f * static_cast<float>(M_PI) / 3.0f);
    } else {
        balance_saturated_ = false;
        balance_clipped_ = false;
    }
    sim_time_ += dt;

    // --- Servos ---
    const std::array<double, 3> cmd_rad = {
        alpha_cmd_deg_[0] * kDeg, alpha_cmd_deg_[1] * kDeg, alpha_cmd_deg_[2] * kDeg};
    const std::array<double, 3> next =
        stepServosOnPlate(tk_, alpha_rad, cmd_rad, design_.servo_tau, dt);
    for (int i = 0; i < 3; ++i)
        alpha_deg_[i] = static_cast<float>(next[i] / kDeg);

    // --- Kinematics ---
    const std::array<double, 3> alpha_now = legsRad();

    // `solve_pose`, not `forward_kinematics`: the constraint equations have a
    // second root with the table folded flat on the base, and the retry this
    // used to do tested only convergence — which the folded root passes, with
    // a residual of 1e-11, in one iteration.  Once the pose fell onto it the
    // next frame was seeded from it and the plate stayed flat for the rest of
    // the session.  See issue #22.
    fk_result_ = tk_.solve_pose(alpha_now, home_);

    // The pose is adopted only on success.  Before, it was assigned whatever
    // the solver last held even when that solve had failed — so a leg triple
    // with no assembly at all still moved the plate, using a pose that solved
    // nothing.  Keeping the previous one is what a mechanism does when it is
    // driven into a singularity: it binds and stops.
    if (fk_result_.converged) {
        pose_ = fk_result_.pose;
        home_ = fk_result_.pose;
    }
    condition_num_ = tk_.condition_number(alpha_now, pose_);
    manipulability_ = tk_.manipulability(alpha_now, pose_);

    // --- Attract mode ---
    // After the controller has read the ball and before the ball is
    // integrated: a disturbance the loop can see on the frame it happens is
    // not a disturbance, it is a feedforward term.
    stepAttract();

    // --- Ball ---
    if (ball_enabled_ && ball_on_plate_) {
        ball_state_ = stepBall(ball_dynamics_, ball_state_, pose_, dt);
        if (!ballOnPlate(ball_state_, tk_.params().R_table, kBallRadius)) {
            ball_on_plate_ = false;
            if (ball_auto_reset_) resetBall();
        }
    }

    // --- Plots ---
    plot_state_.push_time(sim_time_);
    push_series(s_bx_,    static_cast<float>(ball_state_(0) * 1000), plot_state_);
    push_series(s_by_,    static_cast<float>(ball_state_(1) * 1000), plot_state_);
    push_series(s_phi_,   static_cast<float>(pose_.phi / kDeg), plot_state_);
    push_series(s_theta_, static_cast<float>(pose_.theta / kDeg), plot_state_);
    push_series(s_a0_,    alpha_deg_[0], plot_state_);
    push_series(s_a1_,    alpha_deg_[1], plot_state_);
    push_series(s_a2_,    alpha_deg_[2], plot_state_);
    push_series(s_cond_,  static_cast<float>(condition_num_), plot_state_);
    push_series(s_zc_,    static_cast<float>(pose_.z_c * 1000), plot_state_);
}

void PlateView::stepAttract() {
    if (!attract_running_) return;

    // The latch, and the only one: attract mode never resumes.  An idle
    // visitor is still a visitor, and a demo that started kicking the ball
    // again behind their back would be worse than one that stopped.
    if (visitorIsTouchingSomething()) {
        attract_running_ = false;
        return;
    }

    // One kick at a time, and it can only ever BE one: `sim_time_` advances by
    // a fixed 1/60 per frame, so a gap measured in seconds cannot be spanned
    // in a single step however long the frame took in real time.  The count is
    // here to make a kick undeliverable twice, not to catch up on a backlog
    // this clock cannot accumulate.
    if (attractKicksBy(attract_, sim_time_) > attract_kicks_done_) {
        ball_state_ += attractKick(attract_, attract_kicks_done_);
        ++attract_kicks_done_;
    }
}

void PlateView::drawPanels() {
    drawControls();
    draw_time_series_panel("Plate Plots", plot_state_, plots_);
    comparison_panel_.draw();
}

void PlateView::drawControls() {
    ImGui::Begin("Plate Control");

    // --- Ball ---
    ImGui::SeparatorText("Ball");
    ImGui::Checkbox("Simulate ball", &ball_enabled_);
    ImGui::SameLine();
    ImGui::Checkbox("Auto-reset", &ball_auto_reset_);

    if (ball_on_plate_) {
        ImGui::Text("Position: %+7.1f, %+7.1f mm",
                    ball_state_(0) * 1000, ball_state_(1) * 1000);
        ImGui::Text("Velocity: %+7.1f, %+7.1f mm/s",
                    ball_state_(2) * 1000, ball_state_(3) * 1000);
    } else {
        ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "Ball has left the plate");
    }

    if (ImGui::Button("Reset Ball")) resetBall();
    ImGui::SameLine();
    if (ImGui::Button("Nudge +x")) ball_state_(2) += ball_nudge_;
    ImGui::SameLine();
    if (ImGui::Button("Nudge +y")) ball_state_(3) += ball_nudge_;
    ImGui::SameLine();
    // The disturbance the loop is watched against: put the ball a fifth of the
    // plate out and let go.  A nudge tests recovery from a kick; this tests
    // recovery from a position, which is what the design surface is aimed at.
    if (ImGui::Button("Displace")) {
        ball_state_ << 0.06, -0.04, 0.0, 0.0;
        ball_on_plate_ = true;
    }
    ImGui::SliderFloat("Nudge [m/s]", &ball_nudge_, 0.02f, 0.5f, "%.2f");

    drawBalanceControls();

    // --- Servo Angles ---
    ImGui::SeparatorText("Servo Angles");

    const bool driven = loopDriving();
    if (driven) {
        // Not hidden — the sliders become the readout of what u = -Kx is
        // asking for, which is the most direct view of the loop working.
        ImGui::TextDisabled("commanded by the LQR gain");
    }
    ImGui::BeginDisabled(driven);

    ImGui::Checkbox("Link all servos", &link_servos_);

    const float a_min = static_cast<float>(tk_.params().alpha_min / kDeg);
    const float a_max = static_cast<float>(tk_.params().alpha_max / kDeg);

    if (link_servos_) {
        if (ImGui::SliderFloat("All##servo", &alpha_cmd_deg_[0], a_min, a_max, "%.1f deg")) {
            alpha_cmd_deg_[1] = alpha_cmd_deg_[0];
            alpha_cmd_deg_[2] = alpha_cmd_deg_[0];
        }
    } else {
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x80 [deg]", &alpha_cmd_deg_[0], a_min, a_max, "%.1f");
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x81 [deg]", &alpha_cmd_deg_[1], a_min, a_max, "%.1f");
        ImGui::SliderFloat("\xce\xb1\xe2\x82\x82 [deg]", &alpha_cmd_deg_[2], a_min, a_max, "%.1f");
    }

    if (ImGui::Button("Home (45)")) commandAllServos(45.0f);
    ImGui::SameLine();
    if (ImGui::Button("Low (20)")) commandAllServos(20.0f);
    ImGui::SameLine();
    if (ImGui::Button("High (70)")) commandAllServos(70.0f);
    ImGui::EndDisabled();

    // The sliders are the COMMAND; the legs lag behind it.  Without this line
    // the two readouts disagree on screen with nothing saying why.
    ImGui::Text("legs at %.1f, %.1f, %.1f deg (lag \xcf\x84 = %.3f s)",
                alpha_deg_[0], alpha_deg_[1], alpha_deg_[2], design_.servo_tau);

    ImGui::Spacing();
    ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.7f, 0.2f, 0.2f, 1.0f));
    if (ImGui::Button("Reset All", ImVec2(-1, 0))) resetAll();
    ImGui::PopStyleColor();

    // --- Animation ---
    ImGui::SeparatorText("Animation");
    ImGui::BeginDisabled(driven);
    ImGui::Checkbox("Animate", &animate_);
    if (animate_) {
        ImGui::SliderFloat("Speed", &anim_speed_, 0.1f, 5.0f, "%.1f");
        ImGui::SliderFloat("Amplitude [deg]", &anim_amplitude_, 1.0f, 20.0f, "%.1f");
    }
    ImGui::EndDisabled();

    // --- Pause ---
    ImGui::SeparatorText("Simulation");
    if (ImGui::Button(plot_state_.paused ? "Resume" : "Pause", ImVec2(-1, 0))) {
        plot_state_.paused = !plot_state_.paused;
        if (plot_state_.paused) {
            plot_state_.pause_t_max = plot_state_.latest_time();
            plot_state_.pause_t_min = plot_state_.pause_t_max - plot_state_.time_window;
        }
    }

    // --- Table Pose ---
    ImGui::SeparatorText("Table Pose (FK)");

    auto ok_col = [](bool ok) {
        return ok ? ImVec4(0.3f, 1.0f, 0.3f, 1.0f) : ImVec4(1.0f, 0.3f, 0.3f, 1.0f);
    };
    ImGui::TextColored(ok_col(fk_result_.converged),
                       "FK: %s (%d iter)", fk_result_.converged ? "OK" : "FAIL",
                       fk_result_.iterations);
    ImGui::Text("Roll:   %+.2f deg", pose_.phi / kDeg);
    ImGui::Text("Pitch:  %+.2f deg", pose_.theta / kDeg);
    ImGui::Text("Heave:  %.1f mm",   pose_.z_c * 1000);

    // --- Jacobian ---
    ImGui::SeparatorText("Jacobian Analysis");
    auto cond_col = [](double c) -> ImVec4 {
        if (c < 10) return {0.3f, 1.0f, 0.3f, 1.0f};
        if (c < 20) return {1.0f, 0.8f, 0.2f, 1.0f};
        return {1.0f, 0.3f, 0.3f, 1.0f};
    };
    ImGui::TextColored(cond_col(condition_num_),
                       "Condition #:   %.2f", condition_num_);
    ImGui::Text("Manipulability: %.4f", manipulability_);

    const float cond_frac = std::min(static_cast<float>(condition_num_ / 50.0), 1.0f);
    ImGui::ProgressBar(cond_frac, ImVec2(-1, 0),
                       condition_num_ < 10 ? "Good" :
                       condition_num_ < 20 ? "Degraded" : "Poor");

    // --- Display ---
    ImGui::SeparatorText("Display");
    ImGui::Checkbox("Grid", &show_grid_);
    ImGui::SameLine();
    ImGui::Checkbox("Axes", &show_axes_);
    ImGui::SameLine();
    ImGui::Checkbox("Joints", &show_joints_);

    // --- Camera ---
    ImGui::SeparatorText("Camera");
    ImGui::SliderFloat("Azimuth",   &camera_.azimuth,   -180.0f, 180.0f, "%.0f deg");
    ImGui::SliderFloat("Elevation", &camera_.elevation,  5.0f, 89.0f, "%.0f deg");
    ImGui::SliderFloat("Distance",  &camera_.distance,   0.3f, 2.0f, "%.2f m");
    if (ImGui::Button("Reset Camera")) camera_ = defaultCamera();

    // --- Velocity Jacobian (collapsed) ---
    if (ImGui::CollapsingHeader("Velocity Jacobian")) {
        const Eigen::Matrix3d Jv = tk_.velocity_jacobian(legsRad(), pose_);
        ImGui::Text("        servo0    servo1    servo2");
        ImGui::Text("roll   %+.4f   %+.4f   %+.4f", Jv(0,0), Jv(0,1), Jv(0,2));
        ImGui::Text("pitch  %+.4f   %+.4f   %+.4f", Jv(1,0), Jv(1,1), Jv(1,2));
        ImGui::Text("heave  %+.4f   %+.4f   %+.4f", Jv(2,0), Jv(2,1), Jv(2,2));
    }

    ImGui::End();
}

// The one control the whole ticket is about: hand the plate over to the gain
// the model panel just designed.
//
// The three servo sliders and u = -Kx want to write the same array, so they
// cannot both be live.  The loop wins while it is engaged and the sliders go
// read-only rather than disappearing — a disabled slider that keeps moving is
// the clearest possible statement of what the controller is doing, and there
// is nothing left for a manual command to mean.  The setpoint takes over the
// steering job, and it is a BALL position rather than a plate tilt, because
// that is the quantity the design is regulating.
void PlateView::drawBalanceControls() {
    ImGui::SeparatorText("Balance Loop");

    // Read-only here.  `setDesign` is the one place the loop is dropped, so a
    // draw pass cannot disagree with it about whether the loop is running.
    if (!designUsable()) {
        ImGui::BeginDisabled(true);
        bool off = false;
        ImGui::Checkbox("Engage", &off);
        ImGui::EndDisabled();
        ImGui::TextWrapped("unavailable: %s", design_reason_.c_str());
        return;
    }

    ImGui::BeginDisabled(!ball_enabled_);
    ImGui::Checkbox("Engage", &balance_engaged_);
    ImGui::EndDisabled();
    if (!ball_enabled_) {
        ImGui::TextWrapped("unavailable: there is no ball to balance");
        return;
    }
    if (!balance_engaged_) {
        ImGui::TextDisabled("gain ready - %d x %d, legs from %.1f deg",
                            (int)design_.K.rows(), (int)design_.K.cols(),
                            design_.home_leg_rad / kDeg);
        return;
    }

    // A ball at rest anywhere on a flat plate is an equilibrium, so this
    // setpoint needs no feedforward: it enters as a reference STATE and the
    // regulator is already a tracker.  Bounded well inside the plate, since
    // the edge is where the ball leaves.
    const float lim = static_cast<float>((tk_.params().R_table - 0.05) * 1000.0);
    ImGui::SliderFloat("set x [mm]", &sp_x_mm_, -lim, lim, "%.0f");
    ImGui::SliderFloat("set y [mm]", &sp_y_mm_, -lim, lim, "%.0f");
    if (ImGui::Button("Centre setpoint")) { sp_x_mm_ = 0.0f; sp_y_mm_ = 0.0f; }

    const double err = std::hypot(ball_state_(0) - sp_x_mm_ * 1e-3,
                                  ball_state_(1) - sp_y_mm_ * 1e-3);
    ImGui::Text("error: %6.1f mm", err * 1000.0);

    if (balance_clipped_) {
        // A different thing from saturation, and worth its own line: the legs
        // are inside their travel and the controller is still not getting what
        // it asked for, because the pose it wants is not one this mechanism
        // has.  See CONTEXT.md, "Workspace vs. servo box".
        ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                           "command clipped to the plate's workspace");
    }

    if (balance_saturated_) {
        // Not a failure — the legs really do stop at 10 and 80 degrees.  But a
        // tuning that lives against the stops is not the tuning that was
        // designed, and the closed-loop poles on screen stop describing it.
        ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.2f, 1.0f),
                           "leg command saturated");
    }
}

void PlateView::drawScene(float aspect) {
    if (!renderer_) return;

    const Eigen::Matrix4f proj = perspective(45.0f, aspect, 0.01f, 10.0f);
    const Eigen::Matrix4f view = camera_.view_matrix();
    renderer_->set_vp(proj * view);
    renderer_->begin();
    drawMechanism();
    renderer_->flush();
}

void PlateView::drawMechanism() {
    LineRenderer& lr = *renderer_;
    const TableParams& p = tk_.params();

    // Ground grid
    if (show_grid_) {
        const double s = 0.5;
        for (double x = -s; x <= s + 1e-9; x += 0.05) {
            lr.line({x, -s, 0}, {x, s, 0}, col::grid, 1.0f);
            lr.line({-s, x, 0}, {s, x, 0}, col::grid, 1.0f);
        }
    }

    if (show_axes_) lr.axes({0, 0, 0}, 0.08);

    // Ground: filled disc + edge
    lr.disc({0, 0, 0}, p.R_ground, {0, 0, 1}, col::ground_f, 64);
    lr.circle({0, 0, 0}, p.R_ground, {0, 0, 1}, col::ground, 64, 2.0f);

    // Table: filled disc + edge
    const Eigen::Matrix3d R = tk_.table_rotation(pose_.phi, pose_.theta);
    const Eigen::Vector3d tc(0, 0, pose_.z_c);
    const Eigen::Vector3d tn = R * Eigen::Vector3d(0, 0, 1);

    lr.disc(tc, p.R_table, tn, col::table_f, 64);
    lr.circle(tc, p.R_table, tn, col::table_c, 64, 2.5f);

    // Table cross-hairs
    const Eigen::Vector3d tx = R * Eigen::Vector3d(p.R_table * 0.8, 0, 0);
    const Eigen::Vector3d ty = R * Eigen::Vector3d(0, p.R_table * 0.8, 0);
    lr.line(tc - tx, tc + tx, {0.3f, 0.6f, 0.9f, 0.5f}, 1.0f);
    lr.line(tc - ty, tc + ty, {0.3f, 0.6f, 0.9f, 0.5f}, 1.0f);

    // Legs
    for (int i = 0; i < 3; ++i) {
        const double alpha_i = alpha_deg_[i] * kDeg;
        const Eigen::Vector3d G = tk_.ground_point(i);
        const Eigen::Vector3d K = tk_.knee_position(i, alpha_i);
        const Eigen::Vector3d P = tk_.table_point_world(i, pose_);

        lr.line(G, K, col::leg_L1, 3.0f);
        lr.line(K, P, col::leg_L2, 3.0f);

        if (show_joints_) {
            lr.point(G, col::joint_sph, 0.008f, 4.0f);
            lr.point(K, col::joint_rev, 0.008f, 4.0f);
            lr.point(P, col::joint_sph, 0.008f, 4.0f);
        }
        if (show_axes_) {
            const Eigen::Vector3d t = tk_.tangent_dir(i);
            lr.line(G - 0.02 * t, G + 0.02 * t, {0.7f, 0.7f, 0.0f, 0.6f}, 1.5f);
        }
    }

    if (show_joints_) lr.point(tc, col::table_c, 0.006f, 3.0f);

    // Ball.  Its state is in the plate's own frame, so the same rotation that
    // places the leg attachment points places the ball — lifted off the surface
    // by its radius along the plate normal.
    if (ball_enabled_) {
        const Eigen::Vector3d bc =
            tc + R * Eigen::Vector3d(ball_state_(0), ball_state_(1), kBallRadius);
        const auto& edge = ball_on_plate_ ? col::ball : col::ball_off;

        // A wire sphere: a filled disc face-on to the plate, plus three great
        // circles.  The line renderer draws no triangulated solids, and three
        // circles read as a sphere where one circle reads as a hole.
        lr.disc(bc, kBallRadius, tn, col::ball_f, 24);
        lr.circle(bc, kBallRadius, tn, edge, 24, 2.0f);
        lr.circle(bc, kBallRadius, R * Eigen::Vector3d(1, 0, 0), edge, 24, 1.5f);
        lr.circle(bc, kBallRadius, R * Eigen::Vector3d(0, 1, 0), edge, 24, 1.5f);

        // A dropline to the plate centre-plane, so the ball's offset from the
        // middle of the plate is readable from any camera angle.
        const Eigen::Vector3d foot =
            tc + R * Eigen::Vector3d(ball_state_(0), ball_state_(1), 0.0);
        lr.line(tc, foot, {0.98f, 0.45f, 0.09f, 0.35f}, 1.0f);
    }
}

}  // namespace caliburn
