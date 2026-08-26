# Physical Parameter Tuning & Derivation Display — Implementation Plan

> **For agentic workers:** Implement this plan task-by-task, in order. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> *Status: implemented — see commits `190f843` and `e1d01cb`.*

**Goal:** Add physical parameter sliders and ODE→Laplace→State-Space derivation steps to each preset model, so users can tune physical quantities and see the full mathematical derivation.

**Architecture:** Extend `ModelEntry` with physical params, derivation steps, and a builder function. Add a mutable `current_params` copy to `AppState`. Render two new collapsible sections (Physical Parameters, Derivation) in the model panel before the existing State-Space section.

**Tech Stack:** C++17, Eigen 3.4, ImGui (docking), ImPlot

---

### Task 1: Extend Data Model (model_library.h + app_state.h)

**Files:**
- Modify: `src/analysis/model_library.h`
- Modify: `src/app_state.h`

- [ ] **Step 1: Add PhysicalParam, DerivationStep structs and extend ModelEntry in model_library.h**

Replace the entire file with:

```cpp
// src/analysis/model_library.h
#pragma once

#include "../linear_system.h"
#include <Eigen/Core>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace caliburn {

struct PhysicalParam {
    std::string name;      // Display name ("Mass")
    std::string symbol;    // Short symbol ("m")
    std::string unit;      // Unit string ("kg")
    float value;           // Default/current value
    float min_val;         // Slider minimum
    float max_val;         // Slider maximum
    bool logarithmic = false;
};

struct DerivationStep {
    std::string title;     // "Equations of Motion (ODE)"
    std::string content;   // Multi-line explanation with equations
};

struct ModelEntry {
    std::string name;
    std::string description;
    LinearSystem system;

    // Physical modeling (optional — empty means no physical params)
    std::vector<PhysicalParam> params;
    std::vector<DerivationStep> derivation;
    std::function<LinearSystem(const std::vector<PhysicalParam>&)> builder;
};

// Built-in preset models.
std::vector<ModelEntry> getBuiltinModels();

// Parse a matrix from a MATLAB-style string: "0 1; -2 -3"
// Rows separated by ';', columns by whitespace.
// Returns nullopt on parse failure.
std::optional<Eigen::MatrixXd> parseMatrix(const std::string& text);

// Format a matrix as MATLAB-style string: "0 1; -2 -3"
std::string matrixToString(const Eigen::MatrixXd& mat);

}  // namespace caliburn
```

- [ ] **Step 2: Add current_params to AppState in app_state.h**

After line 65 (`std::vector<RootLocusPoint> root_locus;`), insert:

```cpp
    // --- Physical parameters (mutable copy of current preset's params) ---
    std::vector<PhysicalParam> current_params;
```

- [ ] **Step 3: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Clean compilation (no errors). Warnings about unused fields are OK at this stage.

- [ ] **Step 4: Commit**

```bash
git add src/analysis/model_library.h src/app_state.h
git commit -m "feat: add PhysicalParam/DerivationStep structs, extend ModelEntry and AppState"
```

---

### Task 2: Add Physical Params, Builders, and Derivation Text to All Presets

**Files:**
- Modify: `src/analysis/model_library.cpp`

- [ ] **Step 1: Rewrite getBuiltinModels() with params, builders, and derivation for all 6 models**

Replace everything from `std::vector<ModelEntry> getBuiltinModels() {` (line 59) through the closing `}` of the function (line 173) with the following. Keep `parseMatrix` and `matrixToString` unchanged above it.

```cpp
std::vector<ModelEntry> getBuiltinModels() {
    std::vector<ModelEntry> models;

    // =====================================================================
    // 1. First-Order (K=1, tau=1)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "First-Order";
        entry.description = "1/(s+1), \xcf\x84=1s";

        entry.params = {
            {"Gain", "K", "", 1.0f, 0.01f, 10.0f, false},
            {"Time constant", "\xcf\x84", "s", 1.0f, 0.01f, 100.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double K = p[0].value, tau = p[1].value;
            LinearSystem sys;
            sys.A = (Eigen::MatrixXd(1, 1) << -1.0 / tau).finished();
            sys.B = (Eigen::MatrixXd(1, 1) << 1.0 / tau).finished();
            sys.C = (Eigen::MatrixXd(1, 1) << K).finished();
            sys.D = Eigen::MatrixXd::Zero(1, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A first-order system has one energy storage element.\n"
             "Examples:\n"
             "  - RC low-pass filter (capacitor stores charge)\n"
             "  - Thermal system (mass stores heat)\n"
             "  - Hydraulic tank (volume stores fluid)\n\n"
             "A single time constant \xcf\x84 governs the speed of response.\n"
             "DC gain K sets the steady-state output per unit input."},

            {"Governing ODE",
             "The first-order ODE relating input u(t) to output y(t):\n\n"
             "  \xcf\x84 \xc2\xb7 dy/dt + y = K \xc2\xb7 u\n\n"
             "where:\n"
             "  y(t) = output (voltage, temperature, level, ...)\n"
             "  u(t) = input (source, heat flux, flow, ...)\n"
             "  \xcf\x84 = time constant [s]\n"
             "  K = DC gain [-]"},

            {"Laplace Transform",
             "Apply Laplace transform (zero initial conditions):\n\n"
             "  \xcf\x84\xc2\xb7s\xc2\xb7Y(s) + Y(s) = K\xc2\xb7U(s)\n"
             "  Y(s)\xc2\xb7(\xcf\x84s + 1) = K\xc2\xb7U(s)\n\n"
             "Transfer function:\n"
             "  G(s) = Y(s)/U(s) = K / (\xcf\x84s + 1)\n\n"
             "Properties:\n"
             "  Pole: s = -1/\xcf\x84 (stable for \xcf\x84 > 0)\n"
             "  DC gain: G(0) = K\n"
             "  Bandwidth: \xcf\x89_bw = 1/\xcf\x84 rad/s\n"
             "  Step response settles in ~4\xcf\x84 seconds"},

            {"State-Space Form",
             "Define state x such that y = K\xc2\xb7x:\n\n"
             "  dx/dt = -(1/\xcf\x84)\xc2\xb7x + (1/\xcf\x84)\xc2\xb7u\n"
             "  y = K\xc2\xb7x\n\n"
             "Matrix form:\n"
             "  A = [-1/\xcf\x84]    B = [1/\xcf\x84]\n"
             "  C = [K]       D = [0]\n\n"
             "The eigenvalue of A is -1/\xcf\x84 = the pole of G(s)."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 2. Second-Order Mass-Spring-Damper (m=1, b=1.4, k=1)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Second-Order";
        entry.description = "Mass-spring-damper (\xcf\x89n=1, \xce\xb6=0.7)";

        entry.params = {
            {"Mass", "m", "kg", 1.0f, 0.1f, 100.0f, true},
            {"Damping", "b", "N\xc2\xb7s/m", 1.4f, 0.0f, 50.0f, false},
            {"Stiffness", "k", "N/m", 1.0f, 0.01f, 100.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double m = p[0].value, b = p[1].value, k = p[2].value;
            LinearSystem sys;
            sys.A = (Eigen::MatrixXd(2, 2) << 0, 1, -k / m, -b / m).finished();
            sys.B = (Eigen::MatrixXd(2, 1) << 0, 1.0 / m).finished();
            sys.C = (Eigen::MatrixXd(1, 2) << 1, 0).finished();
            sys.D = Eigen::MatrixXd::Zero(1, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A mass m connected to a fixed wall by a spring (stiffness k)\n"
             "and a viscous damper (coefficient b). An external force F\n"
             "is applied to the mass.\n\n"
             "This models many real systems:\n"
             "  - Vehicle suspension element\n"
             "  - Mechanical actuator with compliance\n"
             "  - Loudspeaker cone dynamics\n"
             "  - Structural vibration mode"},

            {"Equations of Motion (ODE)",
             "Newton's second law (sum of forces = m\xc2\xb7a):\n\n"
             "  m\xc2\xb7d\xc2\xb2x/dt\xc2\xb2 = -k\xc2\xb7x - b\xc2\xb7dx/dt + F\n\n"
             "Rearranged:\n"
             "  m\xc2\xb7x'' + b\xc2\xb7x' + k\xc2\xb7x = F\n\n"
             "where:\n"
             "  x = displacement from equilibrium [m]\n"
             "  F = applied force [N]\n"
             "  m = mass [kg]\n"
             "  b = viscous damping coefficient [N\xc2\xb7s/m]\n"
             "  k = spring stiffness [N/m]"},

            {"Laplace Transform",
             "Applying Laplace transform:\n\n"
             "  (m\xc2\xb7s\xc2\xb2 + b\xc2\xb7s + k)\xc2\xb7X(s) = F(s)\n\n"
             "Transfer function:\n"
             "  G(s) = X(s)/F(s) = 1 / (m\xc2\xb7s\xc2\xb2 + b\xc2\xb7s + k)\n\n"
             "Standard form (divide by m):\n"
             "  G(s) = (1/m) / (s\xc2\xb2 + 2\xce\xb6\xcf\x89n\xc2\xb7s + \xcf\x89n\xc2\xb2)\n\n"
             "Natural frequency and damping ratio:\n"
             "  \xcf\x89n = sqrt(k/m)       [rad/s]\n"
             "  \xce\xb6 = b / (2\xc2\xb7sqrt(m\xc2\xb7k))  [-]\n\n"
             "Poles: s = -\xce\xb6\xcf\x89n \xc2\xb1 \xcf\x89n\xc2\xb7sqrt(\xce\xb6\xc2\xb2 - 1)\n"
             "  \xce\xb6 < 1: underdamped (oscillates)\n"
             "  \xce\xb6 = 1: critically damped\n"
             "  \xce\xb6 > 1: overdamped (no oscillation)"},

            {"State-Space Form",
             "Choose states: x1 = position, x2 = velocity\n\n"
             "  x1' = x2\n"
             "  x2' = -(k/m)\xc2\xb7x1 - (b/m)\xc2\xb7x2 + (1/m)\xc2\xb7F\n\n"
             "Matrix form:\n"
             "  A = [  0       1   ]    B = [  0  ]\n"
             "      [-k/m   -b/m   ]        [1/m  ]\n\n"
             "  C = [1  0]              D = [0]\n\n"
             "Eigenvalues of A are the poles of G(s)."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 3. Ball-Balancer (ball rolling on tilting plate)
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Ball-Balancer";
        entry.description = "Ball on plate, 4 states, 2 in, 2 out";

        entry.params = {
            {"Gravity", "g", "m/s\xc2\xb2", 9.81f, 1.0f, 20.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double g = p[0].value;
            constexpr double k = 5.0 / 7.0;  // solid sphere rolling factor
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 2) = 1.0;
            sys.A(1, 3) = 1.0;
            sys.B = Eigen::MatrixXd::Zero(4, 2);
            sys.B(2, 1) = k * g;
            sys.B(3, 0) = k * g;
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 1) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 2);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A ball rolling freely on a flat plate that can tilt\n"
             "in two axes. The plate tilt angles \xce\xb8x and \xce\xb8y are\n"
             "the control inputs. Ball position (x, y) is measured.\n\n"
             "The ball rolls without slipping, so rotational inertia\n"
             "adds to the effective translational mass."},

            {"Equations of Motion (ODE)",
             "For a sphere rolling without slip on an incline:\n\n"
             "  (m + I/r\xc2\xb2)\xc2\xb7a = m\xc2\xb7g\xc2\xb7sin(\xce\xb8)\n\n"
             "For a solid sphere, I = (2/5)\xc2\xb7m\xc2\xb7r\xc2\xb2:\n"
             "  (m + 2m/5)\xc2\xb7a = m\xc2\xb7g\xc2\xb7sin(\xce\xb8)\n"
             "  (7/5)\xc2\xb7a = g\xc2\xb7sin(\xce\xb8)\n\n"
             "Linearize for small angles (sin\xce\xb8 \xe2\x89\x88 \xce\xb8):\n"
             "  x'' = (5g/7)\xc2\xb7\xce\xb8y     (tilt about y moves ball in x)\n"
             "  y'' = (5g/7)\xc2\xb7\xce\xb8x     (tilt about x moves ball in y)"},

            {"Transfer Function",
             "Each axis is a double integrator with gain:\n\n"
             "  G(s) = X(s)/\xce\x98(s) = (5g/7) / s\xc2\xb2\n\n"
             "Two poles at s = 0 (marginally stable).\n"
             "No damping — the ball never stops on its own.\n"
             "Active feedback control is essential."},

            {"State-Space Form",
             "States: [x, y, x', y']\n"
             "Inputs: [\xce\xb8x, \xce\xb8y]\n"
             "Outputs: [x, y]\n\n"
             "  A = [0 0 1 0]    B = [  0      0   ]\n"
             "      [0 0 0 1]        [  0      0   ]\n"
             "      [0 0 0 0]        [  0    5g/7  ]\n"
             "      [0 0 0 0]        [5g/7    0    ]\n\n"
             "  C = [1 0 0 0]    D = [0 0]\n"
             "      [0 1 0 0]        [0 0]\n\n"
             "Double integrator: needs at least PD control."},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 4. Inverted Pendulum on Cart
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Inverted Pendulum";
        entry.description = "Cart-pendulum (unstable)";

        entry.params = {
            {"Cart mass", "M", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Pendulum mass", "m", "kg", 0.1f, 0.01f, 10.0f, true},
            {"Pendulum length", "l", "m", 0.5f, 0.05f, 5.0f, true},
            {"Gravity", "g", "m/s\xc2\xb2", 9.81f, 1.0f, 20.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double M = p[0].value, m = p[1].value;
            double l = p[2].value, g = p[3].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 2) = -m * g / M;
            sys.A(2, 3) = 1.0;
            sys.A(3, 2) = (M + m) * g / (M * l);
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / M, 0,
                     -1.0 / (M * l))
                        .finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "A rigid pendulum (mass m, length l) mounted on a cart\n"
             "(mass M). The cart moves on a frictionless track driven\n"
             "by horizontal force F. The pendulum is balanced in the\n"
             "inverted (upright) position.\n\n"
             "Classic benchmark for control design — inherently unstable."},

            {"Equations of Motion (ODE)",
             "Nonlinear Newton-Euler equations:\n\n"
             "  (M+m)\xc2\xb7x'' + m\xc2\xb7l\xc2\xb7\xce\xb8''\xc2\xb7cos\xce\xb8 - m\xc2\xb7l\xc2\xb7\xce\xb8'\xc2\xb2\xc2\xb7sin\xce\xb8 = F\n"
             "  m\xc2\xb7l\xc2\xb7x''\xc2\xb7cos\xce\xb8 + m\xc2\xb7l\xc2\xb2\xc2\xb7\xce\xb8'' - m\xc2\xb7g\xc2\xb7l\xc2\xb7sin\xce\xb8 = 0\n\n"
             "Linearize about \xce\xb8=0 (upright):\n"
             "  cos\xce\xb8 \xe2\x89\x88 1, sin\xce\xb8 \xe2\x89\x88 \xce\xb8, \xce\xb8'\xc2\xb2 \xe2\x89\x88 0\n\n"
             "  (M+m)\xc2\xb7x'' + m\xc2\xb7l\xc2\xb7\xce\xb8'' = F\n"
             "  x'' + l\xc2\xb7\xce\xb8'' = g\xc2\xb7\xce\xb8\n\n"
             "Solve for accelerations:\n"
             "  x'' = -(m\xc2\xb7g/M)\xc2\xb7\xce\xb8 + (1/M)\xc2\xb7F\n"
             "  \xce\xb8'' = ((M+m)\xc2\xb7g/(M\xc2\xb7l))\xc2\xb7\xce\xb8 - (1/(M\xc2\xb7l))\xc2\xb7F"},

            {"Stability Analysis",
             "The linearized system has a positive eigenvalue:\n"
             "  s = +sqrt((M+m)\xc2\xb7g/(M\xc2\xb7l))\n\n"
             "This pole in the right half-plane makes the system\n"
             "open-loop unstable — the pendulum falls without\n"
             "active control.\n\n"
             "Increasing m/M ratio or l makes the unstable pole\n"
             "slower (easier to control). Decreasing l makes it\n"
             "faster (harder to balance)."},

            {"State-Space Form",
             "States: [x, x', \xce\xb8, \xce\xb8']\n"
             "Input: F (force on cart)\n"
             "Outputs: [x, \xce\xb8]\n\n"
             "  A = [0      1         0            0   ]\n"
             "      [0      0      -m\xc2\xb7g/M         0   ]\n"
             "      [0      0         0            1   ]\n"
             "      [0      0   (M+m)g/(Ml)        0   ]\n\n"
             "  B = [0; 1/M; 0; -1/(Ml)]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 5. Quarter-Car Suspension
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Quarter-Car";
        entry.description = "Suspension (sprung + unsprung)";

        entry.params = {
            {"Sprung mass", "ms", "kg", 300.0f, 50.0f, 2000.0f, true},
            {"Unsprung mass", "mu", "kg", 50.0f, 10.0f, 200.0f, true},
            {"Spring stiffness", "ks", "N/m", 20000.0f, 1000.0f, 200000.0f, true},
            {"Tire stiffness", "kt", "N/m", 200000.0f, 50000.0f, 500000.0f, true},
            {"Damping", "bs", "N\xc2\xb7s/m", 1000.0f, 100.0f, 10000.0f, true},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double ms = p[0].value, mu = p[1].value;
            double ks = p[2].value, kt = p[3].value, bs = p[4].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 0) = -ks / ms;
            sys.A(1, 1) = -bs / ms;
            sys.A(1, 2) = ks / ms;
            sys.A(1, 3) = bs / ms;
            sys.A(2, 3) = 1.0;
            sys.A(3, 0) = ks / mu;
            sys.A(3, 1) = bs / mu;
            sys.A(3, 2) = -(ks + kt) / mu;
            sys.A(3, 3) = -bs / mu;
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 0, 0, kt / mu).finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "Two-mass model of a vehicle suspension:\n"
             "  - Sprung mass ms: car body, passengers, cargo\n"
             "  - Unsprung mass mu: wheel, tire, axle, brakes\n"
             "  - Suspension spring ks and damper bs connect the masses\n"
             "  - Tire stiffness kt connects wheel to road\n\n"
             "Road surface profile zr is the disturbance input."},

            {"Equations of Motion (ODE)",
             "Newton's second law for each mass:\n\n"
             "Sprung mass (car body):\n"
             "  ms\xc2\xb7zs'' = -ks\xc2\xb7(zs - zu) - bs\xc2\xb7(zs' - zu')\n\n"
             "Unsprung mass (wheel):\n"
             "  mu\xc2\xb7zu'' = ks\xc2\xb7(zs - zu) + bs\xc2\xb7(zs' - zu') - kt\xc2\xb7(zu - zr)\n\n"
             "where:\n"
             "  zs = sprung mass displacement [m]\n"
             "  zu = unsprung mass displacement [m]\n"
             "  zr = road profile input [m]"},

            {"Design Trade-offs",
             "Increasing suspension stiffness (ks):\n"
             "  + Reduces body roll, better handling\n"
             "  - Harsher ride, less comfort\n\n"
             "Increasing damping (bs):\n"
             "  + Controls resonance peaks\n"
             "  - Transmits more high-frequency vibration\n\n"
             "Key frequencies:\n"
             "  Body bounce: ~(1/2\xcf\x80)\xc2\xb7sqrt(ks/ms) ~ 1 Hz\n"
             "  Wheel hop:   ~(1/2\xcf\x80)\xc2\xb7sqrt(kt/mu) ~ 10 Hz"},

            {"State-Space Form",
             "States: [zs, zs', zu, zu']\n"
             "Input: zr (road profile)\n"
             "Outputs: [zs, zu]\n\n"
             "  A = [  0        1          0          0     ]\n"
             "      [-ks/ms  -bs/ms     ks/ms      bs/ms   ]\n"
             "      [  0        0          0          1     ]\n"
             "      [ ks/mu   bs/mu  -(ks+kt)/mu  -bs/mu   ]\n\n"
             "  B = [0; 0; 0; kt/mu]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    // =====================================================================
    // 6. Double Mass-Spring-Damper
    // =====================================================================
    {
        ModelEntry entry;
        entry.name = "Double Mass-Spring";
        entry.description = "Two coupled oscillators";

        entry.params = {
            {"Mass 1", "m1", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Mass 2", "m2", "kg", 1.0f, 0.1f, 50.0f, true},
            {"Spring 1", "k1", "N/m", 1.0f, 0.01f, 100.0f, true},
            {"Spring 2", "k2", "N/m", 1.0f, 0.01f, 100.0f, true},
            {"Damper 1", "b1", "N\xc2\xb7s/m", 0.1f, 0.0f, 10.0f, false},
            {"Damper 2", "b2", "N\xc2\xb7s/m", 0.1f, 0.0f, 10.0f, false},
        };

        entry.builder = [](const std::vector<PhysicalParam>& p) -> LinearSystem {
            double m1 = p[0].value, m2 = p[1].value;
            double k1 = p[2].value, k2 = p[3].value;
            double b1 = p[4].value, b2 = p[5].value;
            LinearSystem sys;
            sys.A = Eigen::MatrixXd::Zero(4, 4);
            sys.A(0, 1) = 1.0;
            sys.A(1, 0) = -(k1 + k2) / m1;
            sys.A(1, 1) = -(b1 + b2) / m1;
            sys.A(1, 2) = k2 / m1;
            sys.A(1, 3) = b2 / m1;
            sys.A(2, 3) = 1.0;
            sys.A(3, 0) = k2 / m2;
            sys.A(3, 1) = b2 / m2;
            sys.A(3, 2) = -k2 / m2;
            sys.A(3, 3) = -b2 / m2;
            sys.B = (Eigen::MatrixXd(4, 1) << 0, 1.0 / m1, 0, 0).finished();
            sys.C = Eigen::MatrixXd::Zero(2, 4);
            sys.C(0, 0) = 1.0;
            sys.C(1, 2) = 1.0;
            sys.D = Eigen::MatrixXd::Zero(2, 1);
            return sys;
        };
        entry.system = entry.builder(entry.params);

        entry.derivation = {
            {"Physical System",
             "Two masses connected by springs and dampers in series.\n"
             "Mass 1 is attached to a wall by spring k1 and damper b1,\n"
             "and coupled to mass 2 by spring k2 and damper b2.\n"
             "Force F is applied to mass 1.\n\n"
             "Models: coupled oscillators, multi-story buildings,\n"
             "drivetrain compliance, multi-DOF vibration."},

            {"Equations of Motion (ODE)",
             "Newton's second law for each mass:\n\n"
             "Mass 1:\n"
             "  m1\xc2\xb7x1'' = -k1\xc2\xb7x1 - b1\xc2\xb7x1'\n"
             "             + k2\xc2\xb7(x2 - x1) + b2\xc2\xb7(x2' - x1') + F\n\n"
             "Expand:\n"
             "  m1\xc2\xb7x1'' = -(k1+k2)\xc2\xb7x1 - (b1+b2)\xc2\xb7x1'\n"
             "             + k2\xc2\xb7x2 + b2\xc2\xb7x2' + F\n\n"
             "Mass 2:\n"
             "  m2\xc2\xb7x2'' = -k2\xc2\xb7(x2 - x1) - b2\xc2\xb7(x2' - x1')\n"
             "  m2\xc2\xb7x2'' = k2\xc2\xb7x1 + b2\xc2\xb7x1' - k2\xc2\xb7x2 - b2\xc2\xb7x2'"},

            {"Modal Analysis",
             "The system has two natural modes:\n"
             "  Mode 1 (in-phase): both masses move together\n"
             "  Mode 2 (out-of-phase): masses move opposite\n\n"
             "Natural frequencies depend on mass and stiffness ratios.\n"
             "With equal masses and springs, the modes are well-separated\n"
             "and visible as two resonance peaks in the Bode plot.\n\n"
             "Changing the mass ratio shifts peak separation.\n"
             "Increasing damping reduces peak height and broadens them."},

            {"State-Space Form",
             "States: [x1, x1', x2, x2']\n"
             "Input: F (force on mass 1)\n"
             "Outputs: [x1, x2]\n\n"
             "  A = [    0          1        0       0    ]\n"
             "      [-(k1+k2)/m1 -(b1+b2)/m1 k2/m1 b2/m1]\n"
             "      [    0          0        0       1    ]\n"
             "      [  k2/m2      b2/m2   -k2/m2 -b2/m2  ]\n\n"
             "  B = [0; 1/m1; 0; 0]\n"
             "  C = [1 0 0 0; 0 0 1 0]\n"
             "  D = [0; 0]"},
        };

        models.push_back(std::move(entry));
    }

    return models;
}
```

- [ ] **Step 2: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Clean compilation.

- [ ] **Step 3: Run existing tests to verify nothing is broken**

Run: `cd build && ctest --output-on-failure 2>&1 | tail -10`
Expected: All tests pass (the `ModelEntry` struct change is backwards-compatible since the new fields have default constructors).

- [ ] **Step 4: Commit**

```bash
git add src/analysis/model_library.cpp
git commit -m "feat: add physical params, builders, and derivation text to all 6 presets"
```

---

### Task 3: Add Physical Parameters and Derivation UI to Model Panel

**Files:**
- Modify: `src/panels/model_panel.cpp`

- [ ] **Step 1: Update preset selection to copy current_params**

In `drawModelPanel`, inside the `if (ImGui::Combo("Preset", ...))` block (around line 111), after `extractTFFromSS(state);` and before `state.needs_recompute = true;`, add:

```cpp
        state.current_params = model.params;
```

- [ ] **Step 2: Add Physical Parameters section after description text**

After the description text block (after line 125: `ImGui::TextDisabled(...)`), before the matrix text fields (line 128: `ImGui::InputText("A", ...)`), insert:

```cpp
    // --- Physical parameters ---
    if (!state.current_params.empty()) {
        const auto& preset = presets[state.preset_index];
        if (ImGui::CollapsingHeader("Physical Parameters",
                                     ImGuiTreeNodeFlags_DefaultOpen)) {
            bool param_changed = false;
            for (int i = 0; i < static_cast<int>(state.current_params.size()); ++i) {
                auto& pp = state.current_params[i];
                char label[64];
                if (!pp.unit.empty())
                    std::snprintf(label, sizeof(label), "%s [%s]",
                                  pp.name.c_str(), pp.unit.c_str());
                else
                    std::snprintf(label, sizeof(label), "%s",
                                  pp.name.c_str());

                ImGui::TextUnformatted(label);
                ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                char slider_id[32];
                std::snprintf(slider_id, sizeof(slider_id), "##pp_%d", i);
                ImGuiSliderFlags flags = pp.logarithmic
                                             ? ImGuiSliderFlags_Logarithmic
                                             : 0;
                if (ImGui::SliderFloat(slider_id, &pp.value,
                                       pp.min_val, pp.max_val, "%.4g",
                                       flags)) {
                    param_changed = true;
                }
            }
            if (param_changed && preset.builder) {
                state.plant = preset.builder(state.current_params);
                matrixToTextBuf(state.plant.A, state.A_text,
                                sizeof(state.A_text));
                matrixToTextBuf(state.plant.B, state.B_text,
                                sizeof(state.B_text));
                matrixToTextBuf(state.plant.C, state.C_text,
                                sizeof(state.C_text));
                matrixToTextBuf(state.plant.D, state.D_text,
                                sizeof(state.D_text));
                extractTFFromSS(state);
                state.needs_recompute = true;
            }
        }
    }
```

- [ ] **Step 3: Add Derivation section after Physical Parameters**

Immediately after the Physical Parameters block (before the matrix text fields), insert:

```cpp
    // --- Derivation ---
    if (!presets.empty() &&
        !presets[state.preset_index].derivation.empty()) {
        if (ImGui::CollapsingHeader("Derivation")) {
            ImGui::PushID(state.preset_index);
            const auto& steps = presets[state.preset_index].derivation;
            for (size_t i = 0; i < steps.size(); ++i) {
                if (ImGui::TreeNode(steps[i].title.c_str())) {
                    ImGui::TextWrapped("%s", steps[i].content.c_str());
                    ImGui::TreePop();
                }
            }
            ImGui::PopID();
        }
    }
```

- [ ] **Step 4: Update Reset All to restore current_params**

In the Reset All block (around line 348-357), after `extractTFFromSS(state);` and before `state.needs_recompute = true;`, add:

```cpp
        state.current_params = presets[0].params;
```

- [ ] **Step 5: Build to verify compilation**

Run: `cmake --build build 2>&1 | tail -5`
Expected: Clean compilation.

- [ ] **Step 6: Run existing tests**

Run: `cd build && ctest --output-on-failure 2>&1 | tail -10`
Expected: All tests pass.

- [ ] **Step 7: Commit**

```bash
git add src/panels/model_panel.cpp
git commit -m "feat: add physical parameter sliders and derivation display to model panel"
```

---

### Task 4: Visual Verification

- [ ] **Step 1: Launch the visualizer**

Run: `./build/visualizer`

- [ ] **Step 2: Verify First-Order preset**

1. Select "First-Order" preset
2. Confirm "Physical Parameters" section is open with K and tau sliders
3. Slide K — Bode gain shifts, SS matrix C updates
4. Slide tau — Bode bandwidth shifts, pole moves on pole-zero plot
5. Expand "Derivation" header — confirm 4 tree nodes (Physical System, Governing ODE, Laplace Transform, State-Space Form)
6. Expand each node — confirm readable text with equations

- [ ] **Step 3: Verify Second-Order (MSD) preset**

1. Select "Second-Order" preset
2. Physical Parameters shows m, b, k sliders
3. Slide m — watch omega_n change in TF section, poles move
4. Slide b — watch zeta change, overshoot in time response changes
5. Slide k — watch omega_n shift
6. Derivation shows 4 steps with ODE, Laplace, and SS forms

- [ ] **Step 4: Verify Inverted Pendulum preset**

1. Select "Inverted Pendulum"
2. Physical Parameters shows M, m, l, g sliders
3. Confirm pole-zero plot shows RHP pole (unstable)
4. Slide l — watch unstable pole move
5. Derivation includes "Stability Analysis" step

- [ ] **Step 5: Verify Quarter-Car preset**

1. Select "Quarter-Car"
2. Physical Parameters shows ms, mu, ks, kt, bs sliders (5 params)
3. Slide ks — watch resonance peaks shift in Bode
4. Derivation includes "Design Trade-offs" step

- [ ] **Step 6: Verify Reset All**

1. Modify some physical params on any preset
2. Click "Reset All"
3. Confirm physical params reset to First-Order defaults

- [ ] **Step 7: Verify last-action-wins behavior**

1. Select Second-Order, adjust m slider
2. Adjust TF omega_n slider — physical param m does NOT change (correct: one-way flow)
3. Adjust A matrix slider — physical params do NOT change (correct)
