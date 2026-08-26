# Tier 1: PID + Lead/Lag Controllers — Design Spec

## Summary

Replace the current controller dropdown (None / State-Space / Gain Matrix K) with an expanded set of controller types. Tier 1 adds PID (with all sub-combinations) and Lead/Lag compensators. These are dynamic output-feedback controllers that connect via the existing `seriesConnect` + `feedbackConnect` pipeline.

## Design Decisions (confirmed by user)

- **One-way data flow**: Controller params → SS matrices → closed-loop. No back-propagation.
- **Full-state feedback deferred**: PID and Lead/Lag use output feedback (unity feedback loop). Observer-based state feedback is a future addition.
- **MIMO diagonal PID**: Supported. One PID per channel, block-diagonal SS. UI shows a coupling warning when plant is strongly coupled.
- **Derivative filter mandatory**: Pure Kd·s is improper. Always use filtered derivative: Kd·s / (τf·s + 1).

## Current Architecture (context for implementer)

### ControllerType enum (`src/app_state.h:36`)

```cpp
enum class ControllerType { None, StateSpace, GainMatrix };
```

### Controller connection logic (`src/visualizer.cpp:60-94`)

Two paths exist:
1. **StateSpace**: `seriesConnect(ctrl_ss, plant)` → `feedbackConnect(open_loop)` → closed-loop
2. **GainMatrix**: `stateFeedbackClose(plant, K)` → closed-loop (no controller or open-loop system)

### System array (`src/app_state.h:49`)

```
systems[0] = Plant
systems[1] = Controller (StateSpace only)
systems[2] = Open-Loop = series(Controller, Plant) — computed but hidden in UI
systems[3] = Closed-Loop
```

### Connection functions (`src/analysis/system_connect.h`)

```cpp
LinearSystem seriesConnect(const LinearSystem& sys1, const LinearSystem& sys2);
LinearSystem feedbackConnect(const LinearSystem& open_loop);
LinearSystem stateFeedbackClose(const LinearSystem& plant, const Eigen::MatrixXd& K);
```

### Key files

| File | Role |
|---|---|
| `src/app_state.h` | ControllerType enum, AppState fields |
| `src/panels/model_panel.cpp:282-358` | Controller UI section |
| `src/visualizer.cpp:60-94` | Recompute loop (system composition) |
| `src/analysis/system_connect.h/.cpp` | Series, feedback, state-feedback connections |
| `src/analysis/system_properties.h/.cpp` | Controllability/observability checks |
| `src/linear_system.h` | LinearSystem struct (A, B, C, D) |
| `CMakeLists.txt` | Build config (C++17, Eigen 3.4, ImGui, ImPlot) |

## New ControllerType Enum

```cpp
enum class ControllerType {
    None,
    PID,         // SISO: single PID. MIMO: diagonal PIDs (one per channel)
    LeadLag,     // SISO: single lead, lag, or lead-lag. MIMO: diagonal.
    StateSpace,  // Existing: arbitrary SS controller
    GainMatrix,  // Existing: static state feedback K
};
```

## PID Controller

### State-Space Realization (parallel form — recommended)

Decompose `C(s) = Kp + Ki/s + Kd·s/(τf·s + 1)` into parallel components:

| Type | n | A | B | C | D |
|---|---|---|---|---|---|
| P | 0 | [] | [] | [] | [Kp] |
| PI | 1 | [0] | [1] | [Ki] | [Kp] |
| PD | 1 | [-1/τf] | [1/τf] | [-Kd/τf] | [Kp + Kd/τf] |
| PID | 2 | diag(0, -1/τf) | [1; 1/τf] | [Ki, -Kd/τf] | [Kp + Kd/τf] |

Physical meaning: x1 = integrator accumulation, x2 = derivative filter memory.

### MIMO Diagonal PID (m channels)

Block-diagonal construction. For m-input, m-output plant with full PID on each channel:

```
A_blk = blkdiag(A_1, A_2, ..., A_m)     (2m × 2m)
B_blk = blkdiag(B_1, B_2, ..., B_m)     (2m × m)
C_blk = blkdiag(C_1, C_2, ..., C_m)     (m × 2m)
D_blk = blkdiag(D_1, D_2, ..., D_m)     (m × m)
```

Each channel can independently be P, PI, PD, or PID (different state counts).

### Connection Topology

PID is an output-feedback controller. It uses the same pipeline as the existing StateSpace type:

```
error = reference - output
PID output = C(s) · error
Plant output = G(s) · PID output

series(PID_ss, plant) → feedbackConnect() → closed-loop
```

### AppState Fields for PID

```cpp
// Per-channel PID parameters
struct PIDParams {
    float Kp = 1.0f;
    float Ki = 0.0f;
    float Kd = 0.0f;
    float tau_f = 0.01f;  // derivative filter time constant
};

// In AppState:
std::vector<PIDParams> pid_params;  // one per channel (size = plant.inputs())
```

### Parameter Ranges (sliders)

| Param | Range | Default | Scale |
|---|---|---|---|
| Kp | 0.01 — 100 | 1.0 | Logarithmic |
| Ki | 0 — 50 | 0 | Linear 0—1, then log 1—50 |
| Kd | 0 — 20 | 0 | Linear 0—1, then log 1—20 |
| τf | 0.001 — 1.0 | 0.01 | Logarithmic |

### PID Builder Function

```cpp
// Build a SISO PID controller in state-space form
LinearSystem buildPID(float Kp, float Ki, float Kd, float tau_f);

// Build a MIMO diagonal PID controller
LinearSystem buildDiagonalPID(const std::vector<PIDParams>& params);
```

The builder checks which terms are active (Ki > 0, Kd > 0) and constructs the minimal-state realization:
- P-only (Ki=0, Kd=0): 0 states, D=[Kp]
- PI (Kd=0): 1 state
- PD (Ki=0): 1 state
- PID: 2 states

### MIMO Coupling Warning

When the user selects PID on a MIMO plant, check if the plant is strongly coupled. A simple heuristic:

```cpp
// Check if off-diagonal DC gains are significant compared to diagonal
// DC gain of G(s): G(0) = -C * A^{-1} * B + D
Eigen::MatrixXd dc_gain = -plant.C * plant.A.inverse() * plant.B + plant.D;
// Compare off-diagonal to diagonal magnitudes
```

If off-diagonal DC gains exceed 30% of diagonal, show:
`ImGui::TextColored(ImVec4(1, 0.8f, 0.2f, 1), "Plant appears coupled — MIMO controller recommended");`

## Lead/Lag Compensator

### Transfer Function

```
C(s) = Kc · (Ts + 1) / (αTs + 1)

α < 1 → Lead (phase advance)
α > 1 → Lag (gain boost at low freq)
α = 1 → Pure gain Kc
```

### State-Space Realization (1 state)

```
z = 1/T,  p = 1/(αT)

A = [-p]                    = [-1/(αT)]
B = [1]
C = [Kc · (z - p)]          = [Kc · (α - 1) / (αT)]
D = [Kc]
```

### Lead-Lag (cascade, 2 states)

Series connection of lead section and lag section:

```
z1 = 1/T_lead,  p1 = 1/(α_lead · T_lead)
z2 = 1/T_lag,   p2 = 1/(α_lag · T_lag)

A = [ -p1,       0   ]
    [ (z1-p1),  -p2  ]

B = [ 1 ]
    [ 1 ]

C = Kc · [ (z1-p1),  (z2-p2) ]

D = [ Kc ]
```

### AppState Fields for Lead/Lag

```cpp
enum class CompensatorMode { Lead, Lag, LeadLag };

struct LeadLagParams {
    CompensatorMode mode = CompensatorMode::Lead;
    float Kc = 1.0f;
    float alpha_lead = 0.1f;   // < 1
    float T_lead = 1.0f;
    float alpha_lag = 10.0f;   // > 1
    float T_lag = 1.0f;
};

// In AppState:
std::vector<LeadLagParams> leadlag_params;  // one per channel
```

### Parameter Ranges (sliders)

| Param | Range | Default | Scale |
|---|---|---|---|
| Kc | 0.1 — 100 | 1.0 | Logarithmic |
| α (lead) | 0.01 — 0.99 | 0.1 | Logarithmic |
| α (lag) | 1.01 — 100 | 10.0 | Logarithmic |
| T | 0.001 — 100 | 1.0 | Logarithmic |

### Builder Functions

```cpp
LinearSystem buildLeadLag(const LeadLagParams& params);
LinearSystem buildDiagonalLeadLag(const std::vector<LeadLagParams>& params);
```

### Connection Topology

Same as PID: `seriesConnect(compensator_ss, plant)` → `feedbackConnect()`.

## UI Layout (Model Panel Controller Section)

```
Controller
├── Type: [None ▼ | PID | Lead/Lag | State-Space | Gain Matrix K]
│
├── (PID selected, SISO plant)
│   ├── Kp [====slider====]
│   ├── Ki [====slider====]
│   ├── Kd [====slider====]
│   ├── τf [====slider====]
│   └── Formula: "C(s) = 1.0 + 0.5/s + 0.1s/(0.01s+1)"
│
├── (PID selected, MIMO plant)
│   ├── ⚠ "Plant appears coupled — MIMO controller recommended"  (if coupled)
│   ├── Channel 1:
│   │   ├── Kp, Ki, Kd, τf sliders
│   ├── Channel 2:
│   │   ├── Kp, Ki, Kd, τf sliders
│   └── ...
│
├── (Lead/Lag selected)
│   ├── Mode: [Lead | Lag | Lead-Lag]
│   ├── Kc [====slider====]
│   ├── (Lead) α [====slider====]  T [====slider====]
│   ├── (Lag)  α [====slider====]  T [====slider====]
│   └── Formula: "C(s) = 1.0·(s+1.0)/(s+10.0)"
│
├── (State-Space selected — unchanged from current)
└── (Gain Matrix K selected — unchanged from current)
```

## Recompute Loop Changes (`visualizer.cpp`)

The PID and Lead/Lag controller types follow the same path as the existing StateSpace type. The only change is how `ctrl_ss` is populated:

```cpp
if (state.ctrl_type == ControllerType::PID) {
    // Build SS from PID params
    if (is_siso) {
        state.ctrl_ss = buildPID(params[0].Kp, params[0].Ki, params[0].Kd, params[0].tau_f);
    } else {
        state.ctrl_ss = buildDiagonalPID(state.pid_params);
    }
    // Then fall through to existing StateSpace connection logic
}
```

The rest of the pipeline (seriesConnect → feedbackConnect → analysis) is unchanged.

## Files Changed

| File | Change |
|---|---|
| `src/app_state.h` | Extend ControllerType enum. Add PIDParams, LeadLagParams, CompensatorMode. Add pid_params, leadlag_params to AppState. |
| `src/analysis/controller_builders.h` | New file. buildPID, buildDiagonalPID, buildLeadLag, buildDiagonalLeadLag declarations. |
| `src/analysis/controller_builders.cpp` | New file. Implementations of all builder functions. |
| `src/panels/model_panel.cpp:282-358` | Rewrite controller section with new type dropdown, PID sliders, Lead/Lag sliders, coupling warning. |
| `src/visualizer.cpp:60-94` | Add PID and LeadLag cases to recompute loop (build ctrl_ss, then use existing StateSpace path). |
| `CMakeLists.txt` | Add controller_builders.cpp to analysis_lib sources. |
| `tests/test_controller_builders.cpp` | New test file for PID and Lead/Lag builder functions. |

## Verification

1. Build: `cmake --build build` — clean compile
2. Select First-Order preset, PID controller — see Kp/Ki/Kd/τf sliders
3. Slide Kp — Bode magnitude shifts, closed-loop poles move
4. Add Ki — steady-state error disappears in step response
5. Add Kd — overshoot reduces, phase margin improves
6. Select Lead/Lag, Lead mode — slide α, T — phase lead visible on Bode
7. Switch to Lag mode — low-frequency gain boost visible
8. Select a MIMO preset (Quarter-Car), PID — see per-channel sliders + coupling warning
9. Existing State-Space and Gain Matrix controllers still work unchanged
