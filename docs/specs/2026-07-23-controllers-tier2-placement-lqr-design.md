# Tier 2: Pole Placement + LQR Controllers — Design Spec

## Summary

Add pole placement (Ackermann's formula) and LQR (Linear Quadratic Regulator) controller types. Both are state-feedback controllers that compute a gain matrix K such that the closed-loop system `A - BK` has desired properties. They use the existing `stateFeedbackClose()` pipeline.

## Prerequisites

- Tier 1 (PID + Lead/Lag) should be complete first. Tier 1 expands the ControllerType enum and restructures the controller UI section. Tier 2 adds two more entries to that same enum.

## Design Decisions (confirmed by user)

- **Full-state feedback assumed**: No observer. The user provides desired poles (placement) or Q/R weights (LQR), and the system computes K. This is educational — real systems need observers, but that's a later tier.
- **Controllability check**: Both methods require a controllable (A, B). The existing `checkControllability()` in `system_properties.h` is used as a pre-check.

## Current Architecture (context for implementer)

### State feedback connection (`src/analysis/system_connect.cpp:58-66`)

```cpp
LinearSystem stateFeedbackClose(const LinearSystem& plant, const Eigen::MatrixXd& K) {
    LinearSystem result;
    result.A = plant.A - plant.B * K;
    result.B = plant.B;
    result.C = plant.C - plant.D * K;
    result.D = plant.D;
    return result;
}
```

This produces the closed-loop system with `u = -Kx + r`. No controller system (systems[1]) or open-loop system (systems[2]) — only plant and closed-loop are valid.

### Controllability check (`src/analysis/system_properties.h`)

```cpp
struct PropertyResult {
    Eigen::MatrixXd matrix;
    int rank;
    int required_rank;
    bool pass;  // rank == required_rank
};

PropertyResult checkControllability(const LinearSystem& sys);
```

### Recompute loop for GainMatrix (`src/visualizer.cpp:84-89`)

```cpp
else if (state.ctrl_type == caliburn::ControllerType::GainMatrix &&
         state.ctrl_K.size() > 0) {
    state.system_valid[1] = false;
    state.system_valid[2] = false;
    state.systems[3] = caliburn::stateFeedbackClose(state.plant, state.ctrl_K);
    state.system_valid[3] = true;
}
```

### Key files

| File | Role |
|---|---|
| `src/app_state.h` | ControllerType enum, AppState (ctrl_K, system_valid, etc.) |
| `src/analysis/system_connect.h/.cpp` | `stateFeedbackClose()` |
| `src/analysis/system_properties.h/.cpp` | `checkControllability()` |
| `src/panels/model_panel.cpp` | Controller UI section |
| `src/visualizer.cpp:60-94` | Recompute loop |
| `src/linear_system.h` | LinearSystem struct |
| `CMakeLists.txt` | C++17, Eigen 3.4 via FetchContent |

## Pole Placement (Ackermann's Formula)

### Algorithm (SISO only)

Given desired closed-loop poles `p_1, p_2, ..., p_n`, find K such that `eig(A - BK) = {p_1, ..., p_n}`.

1. Check controllability: `rank([B, AB, ..., A^(n-1)B]) = n`
2. Form desired characteristic polynomial: `φ(s) = (s - p_1)(s - p_2)...(s - p_n) = s^n + a_{n-1}s^{n-1} + ... + a_0`
3. Expand polynomial coefficients from roots
4. Evaluate matrix polynomial: `φ(A) = A^n + a_{n-1}A^{n-1} + ... + a_0·I`
5. Form controllability matrix: `Mc = [B, AB, A²B, ..., A^(n-1)B]`
6. Ackermann's formula: `K = [0 0 ... 0 1] · Mc^{-1} · φ(A)`

### Implementation Notes

- Complex conjugate poles must come in pairs (for real K)
- All current presets have n ≤ 10 — Ackermann is numerically safe
- Polynomial expansion from roots: iterative convolution `(s - p_i)` accumulated
- Matrix polynomial evaluation: Horner's method `φ(A) = (...((I·A + a_{n-1}·I)·A + a_{n-2}·I)·A + ... + a_0·I)`
- Uses only Eigen operations: matrix multiply, QR solve (for Mc inverse)

### Limitations

- SISO only (single-input). For MIMO plants, either:
  - Gray out pole placement and show "SISO only"
  - Or use eigenvalue assignment (Tier 3) which handles MIMO

### User Parameters

The user specifies n desired pole locations. For educational presets:

| Preset | n | Default poles | Rationale |
|---|---|---|---|
| First-Order | 1 | [-2] | Faster than default (-1) |
| Second-Order | 2 | [-1+j, -1-j] | ζ=0.707, ωn=√2 |
| Ball-Balancer | 4 | [-1, -2, -3, -4] | All real, increasing separation |
| Inverted Pendulum | 4 | [-1±j, -3±j] | Two pairs, well-damped |

### UI for Pole Placement

```
Pole Placement
├── ⚠ "SISO only" (if plant.inputs() > 1)
├── ⚠ "Not controllable" (if checkControllability fails)
├── Desired poles:
│   ├── p1: Re [===] Im [===]    (Re=-1.0, Im=0.0)
│   ├── p2: Re [===] Im [===]    (Re=-1.0, Im=1.0)
│   ├── p3: Re [===] Im [===]    (Re=-1.0, Im=-1.0)  ← auto-conjugate
│   └── ...
├── [Compute K]
├── K = [0.50  1.20  0.80  0.30]
└── CL poles: -1.0, -1.0±j
```

Conjugate pairing: when the user edits a pole with nonzero Im, the conjugate is automatically updated. This prevents invalid (complex K) configurations.

### Slider ranges for desired poles

| Param | Range | Default | Scale |
|---|---|---|---|
| Real part | -20 to 0 | -1.0 | Linear |
| Imaginary part | 0 to 20 | 0.0 | Linear |

Only stable poles allowed (Re < 0). Imaginary part ≥ 0, conjugate auto-generated.

## LQR Controller

### Algorithm

Minimize `J = ∫₀^∞ (x'Qx + u'Ru) dt`. Solution: `u = -Kx` where `K = R⁻¹B'P` and P solves the CARE:

```
A'P + PA - PBR⁻¹B'P + Q = 0
```

### CARE Solver: Matrix Sign Function (Drake's method)

**No external libraries needed.** Pure Eigen, ~40 lines of core logic.

```cpp
bool solveCARE(const Eigen::MatrixXd& A,
               const Eigen::MatrixXd& B,
               const Eigen::MatrixXd& Q,
               const Eigen::MatrixXd& R,
               Eigen::MatrixXd& P)
{
    const int n = A.rows();

    // 1. Validate R positive definite via Cholesky
    Eigen::LLT<Eigen::MatrixXd> R_chol(R);
    if (R_chol.info() != Eigen::Success) return false;

    // 2. Form 2n × 2n Hamiltonian
    Eigen::MatrixXd BRinvBt = B * R_chol.solve(B.transpose());
    Eigen::MatrixXd H(2*n, 2*n);
    H << A,    BRinvBt,
         Q,   -A.transpose();

    // 3. Matrix sign function iteration with determinant scaling
    Eigen::MatrixXd Z = H;
    const double tol = 1e-9;
    const int max_iter = 100;
    const double p = static_cast<double>(2 * n);

    for (int iter = 0; iter < max_iter; ++iter) {
        Eigen::MatrixXd Z_old = Z;
        double ck = std::pow(std::abs(Z.determinant()), -1.0 / p);
        Z *= ck;
        Z = Z - 0.5 * (Z - Z.inverse());
        if ((Z - Z_old).norm() < tol) break;
    }

    // 4. Extract P from sign matrix blocks
    Eigen::MatrixXd W11 = Z.block(0, 0, n, n);
    Eigen::MatrixXd W12 = Z.block(0, n, n, n);
    Eigen::MatrixXd W21 = Z.block(n, 0, n, n);
    Eigen::MatrixXd W22 = Z.block(n, n, n, n);

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd lhs(2*n, n), rhs(2*n, n);
    lhs << W12, W22 + I;
    rhs << W11 + I, W21;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        lhs, Eigen::ComputeThinU | Eigen::ComputeThinV);
    P = svd.solve(rhs);
    P = 0.5 * (P + P.transpose());  // symmetrize
    return true;
}
```

Converges quadratically, typically 15-25 iterations. All Eigen operations used (inverse, determinant, JacobiSVD, LLT, block) are available in Eigen 3.4 and Emscripten-compatible.

### Result Structure

```cpp
struct LqrResult {
    Eigen::MatrixXd K;     // m × n gain matrix
    Eigen::MatrixXd P;     // n × n CARE solution
    std::vector<std::complex<double>> closed_loop_poles;
    bool success;
    std::string error;
};

LqrResult computeLQR(const LinearSystem& sys,
                      const Eigen::MatrixXd& Q,
                      const Eigen::MatrixXd& R);
```

### User Parameters

Q and R are typically diagonal. The UI exposes diagonal entries as sliders:

```
LQR
├── ⚠ "Not controllable" (if checkControllability fails)
├── State weights Q (diagonal):
│   ├── q1 (position): [====log slider====]  default: 1.0
│   ├── q2 (velocity): [====log slider====]  default: 0.1
│   └── ...
├── Input weights R (diagonal):
│   ├── r1 (force):    [====log slider====]  default: 1.0
│   └── ...
├── [Compute K]
├── K = [1.00  1.73]
├── CL poles: -0.87±0.87j
└── Cost J = 3.14
```

### Slider Ranges

| Param | Range | Default | Scale |
|---|---|---|---|
| Q diagonal entries | 0.01 — 1000 | 1.0 | Logarithmic |
| R diagonal entries | 0.01 — 1000 | 1.0 | Logarithmic |

### Default Q/R per Preset

| Preset | Q | R |
|---|---|---|
| First-Order (n=1, m=1) | [1] | [1] |
| Second-Order (n=2, m=1) | diag(1, 0.1) | [1] |
| Ball-Balancer (n=4, m=2) | diag(50, 1, 50, 1) | diag(1, 1) |
| Inverted Pendulum (n=4, m=1) | diag(10, 1, 100, 1) | [1] |
| Quarter-Car (n=4, m=1) | diag(1, 1, 1, 1) | [1] |
| Double Mass-Spring (n=4, m=2) | diag(1, 0.1, 1, 0.1) | diag(1, 1) |

### Validation

Before computing LQR:
1. Check controllability (existing function)
2. Check Q is symmetric positive semidefinite
3. Check R is symmetric positive definite
4. After solving, verify all closed-loop poles have Re < 0

### Verification Test

For the second-order mass-spring-damper (m=1, b=1.4, k=1):
- A = [0 1; -1 -1.4], B = [0; 1]
- Q = I, R = [1]
- Expected: K ≈ [0.4142135624, 0.5463882256], CL poles ≈ -0.9731941 ± 0.6834521j
- Cross-check: residual `||A'P + PA - PBR⁻¹B'P + Q|| < 1e-8`

> **Corrected while implementing [#12](https://github.com/caliburn-engineering/caliburn/issues/12).**
> This spec originally gave K ≈ [0.414, 0.318] with poles ≈ -0.86 ± 0.86j. Those
> two numbers are inconsistent with each other and with the plant: K = [0.414,
> 0.318] places the poles at -0.859 ± 0.822j, while poles at -0.86 ± 0.86j
> require K = [0.479, 0.320]. Neither is the CARE solution.
>
> Solving the 2x2 Riccati equation by hand gives `P12 = sqrt(2) - 1` and
> `P22 = -1.4 + sqrt(0.96 + 2*sqrt(2))`, so `K = [P12, P22]` at R = 1; the
> stable eigenvalues of the Hamiltonian, computed with no Riccati solver
> involved, are -0.9731941 ± 0.6834521j. The residual of the corrected P is
> 5e-16. The first gain component was right, the second was not.
>
> Only the residual bound survived unchanged, which is the point the ticket
> makes about it: it holds against the equation itself rather than against a
> transcribed number, so it was the one line here that could not go stale.

## Connection to Existing Code

Both Pole Placement and LQR produce a gain matrix K. They reuse the existing GainMatrix path:

```cpp
if (state.ctrl_type == ControllerType::PolePlacement ||
    state.ctrl_type == ControllerType::LQR) {
    // K is computed by the UI (on param change) and stored in state.ctrl_K
    // Then the existing GainMatrix path handles the rest:
    state.systems[3] = stateFeedbackClose(state.plant, state.ctrl_K);
    state.system_valid[3] = true;
}
```

No changes to system_connect.cpp needed.

## New ControllerType Entries

After Tier 1 expands the enum, Tier 2 adds:

```cpp
enum class ControllerType {
    None,
    PID,            // Tier 1
    LeadLag,        // Tier 1
    PolePlacement,  // Tier 2
    LQR,            // Tier 2
    StateSpace,     // Existing
    GainMatrix,     // Existing
};
```

## Files Changed

| File | Change |
|---|---|
| `src/app_state.h` | Add PolePlacement and LQR to ControllerType. Add desired_poles, lqr_Q_diag, lqr_R_diag fields. |
| `src/analysis/pole_placement.h` | New file. `Eigen::MatrixXd ackermannPlace(const LinearSystem& sys, const std::vector<std::complex<double>>& desired_poles)` |
| `src/analysis/pole_placement.cpp` | New file. Ackermann's formula implementation. |
| `src/analysis/lqr.h` | New file. LqrResult struct, `computeLQR()` declaration. |
| `src/analysis/lqr.cpp` | New file. CARE solver (matrix sign function) + LQR computation. |
| `src/panels/model_panel.cpp` | Add PolePlacement and LQR UI sections within controller dropdown. |
| `src/visualizer.cpp` | Add PolePlacement and LQR cases to recompute loop (same path as GainMatrix). |
| `CMakeLists.txt` | Add pole_placement.cpp and lqr.cpp to analysis_lib. |
| `tests/test_pole_placement.cpp` | New test file. |
| `tests/test_lqr.cpp` | New test file. Verify CARE solution, K computation, CL stability. |

## Verification

1. Build: `cmake --build build` — clean compile
2. Select Second-Order preset, Pole Placement — specify desired poles [-2+j, -2-j]
3. Click Compute K — K displayed, CL poles shown, pole-zero plot shows new locations
4. Change desired poles — K and all plots update
5. Try on MIMO preset — see "SISO only" warning
6. Select First-Order preset, LQR — slide Q/R diagonal — K updates, CL poles move
7. Select Inverted Pendulum, LQR — verify it stabilizes (all CL poles in LHP)
8. All existing controller types (None, PID, Lead/Lag, SS, GainMatrix) still work
