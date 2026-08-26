# Tier 3: Eigenvalue & Eigenstructure Assignment — Design Spec

## Summary

Add MIMO-capable pole placement via eigenvalue assignment (Kautsky-Nichols-Van Dooren method) and eigenstructure assignment. These generalize the SISO Ackermann method from Tier 2 to multi-input systems, giving the user control over both pole locations and (optionally) the response mode shapes.

## Prerequisites

- Tier 1 (PID + Lead/Lag) and Tier 2 (Pole Placement + LQR) should be complete.
- The ControllerType enum and controller UI structure from Tiers 1-2 are assumed.

## Design Decisions (confirmed by user)

- **Full-state feedback assumed** (same as Tier 2). Observer-based output feedback is a future addition.
- **Eigenstructure (eigenvector specification) UI deferred**: The algorithm supports it, but the eigenvector input UI is complex. Phase 1 implements eigenvalue assignment only. Phase 2 adds eigenvector direction specification.

## Current Architecture (context for implementer)

### Relevant existing code

| File | Role |
|---|---|
| `src/analysis/system_connect.h/.cpp` | `stateFeedbackClose(plant, K)` — produces A-BK closed-loop |
| `src/analysis/system_properties.h/.cpp` | `checkControllability(sys)` — rank check |
| `src/analysis/pole_placement.h/.cpp` | Ackermann's formula (Tier 2, SISO only) |
| `src/analysis/lqr.h/.cpp` | LQR solver (Tier 2) |
| `src/app_state.h` | ControllerType enum, desired_poles, ctrl_K |
| `src/linear_system.h` | LinearSystem struct (A, B, C, D with n states, m inputs, p outputs) |
| `CMakeLists.txt` | C++17, Eigen 3.4 |

### Connection path (same as GainMatrix / Pole Placement / LQR)

```cpp
state.ctrl_K = computedK;
state.systems[3] = stateFeedbackClose(state.plant, state.ctrl_K);
```

## Eigenvalue Assignment (MIMO Pole Placement)

### Problem Statement

Given plant (A, B) with n states and m inputs (m > 1), find gain matrix K (m × n) such that `eig(A - BK) = {λ_1, λ_2, ..., λ_n}` (desired poles).

### Algorithm: Kautsky-Nichols-Van Dooren (KNV)

This is the standard method for MIMO eigenvalue assignment. It uses SVD to find the freedom in eigenvector selection when m > 1.

For each desired eigenvalue λ_i:

1. Form the matrix `[λ_i·I - A | B]` (n × (n+m))
2. Compute its null space via SVD — this gives a basis for vectors `[v; w]` such that `(λ_i·I - A)v = Bw`
3. The null space has dimension m (for a controllable system with m inputs)
4. Choose one vector from this null space: the eigenvector v_i is the top n components, and k_i = -w_i (the bottom m components give the feedback direction)
5. After processing all n eigenvalues, form: `K = W · V^{-1}` where V = [v_1 ... v_n] and W = [w_1 ... w_n]

### Complex Conjugate Handling

For real K, complex eigenvalues must come in conjugate pairs. When processing λ = σ + jω:
- Compute the null space of `[λI - A | B]` using complex arithmetic
- The eigenvector for λ* = σ - jω is the conjugate of the eigenvector for λ
- Use real and imaginary parts: `v_real = Re(v_complex)`, `v_imag = Im(v_complex)`
- These form two columns of V (and corresponding columns of W)

### Eigen API Usage

All operations are built into Eigen 3.4:

```cpp
// Complex null space via SVD
Eigen::MatrixXcd M(n, n + m);
M << lambda * Eigen::MatrixXcd::Identity(n, n) - A.cast<std::complex<double>>(),
     B.cast<std::complex<double>>();
Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);
// Last m columns of V are the null space basis

// For eigenvector selection (eigenvalue assignment without eigenvector preference):
// Choose the null space vector that maximizes the angle to previously chosen eigenvectors
// (KNV robustness heuristic)
```

### Implementation Structure

```cpp
// src/analysis/eigen_assignment.h
#pragma once

#include "../linear_system.h"
#include <Eigen/Core>
#include <complex>
#include <vector>
#include <string>

namespace caliburn {

struct EigenAssignResult {
    Eigen::MatrixXd K;                                // m × n gain matrix
    std::vector<std::complex<double>> achieved_poles;  // actual eig(A-BK)
    Eigen::MatrixXcd eigenvectors;                     // achieved eigenvectors
    bool success;
    std::string error;
};

// MIMO eigenvalue assignment using KNV method
// desired_poles: n complex numbers (conjugate pairs for real K)
// desired_vectors: optional partial eigenvector directions (empty = auto)
EigenAssignResult eigenAssign(
    const LinearSystem& sys,
    const std::vector<std::complex<double>>& desired_poles,
    const std::vector<Eigen::VectorXcd>& desired_vectors = {});

}  // namespace caliburn
```

### Relation to Ackermann (Tier 2)

Eigenvalue assignment subsumes Ackermann:
- For SISO (m=1): the null space has dimension 1, so there is no eigenvector freedom. The KNV method produces the same K as Ackermann.
- For MIMO (m>1): the null space has dimension m, giving freedom in eigenvector selection. KNV resolves this freedom optimally.

Option: Replace Ackermann entirely with KNV (it handles both SISO and MIMO). Or keep Ackermann for SISO as a simpler/faster path and use KNV only for MIMO.

### Controllability Requirement

Same as Ackermann: `rank([B, AB, ..., A^(n-1)B]) = n`. Use existing `checkControllability()`.

## Eigenstructure Assignment (Phase 2)

### Problem Statement

Same as eigenvalue assignment, but additionally specify desired eigenvector directions. The user wants not just WHERE the poles are, but HOW the modes couple to the states.

Example: In a 2-DOF system (x, θ), the user might want:
- Mode 1 (λ_1): mostly position (eigenvector ≈ [1, 0, *, *])
- Mode 2 (λ_2): mostly angle (eigenvector ≈ [0, 1, *, *])

### Algorithm Extension

The KNV null space gives m candidate directions for each eigenvector. With a desired direction d_i:

1. Compute null space N_i of `[λ_i·I - A | B]` (dimension m)
2. Extract the top n rows of each null space vector to get achievable eigenvector directions
3. Project d_i onto the space spanned by these achievable directions
4. Choose the achievable eigenvector closest to d_i (minimize angle)

The achieved eigenvector will exactly match d_i only if d_i lies within the achievable subspace. Otherwise, it's the closest achievable approximation.

### UI Challenge

Eigenvector specification requires:
- n eigenvectors, each with n complex components
- User needs to understand what each component means (which state)
- Conjugate pairing for complex modes

This is significantly more complex than just specifying pole locations.

**Recommended approach**:
- Start with eigenvalue assignment only (no eigenvector preference)
- Add eigenvector UI later as a "Mode Shape" section:
  - Show achieved eigenvectors as bar charts (magnitude per state)
  - Let user drag bars to specify desired directions
  - Real-time update showing achievable vs desired

### Deferral

Eigenstructure assignment uses the same `eigenAssign()` function — it just passes non-empty `desired_vectors`. The algorithm doesn't change; only the UI needs work. This makes it safe to implement the algorithm now and add the UI later.

## AppState Additions

```cpp
// In AppState (extending Tier 2's additions):
// desired_poles already exists from PolePlacement (Tier 2)
// ctrl_K already exists from GainMatrix

// Eigenvalue assignment reuses desired_poles and ctrl_K
// No new AppState fields needed for Phase 1

// Phase 2 (eigenstructure) would add:
// std::vector<Eigen::VectorXcd> desired_eigenvectors;
```

## ControllerType Extension

```cpp
enum class ControllerType {
    None,
    PID,               // Tier 1
    LeadLag,           // Tier 1
    PolePlacement,     // Tier 2 (SISO, Ackermann)
    LQR,               // Tier 2
    EigenAssignment,   // Tier 3 (MIMO, KNV)
    StateSpace,        // Existing
    GainMatrix,        // Existing
};
```

### UI Routing Logic

```
if (plant.inputs() == 1) {
    // SISO: show PolePlacement (Ackermann) — simpler, direct
} else {
    // MIMO: show EigenAssignment (KNV) — handles multi-input
}
```

Or let the user choose either. For SISO plants, EigenAssignment degenerates to the same result as PolePlacement.

## UI Layout

```
Eigenvalue Assignment
├── ⚠ "Not controllable" (if checkControllability fails)
├── Desired poles:
│   ├── λ1: Re [===] Im [===]    (with conjugate auto-pairing)
│   ├── λ2: Re [===] Im [===]
│   └── ...
├── [Compute K]
├── K = [matrix display]
├── Achieved poles: λ1=..., λ2=..., ...
├── Achieved eigenvectors: [matrix display or bar chart]  ← Phase 2
└── Mode shapes: [visual display]  ← Phase 2
```

The desired poles UI is identical to Tier 2's PolePlacement. The difference is the algorithm (KNV instead of Ackermann) and that it works for MIMO plants.

## Files Changed

| File | Change |
|---|---|
| `src/app_state.h` | Add EigenAssignment to ControllerType enum. |
| `src/analysis/eigen_assignment.h` | New file. EigenAssignResult struct, eigenAssign() declaration. |
| `src/analysis/eigen_assignment.cpp` | New file. KNV algorithm with SVD-based null space computation. |
| `src/panels/model_panel.cpp` | Add EigenAssignment UI section (reuse desired poles UI from Tier 2). |
| `src/visualizer.cpp` | Add EigenAssignment case to recompute loop (same path as GainMatrix). |
| `CMakeLists.txt` | Add eigen_assignment.cpp to analysis_lib. |
| `tests/test_eigen_assignment.cpp` | New test: SISO matches Ackermann, MIMO places poles correctly. |

## Test Cases

### SISO Verification (should match Ackermann)

```
Plant: Second-Order (A=[0 1; -1 -1.4], B=[0; 1])
Desired poles: [-2+j, -2-j]
Expected K: same as Ackermann result from Tier 2
```

### MIMO Verification

```
Plant: Ball-Balancer (n=4, m=2)
Desired poles: [-1, -2, -3, -4]
Verify: eig(A - BK) ≈ {-1, -2, -3, -4}
Verify: K is 2×4 real matrix
```

### Robustness

```
Plant: Double Mass-Spring (n=4, m=2)
Desired poles: [-1±j, -2±j]
Verify: complex conjugate poles produce real K
Verify: achieved poles match desired within 1e-6
```

## Verification

1. Build: `cmake --build build` — clean compile
2. Select Ball-Balancer (MIMO), Eigenvalue Assignment — specify 4 desired poles
3. Compute K — 2×4 gain matrix displayed
4. Pole-zero plot shows poles at desired locations
5. Step response shows stable, well-damped behavior
6. Compare SISO case with Ackermann — same K produced
7. Change desired poles interactively — K and all plots update
8. All existing controller types still work unchanged

## Future: Observer Integration

Both Tier 2 and Tier 3 state feedback methods assume full-state access. A future tier would add:

1. **Luenberger Observer**: Choose observer poles (same eigenAssign algorithm applied to `(A', C')`)
2. **Observer gain L**: computed by duality — `eigenAssign` on transposed system, then transpose result
3. **Augmented closed-loop**: 2n states (plant + observer error dynamics)
4. **Separation principle**: Controller and observer designed independently, combined system is stable if both are individually stable

The observer uses the same eigenvalue assignment algorithm — just applied to the dual system. The `eigenAssign()` function from this tier is directly reusable.
