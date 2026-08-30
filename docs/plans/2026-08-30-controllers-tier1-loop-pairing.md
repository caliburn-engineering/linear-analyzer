# Tier 1 Controllers on a Loop-Pairing Model — Implementation Plan

> **For agentic workers:** Implement this plan task-by-task, in order. Steps use checkbox (`- [ ]`) syntax for tracking.
>
> *Status: not started.*

**Goal:** Ship PID and Lead/Lag controllers built on a **loop-pairing** data model that works for square, SIMO and MISO plants — with the pairing diagnostic (RGA / Channel Share), rule-C channel indexing, and a κ-on-compensator root locus that responds to the tuning.

**Architecture:** The controller is a `std::vector<Loop>`, each loop pairing one plant output to one plant input and carrying both a `PIDParams` and a `LeadLagParams` block. A new `analysis_lib` module assembles the m-out/p-in controller from that list; the existing `seriesConnect` → `feedbackConnect` pipeline is unchanged. The model panel's controller section becomes a **pairing grid** — the p×m diagnostic matrix *is* the editor.

**Tech Stack:** C++17, Eigen 3.4, ImGui (docking), ImPlot, CMake with FetchContent. Desktop + Emscripten.

**Provenance:** Every decision below is settled on [map #1](https://github.com/caliburn-engineering/caliburn/issues/1) and its closed tickets — [#2 loop pairing data model](https://github.com/caliburn-engineering/caliburn/issues/2), [#3 builder API](https://github.com/caliburn-engineering/caliburn/issues/3), [#4 RGA](https://github.com/caliburn-engineering/caliburn/issues/4), [#5 zero-state audit](https://github.com/caliburn-engineering/caliburn/issues/5), [#6 controller UI](https://github.com/caliburn-engineering/caliburn/issues/6), [#8 Channel Share](https://github.com/caliburn-engineering/caliburn/issues/8), [#9 rule C](https://github.com/caliburn-engineering/caliburn/issues/9), [#10 loop locus](https://github.com/caliburn-engineering/caliburn/issues/10). **Do not re-open them.** Where this plan makes a call no ticket covered, it is marked **[plan-level]** with the reasoning inline.

This plan supersedes the MIMO half of `docs/specs/2026-07-23-controllers-tier1-pid-leadlag-design.md`, whose per-channel diagonal model assumed square plants. The spec's **state-space realization tables are correct and reused verbatim**; its `std::vector<PIDParams>`-sized-by-`plant.inputs()` model, its `buildDiagonalPID` API, its DC-gain coupling heuristic and its verification list are not.

---

## Read this first — invariants

Nine facts hold across every task. Contradicting one means you have gone off the plan.

1. **The controller is m-out/p-in.** For a plant with `p` outputs and `m` inputs, the controller reads `p` errors and drives `m` inputs, so `seriesConnect(C, G)` is `p×p` and `feedbackConnect` needs no change. Unpaired reference channels are **inert** (a zero column in `C`), not illegal.
2. **Minimal realization is a correctness requirement, not tidiness.** A state whose output coefficient is zero is unobservable, survives `seriesConnect`/`feedbackConnect` untouched, and appears as a **phantom pole** in the closed-loop pole plot. Emit a state only when its output coefficient is nonzero, tested at **exact zero, no epsilon**.
3. **n = 0 controllers are legal.** A P-only loop contributes zero states. The combined path is numerically exact; only *standalone*-controller analysis faults, and Task 1's three guards close that.
4. **Rule C — three indices.** `output_i ∈ [0,p)`, `input_j ∈ [0,m)`, `ref_r ∈ [0,p)`. Plant is indexed `(i,j)`; controller `(j,i)` **transposed**; open/closed loop `(i,r)`. Never apply a plant-shaped index pair to a system of another shape.
5. **Liveness is read off the loop list, never off the matrices.** Controller channel `(j,i)` is live ⟺ `∃k: out_k == i ∧ in_k == j`. Reference `r` is live ⟺ `∃k: out_k == r`. Dead channels are **suppressed with a stated reason** — the existing fallbacks lie (`arg(0)` draws a flat 0° phase; the pole-zero pencil returns `eig(A)` as coincident zeros, reading as total cancellation).
6. **RGA uses `.transpose()`, never `.adjoint()`.** `Λ = G ∘ (G⁻¹)ᵀ` on the **complex** `G(jω)`. No real-matrix test can catch this bug — Task 6's row-sum test is what catches it.
7. **The RGA exists only on a square paired sub-matrix of size ≥ 2.** Everywhere else — single-input plant, multi-input plant with one loop, MISO fan-out — the readout is **Channel Share**, never labelled or coloured as RGA.
8. **κ is never 0.** `κ ∈ [0.01, 100]`, log-spaced. κ = 0 would zero a loop's gains, the minimal realization would then emit no states for it, the closed-loop dimension would shrink mid-sweep, and `matchPoles` would read `prev[i]` out of bounds. Constant dimension makes both `matchPoles` and the panel's drawing loop safe **by construction**.
9. **Hard sequencing:** the three `n == 0` guards and the `CMAKE_BUILD_TYPE` flip **land together, first** (Task 1). The flip makes `NDEBUG` live, and under `NDEBUG` `Eigen::EigenSolver` on a 0×0 matrix *segfaults* rather than asserting. `tests/test_controller_builders.cpp` also cannot pass in a debug build until guard 3 exists.

### Starting point

Branch `docs/reorganize` of the **linear-analyzer** repo (`projects/linear-analyzer`, its own git repo — `projects/` is gitignored in caliburn). All line numbers below are against that branch's **working tree**, which carries uncommitted changes to `CMakeLists.txt` and five panel files (the Emscripten/web build). The tickets were resolved against this tree; commit or stash-and-reapply it before starting, but do not rebase the plan onto `HEAD` — the line numbers will not match.

### Facts established during charting — do not re-derive

- **Preset shapes.** Six presets: `0` First-Order (1×1, n=1), `1` Second-Order (1×1, n=2), `2` Ball-Balancer (**2×2**, n=4), `3` Inverted Pendulum (1-in/2-out, n=4), `4` Quarter-Car (1-in/2-out, n=4), `5` Double Mass-Spring (1-in/2-out, n=4). **Ball-Balancer is the only genuinely multivariable preset.** Quarter-Car is SIMO, which is why the spec's verification step 8 is unachievable. *(The map's Notes say "2 of 5 presets are square"; the accurate count is 3 of 6 — Second-Order is square too. The substantive claim is unchanged.)*
- **Ball-Balancer's `A` is pure integrators — singular.** `−C·A⁻¹·B + D` is undefined for it; its DC gain is genuinely infinite. Its coupling is a channel **swap**: `B(2,1)` and `B(3,0)` mean input 0 drives output 1's chain. `G = (5g/7)/s² · [[0,1],[1,0]]`, so `Λ = [[0,1],[1,0]]` exactly at every ω > 0, `cond(G) = 1` throughout. RGA number **4.00** for the seeded identity pairing, **0.00** for the swap.
- **Build cost.** `CMAKE_BUILD_TYPE` is empty and `CMakeLists.txt` sets no flags, so all 15 of the app's own TUs compile with no `-O` and no `NDEBUG`. A 200-point locus sweep at 10 closed-loop states costs **321 ms** as built today versus **2.4 ms** at `-O2 -DNDEBUG`.
- **Layering.** `app_state.h:11` includes `imgui.h`; `analysis_lib` and all test targets link only Eigen. Nothing under `src/analysis/` may include `app_state.h`. This is why the parameter structs move *down* into the analysis header rather than the builders reaching *up*.
- **Emscripten.** `analysis_lib` is not gated on `if(NOT EMSCRIPTEN)`, so every new analysis source joins the web build automatically. The panel/UI side is unverified and stays out of scope.

---

### Task 1: Zero-state guards and the Release default

**Files:**
- Modify: `src/analysis/pole_zero.cpp`
- Modify: `src/analysis/pole_zero.h`
- Modify: `src/analysis/frequency_response.cpp`
- Modify: `CMakeLists.txt`

Nothing else may land before this task. The guards are provably behaviour-preserving for every `n > 0` case — the audit measured a worst delta of exactly 0.0 against current behaviour, with both zero branches exercised and `is_stable` unchanged.

- [ ] **Step 1: Guard the poles `EigenSolver` in `computePoleZero`**

In `src/analysis/pole_zero.cpp`, replace lines 38–47 (the "Poles = eigenvalues of A" block) with:

```cpp
    // Poles = eigenvalues of A.
    // n == 0 is a static gain: no poles, trivially stable.  Eigen::EigenSolver
    // on a 0x0 matrix is a real memory fault in BOTH build configs (SIGSEGV
    // under NDEBUG, SIGABRT in debug) — not a disabled sanity check — so this
    // guard is a correctness requirement.  See issue #5.
    result.is_stable = true;
    if (n > 0) {
        Eigen::EigenSolver<Eigen::MatrixXd> es(sys.A, false);
        result.poles.resize(n);
        for (int i = 0; i < n; ++i) {
            result.poles[i] = es.eigenvalues()(i);
            if (result.poles[i].real() >= -1e-10) {
                result.is_stable = false;
            }
        }
    }
```

- [ ] **Step 2: Guard the zeros `EigenSolver` in `computePoleZero`**

Still in `computePoleZero`, wrap the `|D_ij| > 1e-14` branch body (lines 52–59). Change:

```cpp
    if (std::abs(D_ij) > 1e-14) {
        Eigen::MatrixXd A_z = sys.A -
```

to:

```cpp
    if (std::abs(D_ij) > 1e-14) {
      // A static gain has no zeros; leave result.zeros empty rather than
      // pushing a 0x0 matrix into EigenSolver.  See issue #5.
      if (n > 0) {
        Eigen::MatrixXd A_z = sys.A -
```

and close the new scope after the `result.zeros[i] = zes.eigenvalues()(i);` loop, before the `} else {`.

**Leave the `else` (pencil) branch alone.** It is already correct at `n = 0`: `P` and `Q` degenerate to 1×1 with `Q = [0]`, `GeneralizedEigenSolver` handles the singular pencil, `beta = 0`, and the `|beta| > 1e-10` filter drops it — an empty zero set, which is right for a static gain.

- [ ] **Step 3: Early-return in `evalTransferFunction`**

In `src/analysis/frequency_response.cpp`, insert immediately after `int n = sys.states();` (line 11):

```cpp
    // n == 0 is a static gain.  ColPivHouseholderQR asserts on a 0x0 matrix in
    // a debug build; returning D directly is also exactly right.  See issue #5.
    if (n == 0) {
        return {sys.D(output_i, input_j), 0.0};
    }
```

**Do not swap the decomposition instead.** `partialPivLu` handles 0×0 cleanly and looks tidier, but it returns **NaN** at a pole and **−1e24** at 1e−12 off one, where `colPivHouseholderQr` degrades gracefully to a finite value. On a pure-integrator plant — the Ball-Balancer shape, poles at the origin — that is a live case.

- [ ] **Step 4: Guard both locus functions**

Belt-and-braces: neither is reachable with `n == 0` today (both run on the plant or the closed loop, which always carry the plant's states), but the guard costs one branch and the panel already handles an empty locus. In `computeRootLocus` and `computeStateFeedbackLocus`, insert after `int n = sys.states();`:

```cpp
    if (n == 0) return {};
```

- [ ] **Step 5: Expose `matchPoles`** **[plan-level]**

`computeLoopLocus` (Task 7) needs the same pole-continuity matching, and duplicating it would let the two copies drift. Move `matchPoles` out of the anonymous namespace in `pole_zero.cpp` into the `caliburn` namespace, and declare it in `pole_zero.h` above `computePoleZero`:

```cpp
// Reorder `poles` so each entry continues the nearest entry of `prev`.
// Requires poles.size() == prev.size() — a locus whose pole count changes
// between points would read prev[i] out of bounds.
void matchPoles(std::vector<std::complex<double>>& poles,
                const std::vector<std::complex<double>>& prev);
```

This keeps `pole_zero` a general module over any `LinearSystem` with no knowledge of controllers — the property #5's audit and #9's rule C both reason about.

- [ ] **Step 6: Default to a Release build**

In `CMakeLists.txt`, insert after `set(CMAKE_EXPORT_COMPILE_COMMANDS ON)` (line 7):

```cmake
# Default to an optimized build.  Measured (#10): a 200-point root-locus sweep
# at 10 closed-loop states costs 321 ms as this repo built before this line,
# versus 2.4 ms at -O2 -DNDEBUG — and needs_recompute is set from 26 sites
# including every slider-drag frame.  This makes NDEBUG live, which is only
# safe once the n == 0 guards above exist, so the two land together.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
```

An explicit `-DCMAKE_BUILD_TYPE=Debug` still wins — the `if` does not fire when the variable is already set.

- [ ] **Step 7: Reconfigure and verify**

The existing `build/` directory has an empty `CMAKE_BUILD_TYPE` cached, so it must be reconfigured:

```bash
cmake -S . -B build
cmake --build build -j8
ctest --test-dir build --output-on-failure
```

Confirm `grep CMAKE_BUILD_TYPE build/CMakeCache.txt` now reads `Release`, and that all six existing suites pass. Tests are unaffected by `NDEBUG` — `tests/test_helpers.h` uses `std::exit(1)`, not `assert()`.

- [ ] **Step 8: Commit**

```bash
git add src/analysis/pole_zero.cpp src/analysis/pole_zero.h \
        src/analysis/frequency_response.cpp CMakeLists.txt
git commit -m "fix: guard analysis primitives against 0-state systems; default to Release

Eigen::EigenSolver on a 0x0 matrix segfaults in both build configs and
ColPivHouseholderQR asserts in debug, so a 0-state controller could not be
analyzed standalone. Three guards close that. The Release default is folded
in here because it makes NDEBUG live, which is only safe behind the guards.

Refs #5, #10"
```

---

### Task 2: Controller builders module and its tests

**Files:**
- Create: `src/analysis/controller_builders.h`
- Create: `src/analysis/controller_builders.cpp`
- Create: `tests/test_controller_builders.cpp`
- Modify: `CMakeLists.txt`

No UI, no `AppState`. This task is self-contained and testable.

- [ ] **Step 1: Create `src/analysis/controller_builders.h`**

```cpp
// src/analysis/controller_builders.h
#pragma once

#include "../linear_system.h"
#include <Eigen/Core>
#include <vector>

namespace caliburn {

// Parameter blocks live here, not in app_state.h: app_state.h includes imgui.h,
// which analysis_lib and every test target cannot see.  app_state.h already
// pulls in six analysis headers, so it takes a seventh.  Fields stay `float`
// so ImGui sliders bind them directly; the builders widen to double.

struct PIDParams {
    float Kp = 1.0f;
    float Ki = 0.0f;
    float Kd = 0.0f;
    float tau_f = 0.01f;  // derivative filter time constant.  Precondition: > 0.
};

enum class CompensatorMode { Lead, Lag, LeadLag };

struct LeadLagParams {
    CompensatorMode mode = CompensatorMode::Lead;
    float Kc = 1.0f;
    float alpha_lead = 0.1f;  // < 1.  Precondition: > 0.
    float T_lead = 1.0f;      // Precondition: > 0.
    float alpha_lag = 10.0f;  // > 1.  Precondition: > 0.
    float T_lag = 1.0f;       // Precondition: > 0.
};

// One control loop: reads e[out] = r[out] - y[out], drives u[in].
// Both parameter blocks are carried inline and LoopKind selects which is read,
// so switching kind loses neither tuning nor the pairing — the pairing being
// the expensive decision.  See issue #2 for the rejected alternatives.
struct Loop {
    int out = 0;
    int in = 0;
    PIDParams pid{};
    LeadLagParams leadlag{};
};

// Not ControllerType: that enum lives in app_state.h and carries StateSpace and
// GainMatrix, which mean nothing to the assembler.
enum class LoopKind { PID, LeadLag };

// C(s) = Kp + Ki/s + Kd*s/(tau_f*s + 1), minimal realization.
LinearSystem buildPID(const PIDParams& p);

// C(s) = Kc * (s + 1/T) / (s + 1/(alpha*T)), minimal realization; the cascade
// mode is the two sections in series.  See the note below on the constant.
LinearSystem buildLeadLag(const LeadLagParams& p);

// Assemble the m-out / p-in controller for a loop list against `plant`.
// The plant is passed rather than (int p, int m): two same-typed ints transpose
// silently and build a wrongly-shaped controller, and this matches the existing
// stateFeedbackClose(const LinearSystem&, const MatrixXd&).
LinearSystem buildLoopController(const std::vector<Loop>& loops,
                                 LoopKind kind,
                                 const LinearSystem& plant);

}  // namespace caliburn
```

> **Note on the Lead/Lag constant.** The spec writes `C(s) = Kc·(Ts + 1)/(αTs + 1)` at line 193 but gives the realization `A = [−1/(αT)], B = [1], C = [Kc(α−1)/(αT)], D = [Kc]`, whose transfer function is `Kc·(s + 1/T)/(s + 1/(αT))` — the two differ by a factor of α, which `Kc` absorbs. The **realization is what was verified**, so it is what ships, and the pole-zero form above is the contract the tests assert. Both forms agree at α = 1 (pure gain `Kc`).

- [ ] **Step 2: Create `src/analysis/controller_builders.cpp`**

```cpp
// src/analysis/controller_builders.cpp
#include "controller_builders.h"
#include "system_connect.h"

namespace caliburn {

LinearSystem buildPID(const PIDParams& p) {
    // Emit a state only when its output coefficient is nonzero.  A state with a
    // zero entry in C is unobservable, survives seriesConnect and
    // feedbackConnect untouched, and shows up as a phantom pole in the
    // closed-loop pole plot: PI with Ki = 0 puts a marginally-stable pole at
    // the origin, PD with Kd = 0 puts one at -1/tau_f.  Emitting the state does
    // not waste space, it lies to the plots.
    //
    // The threshold is exact zero, no epsilon.  The integrator's A block is [0]
    // regardless of Ki (Ki enters only C), so a tiny Ki is numerically harmless
    // and genuinely is a near-integrator; a dead band would silently contradict
    // a slider reading 1e-9.  The left slider stop writes exact 0.0f.
    //
    // No defensive guard on tau_f <= 0: it is a documented precondition, the
    // sliders enforce the range, and a silent clamp here would hide a
    // programming error at the one place it is cheap to catch.
    const bool has_i = (p.Ki != 0.0f);
    const bool has_d = (p.Kd != 0.0f);
    const int n = (has_i ? 1 : 0) + (has_d ? 1 : 0);

    const double Kp = p.Kp, Ki = p.Ki, Kd = p.Kd, tau_f = p.tau_f;

    LinearSystem c;
    c.A = Eigen::MatrixXd::Zero(n, n);
    c.B = Eigen::MatrixXd::Zero(n, 1);
    c.C = Eigen::MatrixXd::Zero(1, n);
    c.D = Eigen::MatrixXd::Constant(1, 1, Kp);

    int k = 0;
    if (has_i) {  // integrator: Ki / s
        c.A(k, k) = 0.0;
        c.B(k, 0) = 1.0;
        c.C(0, k) = Ki;
        ++k;
    }
    if (has_d) {  // filtered derivative: Kd*s / (tau_f*s + 1)
        c.A(k, k) = -1.0 / tau_f;
        c.B(k, 0) = 1.0 / tau_f;
        c.C(0, k) = -Kd / tau_f;
        c.D(0, 0) += Kd / tau_f;
    }
    return c;
}

namespace {

// One first-order section (s + z) / (s + p), z = 1/T, p = 1/(alpha*T).
// Same observability rule as buildPID: the output coefficient is z - p, which
// vanishes at alpha == 1, so that case drops to a pure gain with no state.
// Unreachable through the sliders (alpha_lead is clamped to [0.01, 0.99],
// alpha_lag to [1.01, 100]) — handled anyway, because a builder that is total
// in alpha does not depend on slider clamps staying where they are.
LinearSystem leadLagSection(double alpha, double T) {
    const double z = 1.0 / T;
    const double pole = 1.0 / (alpha * T);
    const double coeff = z - pole;
    const int n = (coeff != 0.0) ? 1 : 0;

    LinearSystem s;
    s.A = Eigen::MatrixXd::Zero(n, n);
    s.B = Eigen::MatrixXd::Zero(n, 1);
    s.C = Eigen::MatrixXd::Zero(1, n);
    s.D = Eigen::MatrixXd::Constant(1, 1, 1.0);
    if (n == 1) {
        s.A(0, 0) = -pole;
        s.B(0, 0) = 1.0;
        s.C(0, 0) = coeff;
    }
    return s;
}

}  // anonymous namespace

LinearSystem buildLeadLag(const LeadLagParams& p) {
    // Kc == 0 zeroes every output coefficient at once, so the same rule says
    // emit no states.  Also unreachable through the sliders (Kc >= 0.1).
    if (p.Kc == 0.0f) {
        LinearSystem z;
        z.A = Eigen::MatrixXd::Zero(0, 0);
        z.B = Eigen::MatrixXd::Zero(0, 1);
        z.C = Eigen::MatrixXd::Zero(1, 0);
        z.D = Eigen::MatrixXd::Zero(1, 1);
        return z;
    }

    LinearSystem c;
    switch (p.mode) {
        case CompensatorMode::Lead:
            c = leadLagSection(p.alpha_lead, p.T_lead);
            break;
        case CompensatorMode::Lag:
            c = leadLagSection(p.alpha_lag, p.T_lag);
            break;
        case CompensatorMode::LeadLag:
            // Each section drops independently: alpha_lead == 1 gives 1 state,
            // both == 1 gives 0.  seriesConnect handles n1 == 0 exactly (#5).
            c = seriesConnect(leadLagSection(p.alpha_lead, p.T_lead),
                              leadLagSection(p.alpha_lag, p.T_lag));
            break;
    }
    c.C *= static_cast<double>(p.Kc);
    c.D *= static_cast<double>(p.Kc);
    return c;
}

LinearSystem buildLoopController(const std::vector<Loop>& loops,
                                 LoopKind kind,
                                 const LinearSystem& plant) {
    const int p = plant.outputs();
    const int m = plant.inputs();

    // Build each loop's SISO block first — the total state count is not known
    // until every block has decided how minimal it is.  The assembler calls the
    // SISO builders rather than inlining the realization tables a second time,
    // so the observability rule lives in exactly one place.
    std::vector<LinearSystem> blocks;
    std::vector<const Loop*> kept;
    blocks.reserve(loops.size());
    kept.reserve(loops.size());
    int n_total = 0;

    for (const auto& l : loops) {
        // Out-of-range index: skip this loop and build the rest.  Unreachable
        // in principle — the (p, m) dimension guard reseeds wholesale before
        // any recompute — so this only decides what happens when it occurs
        // anyway.  Asserting would crash a UI that must keep rendering while
        // shipped builds scattered out of bounds; clamping would silently
        // control a different physical channel while carrying the old
        // channel's gains; bailing to all-zero would discard every other
        // loop's tuning.  Skipping degrades locally and behaves identically in
        // both build configs.
        if (l.out < 0 || l.out >= p || l.in < 0 || l.in >= m) continue;
        blocks.push_back(kind == LoopKind::PID ? buildPID(l.pid)
                                               : buildLeadLag(l.leadlag));
        kept.push_back(&l);
        n_total += blocks.back().states();
    }

    // Block-diagonal A; B row-block k scattered into column out_k; C
    // column-block k into row in_k; D accumulated at (in_k, out_k).  Only D
    // accumulates — the B and C scatters are unambiguous because each state
    // block belongs to exactly one loop.  An empty (or entirely skipped) list
    // yields the all-zero m x p system with states() == 0, which is issue #2's
    // "empty loop list is equivalent to ControllerType::None".
    LinearSystem c;
    c.A = Eigen::MatrixXd::Zero(n_total, n_total);
    c.B = Eigen::MatrixXd::Zero(n_total, p);
    c.C = Eigen::MatrixXd::Zero(m, n_total);
    c.D = Eigen::MatrixXd::Zero(m, p);

    int off = 0;
    for (std::size_t k = 0; k < blocks.size(); ++k) {
        const LinearSystem& b = blocks[k];
        const int nk = b.states();
        const int out = kept[k]->out;
        const int in = kept[k]->in;

        if (nk > 0) {
            c.A.block(off, off, nk, nk) = b.A;
            c.B.block(off, out, nk, 1) = b.B;
            c.C.block(in, off, 1, nk) = b.C;
        }
        c.D(in, out) += b.D(0, 0);
        off += nk;
    }
    return c;
}

}  // namespace caliburn
```

- [ ] **Step 3: Create `tests/test_controller_builders.cpp`**

Two oracles: `evalTransferFunction` for *what the realization computes*, and `states()` for *how minimally*. The transfer function cannot see an unobservable mode — a PI built with `Ki = 0` matches its transfer function perfectly — so the state count needs its own direct assertion.

```cpp
// tests/test_controller_builders.cpp
#include "analysis/controller_builders.h"
#include "analysis/frequency_response.h"
#include "analysis/system_connect.h"
#include "test_helpers.h"

#include <complex>
#include <cstdio>
#include <vector>

using namespace caliburn;
using cd = std::complex<double>;

namespace {

// Evaluate off the poles: a PI block has a pole at s = 0, so the DC check used
// by test_system_connect.cpp is invalid here.  This set spans below, at and
// above the derivative-filter corner (tau_f = 0.1 -> 10 rad/s).
const cd kS[3] = {cd(0.0, 0.1), cd(0.0, 1.0), cd(0.0, 10.0)};
constexpr double kTol = 1e-10;  // matches the algebraic-identity tests

void assertChannel(const LinearSystem& sys, int i, int j, cd expected, cd s) {
    const cd actual = evalTransferFunction(sys, i, j, s);
    ASSERT_NEAR(actual.real(), expected.real(), kTol);
    ASSERT_NEAR(actual.imag(), expected.imag(), kTol);
}

cd pidTF(const PIDParams& p, cd s) {
    return cd(p.Kp, 0.0) + cd(p.Ki, 0.0) / s
         + cd(p.Kd, 0.0) * s / (cd(p.tau_f, 0.0) * s + 1.0);
}

cd sectionTF(double alpha, double T, cd s) {
    return (s + 1.0 / T) / (s + 1.0 / (alpha * T));
}

cd leadLagTF(const LeadLagParams& p, cd s) {
    cd v(p.Kc, 0.0);
    switch (p.mode) {
        case CompensatorMode::Lead:
            v *= sectionTF(p.alpha_lead, p.T_lead, s); break;
        case CompensatorMode::Lag:
            v *= sectionTF(p.alpha_lag, p.T_lag, s); break;
        case CompensatorMode::LeadLag:
            v *= sectionTF(p.alpha_lead, p.T_lead, s)
               * sectionTF(p.alpha_lag, p.T_lag, s); break;
    }
    return v;
}

// Hand-built shape fixtures, not getBuiltinModels().  Preset numerics belong to
// test_model_library.cpp; coupling these tests to them means a physical-
// parameter tweak breaks an unrelated suite.  Shape is what the assembler
// cares about.
LinearSystem makePlant(int p, int m) {
    LinearSystem g;
    g.A = Eigen::MatrixXd::Constant(1, 1, -1.0);
    g.B = Eigen::MatrixXd::Ones(1, m);
    g.C = Eigen::MatrixXd::Ones(p, 1);
    g.D = Eigen::MatrixXd::Zero(p, m);
    return g;
}

// --------------------------------------------------------------- SISO builders

void test_pid_transfer_functions() {
    PIDParams p_only{2.5f, 0.0f, 0.0f, 0.1f};
    PIDParams pi{1.0f, 3.0f, 0.0f, 0.1f};
    PIDParams pd{1.0f, 0.0f, 0.25f, 0.1f};
    PIDParams pid{1.5f, 3.0f, 0.25f, 0.1f};

    for (const auto& p : {p_only, pi, pd, pid}) {
        const LinearSystem c = buildPID(p);
        ASSERT_EQ(c.outputs(), 1);
        ASSERT_EQ(c.inputs(), 1);
        for (const cd& s : kS) assertChannel(c, 0, 0, pidTF(p, s), s);
    }
}

void test_pid_minimal_realization() {
    // The rule the transfer-function oracle cannot see.
    ASSERT_EQ(buildPID(PIDParams{2.5f, 0.0f, 0.0f, 0.1f}).states(), 0);
    ASSERT_EQ(buildPID(PIDParams{1.0f, 3.0f, 0.0f, 0.1f}).states(), 1);
    ASSERT_EQ(buildPID(PIDParams{1.0f, 0.0f, 0.25f, 0.1f}).states(), 1);
    ASSERT_EQ(buildPID(PIDParams{1.5f, 3.0f, 0.25f, 0.1f}).states(), 2);
}

void test_leadlag_transfer_functions() {
    LeadLagParams lead{CompensatorMode::Lead, 2.0f, 0.1f, 1.0f, 10.0f, 1.0f};
    LeadLagParams lag{CompensatorMode::Lag, 2.0f, 0.1f, 1.0f, 10.0f, 0.5f};
    LeadLagParams both{CompensatorMode::LeadLag, 1.5f, 0.2f, 0.3f, 5.0f, 2.0f};

    for (const auto& p : {lead, lag, both}) {
        const LinearSystem c = buildLeadLag(p);
        for (const cd& s : kS) assertChannel(c, 0, 0, leadLagTF(p, s), s);
    }

    ASSERT_EQ(buildLeadLag(lead).states(), 1);
    ASSERT_EQ(buildLeadLag(lag).states(), 1);
    ASSERT_EQ(buildLeadLag(both).states(), 2);
}

void test_leadlag_degenerate_alpha() {
    // alpha == 1 zeroes the output coefficient -> pure gain, no state.
    LeadLagParams unity{CompensatorMode::Lead, 3.0f, 1.0f, 1.0f, 10.0f, 1.0f};
    ASSERT_EQ(buildLeadLag(unity).states(), 0);
    for (const cd& s : kS) assertChannel(buildLeadLag(unity), 0, 0, cd(3.0, 0.0), s);

    LeadLagParams half{CompensatorMode::LeadLag, 1.0f, 1.0f, 1.0f, 5.0f, 2.0f};
    ASSERT_EQ(buildLeadLag(half).states(), 1);

    LeadLagParams none{CompensatorMode::LeadLag, 1.0f, 1.0f, 1.0f, 1.0f, 2.0f};
    ASSERT_EQ(buildLeadLag(none).states(), 0);
}

// ------------------------------------------------------------------ assembler

void test_assembler_shape_and_placement() {
    const LinearSystem plant = makePlant(2, 2);
    PIDParams a{1.0f, 2.0f, 0.0f, 0.1f};   // 1 state
    PIDParams b{3.0f, 0.0f, 0.5f, 0.1f};   // 1 state
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{1, 1, b, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), plant.inputs());   // m
    ASSERT_EQ(c.inputs(), plant.outputs());   // p
    ASSERT_EQ(c.states(), 2);

    // Loop k occupies channel (in_k, out_k) — the controller is transposed
    // relative to the plant.
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 1, 1, pidTF(b, s), s);
        assertChannel(c, 0, 1, cd(0.0, 0.0), s);
        assertChannel(c, 1, 0, cd(0.0, 0.0), s);
    }
}

void test_assembler_unpaired_channel_is_inert() {
    // SIMO: 1 input, 2 outputs, one loop y0 -> u0.  Channel (0,1) of the 1x2
    // controller must be exactly zero at every frequency — this is the "inert
    // column" the whole loop-pairing model rests on, and nothing else tests it.
    const LinearSystem plant = makePlant(2, 1);
    PIDParams a{2.0f, 1.0f, 0.0f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), 1);
    ASSERT_EQ(c.inputs(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 0, 1, cd(0.0, 0.0), s);
    }
}

void test_assembler_miso_fan_in_sums() {
    // Two loops into one input — the primary SIMO control structure, e.g. the
    // Inverted Pendulum's u = PID_angle(e1) + PID_pos(e0).
    const LinearSystem plant = makePlant(2, 1);
    PIDParams a{1.0f, 2.0f, 0.0f, 0.1f};
    PIDParams b{4.0f, 0.0f, 0.5f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{1, 0, b, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.outputs(), 1);
    ASSERT_EQ(c.inputs(), 2);
    ASSERT_EQ(c.states(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, pidTF(a, s), s);
        assertChannel(c, 0, 1, pidTF(b, s), s);
    }
}

void test_assembler_duplicates_sum() {
    // Unreachable through the pairing grid, but correct for a loop list built
    // programmatically — and it pins the += on D.
    const LinearSystem plant = makePlant(1, 1);
    PIDParams a{1.5f, 2.0f, 0.25f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{0, 0, a, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 4);
    for (const cd& s : kS) assertChannel(c, 0, 0, 2.0 * pidTF(a, s), s);
}

void test_assembler_empty_and_out_of_range() {
    const LinearSystem plant = makePlant(2, 1);

    const LinearSystem empty = buildLoopController({}, LoopKind::PID, plant);
    ASSERT_EQ(empty.outputs(), 1);
    ASSERT_EQ(empty.inputs(), 2);
    ASSERT_EQ(empty.states(), 0);
    ASSERT_NEAR(empty.D.cwiseAbs().sum(), 0.0, 1e-15);

    PIDParams a{2.0f, 1.0f, 0.0f, 0.1f};
    PIDParams bad{9.0f, 9.0f, 0.0f, 0.1f};
    std::vector<Loop> loops = {Loop{0, 0, a, {}}, Loop{5, 0, bad, {}},
                               Loop{0, 3, bad, {}}};
    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 1);           // only the valid loop contributed
    for (const cd& s : kS) assertChannel(c, 0, 0, pidTF(a, s), s);
}

void test_assembler_leadlag_kind() {
    const LinearSystem plant = makePlant(2, 2);
    LeadLagParams l0{CompensatorMode::Lead, 2.0f, 0.1f, 1.0f, 10.0f, 1.0f};
    LeadLagParams l1{CompensatorMode::Lag, 1.0f, 0.1f, 1.0f, 20.0f, 0.5f};
    std::vector<Loop> loops = {Loop{0, 0, {}, l0}, Loop{1, 1, {}, l1}};

    const LinearSystem c = buildLoopController(loops, LoopKind::LeadLag, plant);
    ASSERT_EQ(c.states(), 2);
    for (const cd& s : kS) {
        assertChannel(c, 0, 0, leadLagTF(l0, s), s);
        assertChannel(c, 1, 1, leadLagTF(l1, s), s);
    }
}

// ------------------------------------------------------------ integration

void test_series_connect_is_square() {
    // The map's "closed loop stays p x p" claim, end to end, on the SIMO shape
    // where it is least obvious.  Also exercises the n = 0 path: the P-only
    // loop contributes no states.
    const LinearSystem plant = makePlant(2, 1);
    std::vector<Loop> loops = {Loop{0, 0, PIDParams{2.0f, 0.0f, 0.0f, 0.1f}, {}}};

    const LinearSystem c = buildLoopController(loops, LoopKind::PID, plant);
    ASSERT_EQ(c.states(), 0);

    const LinearSystem open = seriesConnect(c, plant);
    ASSERT_EQ(open.outputs(), 2);
    ASSERT_EQ(open.inputs(), 2);
    ASSERT_EQ(open.outputs(), open.inputs());

    const LinearSystem closed = feedbackConnect(open);
    ASSERT_EQ(closed.outputs(), 2);
    ASSERT_EQ(closed.inputs(), 2);
}

}  // anonymous namespace

int main() {
    test_pid_transfer_functions();
    test_pid_minimal_realization();
    test_leadlag_transfer_functions();
    test_leadlag_degenerate_alpha();
    test_assembler_shape_and_placement();
    test_assembler_unpaired_channel_is_inert();
    test_assembler_miso_fan_in_sums();
    test_assembler_duplicates_sum();
    test_assembler_empty_and_out_of_range();
    test_assembler_leadlag_kind();
    test_series_connect_is_square();
    std::printf("All controller_builders tests passed.\n");
    return 0;
}
```

- [ ] **Step 4: Wire into CMake**

Add `src/analysis/controller_builders.cpp` to the `analysis_lib` source list (after `src/analysis/system_connect.cpp`), and add a test target inside the `if(NOT EMSCRIPTEN)` block, cloned from the existing six:

```cmake
    add_executable(test_controller_builders tests/test_controller_builders.cpp)
    target_link_libraries(test_controller_builders PRIVATE analysis_lib)
    target_include_directories(test_controller_builders PRIVATE tests)
    add_test(NAME test_controller_builders COMMAND test_controller_builders)
```

- [ ] **Step 5: Build and test, in both configs**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug
cmake --build build-debug -j8 && ctest --test-dir build-debug --output-on-failure
```

The debug run is the one that proves Task 1's guard 3: without it, `test_series_connect_is_square` and `test_pid_minimal_realization` hit the `ColPivHouseholderQR` assert on the P-only case. Delete `build-debug/` afterwards.

- [ ] **Step 6: Commit**

```bash
git add src/analysis/controller_builders.h src/analysis/controller_builders.cpp \
        tests/test_controller_builders.cpp CMakeLists.txt
git commit -m "feat: PID and Lead/Lag builders over a loop-pairing list

buildPID/buildLeadLag emit minimal realizations — a state only where its
output coefficient is nonzero, so a zero gain cannot plant a phantom pole in
the closed-loop plot. buildLoopController scatters each loop's SISO block
into (in_k, out_k) of an m-out/p-in controller, so square, SIMO and MISO
plants share one data model.

Refs #2, #3"
```

---

### Task 3: AppState — loop list, scales, rule-C indices, enums

**Files:**
- Modify: `src/app_state.h`

Pure data. The app still builds and behaves exactly as before after this task; nothing reads the new fields yet.

- [ ] **Step 1: Include the builders header and `<algorithm>`**

In `src/app_state.h`, add to the include block (line 4-10):

```cpp
#include "analysis/controller_builders.h"
```

and `#include <algorithm>` beside `<cstring>`.

- [ ] **Step 2: Extend `ControllerType` and replace `PZMode`**

Replace line 36-38:

```cpp
enum class ControllerType { None, PID, LeadLag, StateSpace, GainMatrix };
enum class InputType { Step, Impulse, Ramp };

// Flat 4-way; two entries are conditionally valid.  PlantLocus is the old
// UnityFB, unchanged and honestly relabelled — it sweeps a proportional gain on
// the bare plant channel and does not put the controller in the path.
// LoopLocus sweeps kappa on the whole compensator.  See issue #10.
enum class PZMode { PoleZero, PlantLocus, LoopLocus, StateFB };
```

**Both enums are reordered, and every `static_cast<int>` / `static_cast<Enum>` round-trip through a UI index must be re-checked** — `model_panel.cpp:284-287` and `pole_zero_panel.cpp:29-36` both cast raw combo/radio indices. Tasks 7 and 8 rewrite those call sites.

- [ ] **Step 3: Add the loop list and its dimension cache**

In `AppState`, replace the `// --- Model ---` block's tail and the `// --- Channel selection ---` block:

```cpp
    // --- Controller: loop list ---
    // Maintained for the current plant regardless of ctrl_type, so
    // None -> PID -> None -> PID round-trips without losing tuning and picking
    // PID needs no seeding path of its own.  See issue #2.
    std::vector<Loop> loops;
    int cached_p = -1;  // plant dimensions the loop list was seeded for
    int cached_m = -1;
    int selected_loop = -1;  // last cell clicked in the pairing grid; -1 = none

    // Per-output scale: the maximum allowed deviation for that output.  Channel
    // Share is not scale-invariant and S&P's scaling is engineering intent, not
    // a model property, so it comes from the user.  Reseeded to 1.0 by the same
    // (p, m) guard as `loops`.  See issue #8.
    std::vector<float> output_scales;

    // --- Channel selection (rule C, issue #9) ---
    int output_i = 0;  // [0, p) — plant output, controller input, loop output
    int input_j = 0;   // [0, m) — plant input, controller output
    int ref_r = 0;     // [0, p) — reference channel of the open/closed loop
    bool show_all_channels = false;

    // Set by the recompute path: false when the selected channel has no loop on
    // it.  Dead channels are suppressed with a stated reason rather than drawn
    // with a fabricated flat 0 phase.
    bool channel_live[NUM_SYSTEMS] = {true, true, true, true};
    std::string channel_reason[NUM_SYSTEMS];
```

- [ ] **Step 4: Add the diagnostic controls**

Beside the Bode settings:

```cpp
    // --- Pairing diagnostic (RGA / Channel Share) ---
    // One omega, plus a collapsed RGA-number-vs-omega sweep across the existing
    // Bode grid.  S&P rule 2 gets no separate steady-state frequency: it is a
    // DIC theorem requiring a stable plant, and neither Ball-Balancer (poles at
    // the origin) nor Inverted Pendulum qualifies.  See issues #4, #6.
    float diag_omega_hz = 1.0f;
    bool diag_sweep = false;
```

- [ ] **Step 5: Add the kappa triple and retire `time_input_j`**

In the pole-zero / root locus block, add:

```cpp
    // Loop Locus.  A third explicit triple beside rl_k_* and rl_alpha_*: the
    // three scales are genuinely different (kappa 0.01-100 log, K 0-100 linear,
    // alpha 0-2 linear), so a shared triple would re-default on nearly every
    // mode switch.  kappa_min > 0 is a correctness constraint, not a taste —
    // see invariant 8.
    float rl_kappa_min = 0.01f;
    float rl_kappa_max = 100.0f;
    float rl_kappa_current = 1.0f;
    double rl_loop_gain_margin = -1.0;  // kappa*, or < 0 for no crossing
```

**Delete `int time_input_j = 0;`** (line 96). The time panel uses rule C's own indices: `input_j` for the plant, `ref_r` for the closed loop. Today's field is bounded by `m` but fed to both a p×m and a p×p system — on MISO it writes past the end of a length-p vector, and on SIMO it caps below `p`, so the unpaired reference can never be excited, which is exactly the "is this loop actually closed?" test.

- [ ] **Step 6: Add the seeding helper**

At the bottom of the file, beside `matrixToTextBuf`:

```cpp
// Seed the identity pairing y[k] -> u[k] for k < min(p, m), and reset the
// per-output scales.  Applied unconditionally on a dimension change.
//
// On Ball-Balancer this seeds a structurally DEAD pairing, deliberately: B(2,1)
// and B(3,0) mean input 0 drives output 1's chain, so G is exactly
// anti-diagonal and identity pairing leaves both loops open at every frequency.
// That is the pairing diagnostic's demonstration case — RGA of an anti-diagonal
// G is [[0,1],[1,0]], the loudest possible signal — and smart seeding would
// delete the lesson.  See issue #2.
inline void seedIdentityLoops(AppState& state) {
    state.loops.clear();
    const int p = state.plant.outputs();
    const int m = state.plant.inputs();
    for (int k = 0; k < std::min(p, m); ++k) {
        state.loops.push_back(Loop{k, k, {}, {}});
    }
    state.output_scales.assign(static_cast<std::size_t>(p), 1.0f);
    state.selected_loop = state.loops.empty() ? -1 : 0;
}
```

- [ ] **Step 7: Build**

```bash
cmake --build build -j8
```

It will fail on `time_input_j` in `src/visualizer.cpp` and `src/panels/time_response_panel.cpp`, and on the `PZMode::UnityFB` name in `src/panels/pole_zero_panel.cpp`. Those are Tasks 4, 5 and 7 — do not patch them here beyond what is needed to keep the tree compiling if you prefer to commit per task; otherwise carry the break into Task 4 and commit the two together.

---

### Task 4: Recompute path — dimension guard, seeding, controller assembly

**Files:**
- Modify: `src/visualizer.cpp`

- [ ] **Step 1: Add the dimension guard at the top of the recompute block**

In `render_frame`, immediately after `state.needs_recompute = false;` (line 62), insert:

```cpp
        // Plant dimensions changed: reseed the loop list wholesale.
        //
        // Reseed, not clamp and not partial-drop.  Clamping is the dangerous
        // option — a loop clamped from y1 to y0 controls a different physical
        // quantity while still carrying gains tuned for the old one, with
        // nothing on screen saying so.  Gains are lost, but a Kp tuned for a
        // quarter-car is meaningless on a pendulum.
        //
        // Keyed on dimensions, not plant equality, so the physical-parameter
        // sliders never disturb tuning mid-sweep.  Located here rather than in
        // the preset handler because FOUR sites mutate state.plant — the preset
        // combo (model_panel.cpp:110), Apply Plant (:181), the physical-param
        // sliders (:152) and applyTFToSS — and only a central guard covers all
        // of them, present and future.  See issue #2.
        if (state.plant.outputs() != state.cached_p ||
            state.plant.inputs() != state.cached_m) {
            state.cached_p = state.plant.outputs();
            state.cached_m = state.plant.inputs();
            caliburn::seedIdentityLoops(state);
            state.output_i = 0;
            state.input_j = 0;
            state.ref_r = 0;
        }
```

- [ ] **Step 2: Replace the controller composition block**

Replace lines 64-96 (from `state.systems[0] = state.plant;` through the closing `}` of the final `else`) with:

```cpp
        state.systems[0] = state.plant;
        state.system_valid[0] = true;

        const bool loop_backed =
            state.ctrl_type == caliburn::ControllerType::PID ||
            state.ctrl_type == caliburn::ControllerType::LeadLag;
        const caliburn::LoopKind loop_kind =
            (state.ctrl_type == caliburn::ControllerType::LeadLag)
                ? caliburn::LoopKind::LeadLag
                : caliburn::LoopKind::PID;

        if (loop_backed && !state.loops.empty()) {
            // The controller is m-out/p-in by construction, so seriesConnect's
            // precondition holds and the open loop is p x p — feedbackConnect
            // needs no shape test.  systems[1] is written directly rather than
            // through state.ctrl_ss, so switching to PID does not clobber a
            // hand-entered State-Space controller.
            state.systems[1] = caliburn::buildLoopController(
                state.loops, loop_kind, state.plant);
            state.system_valid[1] = true;
            state.systems[2] =
                caliburn::seriesConnect(state.systems[1], state.plant);
            state.system_valid[2] = true;
            state.systems[3] = caliburn::feedbackConnect(state.systems[2]);
            state.system_valid[3] = true;
        } else if (state.ctrl_type == caliburn::ControllerType::StateSpace &&
                   state.ctrl_ss.D.size() > 0) {
            // Was `ctrl_ss.A.size() > 0`, which conflated "is a controller
            // configured" with "is it safe to analyze".  Task 1's guards make
            // the second question moot, and a 0-state controller is a
            // legitimate controller — D is the matrix that is always present
            // once Apply Controller has run.  See issue #5.
            state.systems[1] = state.ctrl_ss;
            state.system_valid[1] = true;
            if (state.ctrl_ss.outputs() == state.plant.inputs()) {
                state.systems[2] =
                    caliburn::seriesConnect(state.ctrl_ss, state.plant);
                state.system_valid[2] = true;
                if (state.systems[2].outputs() == state.systems[2].inputs()) {
                    state.systems[3] =
                        caliburn::feedbackConnect(state.systems[2]);
                    state.system_valid[3] = true;
                } else {
                    state.system_valid[3] = false;
                }
            } else {
                state.system_valid[2] = false;
                state.system_valid[3] = false;
            }
        } else if (state.ctrl_type == caliburn::ControllerType::GainMatrix &&
                   state.ctrl_K.size() > 0) {
            state.system_valid[1] = false;
            state.system_valid[2] = false;
            state.systems[3] =
                caliburn::stateFeedbackClose(state.plant, state.ctrl_K);
            state.system_valid[3] = true;
        } else {
            // Includes a loop-backed type with an empty loop list, which is
            // equivalent to ControllerType::None (issue #2).  An all-zero m x p
            // controller connects cleanly but gives an open loop of -inf dB at
            // every frequency, which is not worth handing to Bode and Nyquist.
            state.system_valid[1] = false;
            state.system_valid[2] = false;
            state.system_valid[3] = false;
        }
```

`loop_backed` and `loop_kind` stay in scope for Steps in Tasks 5 and 7.

- [ ] **Step 3: Build and run**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
```

`time_input_j` still breaks `visualizer.cpp:126-141` and `time_response_panel.cpp:60`, and `PZMode::UnityFB` breaks `visualizer.cpp:150` and the pole-zero panel. Task 5 and Task 7 fix those; if you want a green build here, temporarily point the time loop at `state.input_j` and the locus at `PZMode::PlantLocus`, both of which are the final answers anyway.

---

### Task 5: Rule C — channel indexing across the panels

**Files:**
- Modify: `src/visualizer.cpp`
- Modify: `src/panels/bode_panel.cpp`
- Modify: `src/panels/time_response_panel.cpp`
- Modify: `src/panels/model_panel.cpp`

The reported crash — plant indices applied to `systems[1]` — was one symptom of a wider defect: the **same index pair applied to four systems of three shapes**. The closed loop carries the identical defect (on a MISO plant `input_j` reaches 1 while `systems[3].D` is 1×1) and so does the time panel.

- [ ] **Step 1: Add the channel map and liveness predicates**

In `render_frame`, after the controller composition block:

```cpp
        // Rule C.  output_i in [0,p), input_j in [0,m), ref_r in [0,p).
        //   plant           p x m -> (i, j)
        //   controller      m x p -> (j, i)   transposed
        //   open/closed     p x p -> (i, r)
        // The output index is genuinely shared: plant output i, controller
        // input i and loop-system output i are the same physical signal — so
        // this is three indices, not six, and every system is in range by
        // construction.  See issue #9.
        auto channelFor = [&](int s, int& out, int& in) {
            switch (s) {
                case 1:  out = state.input_j;  in = state.output_i; break;
                case 2:
                case 3:  out = state.output_i; in = state.ref_r;    break;
                default: out = state.output_i; in = state.input_j;  break;
            }
        };

        // Liveness is read off the loop list, never by inspecting matrices.
        auto loopOn = [&](int out, int in) {
            for (const auto& l : state.loops)
                if (l.out == out && l.in == in) return true;
            return false;
        };
        auto refLive = [&](int r) {
            for (const auto& l : state.loops)
                if (l.out == r) return true;
            return false;
        };

        for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
            state.channel_live[s] = true;
            state.channel_reason[s].clear();
        }
        if (loop_backed && !state.loops.empty()) {
            char buf[96];
            if (!loopOn(state.output_i, state.input_j)) {
                state.channel_live[1] = false;
                std::snprintf(buf, sizeof(buf), "no loop pairs y%d \xe2\x86\x92 u%d",
                              state.output_i, state.input_j);
                state.channel_reason[1] = buf;
            }
            if (!refLive(state.ref_r)) {
                state.channel_live[2] = state.channel_live[3] = false;
                std::snprintf(buf, sizeof(buf),
                              "reference r%d drives no loop", state.ref_r);
                state.channel_reason[2] = state.channel_reason[3] = buf;
            }
        }
```

**Why suppress rather than plot.** The current behaviour is half-defined, which is worse than either: `frequency_response.cpp:47` clamps `mag == 0` to −300 dB so the floor trace is free, but phase is `arg(0) == 0`, drawing a **fabricated flat 0° line**; and `computePoleZero` takes the pencil branch with `B.col(j)` and `C.row(i)` both zero, so the finite generalized eigenvalues come out as `eig(A)` again — every pole reported with a zero exactly on it, which reads as total pole-zero cancellation. Under rule C, the primary way a wrong pairing is discovered is that you select the channel you believe is closed and the compensator trace is not there; a plausible flat line hides exactly that.

- [ ] **Step 2: Rewrite the single-channel Bode/PZ loop**

Replace lines 98-106 (the `for (int s = 0; ...)` computing `bode[s]` and `pole_zero[s]`):

```cpp
        for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
            if (!state.system_valid[s] || s == 2) {
                state.bode[s] = {};
                state.pole_zero[s] = {};
                continue;
            }
            if (!state.channel_live[s]) {
                state.bode[s] = {};
                state.pole_zero[s] = {};
                continue;
            }
            int out = 0, in = 0;
            channelFor(s, out, in);
            state.bode[s] = caliburn::computeBode(
                state.systems[s], out, in,
                state.freq_min_hz, state.freq_max_hz, state.num_freq_points);
            state.pole_zero[s] = caliburn::computePoleZero(
                state.systems[s], out, in);
        }
```

**The controller's pole plot stays unfiltered.** `pole_zero.cpp:38-47` takes poles as `eig(A)` for the whole system regardless of channel; only zeros are channel-specific. That is the correct convention and what the plant already does. It reads oddly here — the controller's `A` is block-diagonal over the loops, so selecting loop 0's channel still shows loop 1's integrator pole beside loop 0's zero — but filtering would require `buildLoopController` to report each loop's state range, a **new return contract on a closed decision**, and would make the controller and the plant behave differently for the same plot. Leave it.

- [ ] **Step 3: Split the channel grid by role**

Replace the `if (state.show_all_channels)` block (lines 108-124). The writer was right and the reader was wrong: `visualizer.cpp` sizes each `bode_grid[s]` from **that system's own** `p, m`, while `bode_panel.cpp:281-283` lays the table out from **`state.plant`**'s and indexes `idx = i*m + j` with the plant's `m`. Overlaying four systems in one cell only ever worked when all four were p×m — true exactly when the plant is square.

```cpp
        if (state.show_all_channels) {
            const int p = state.plant.outputs();
            const int m = state.plant.inputs();
            // Two grids, by role.  Rule C applied to the grid, so the grid needs
            // no rule of its own and each table is indexed by the shape it
            // actually holds.
            //   open grid  p x m : plant at (i,j), controller at (j,i)
            //   loop grid  p x p : open/closed loop at (i,r)
            for (int s = 0; s < caliburn::NUM_SYSTEMS; ++s) {
                if (!state.system_valid[s] || s == 2) {
                    state.bode_grid[s].clear();
                    continue;
                }
                const bool is_loop_system = (s == 3);
                const int rows = p;
                const int cols = is_loop_system ? p : m;
                state.bode_grid[s].assign(static_cast<std::size_t>(rows * cols),
                                          {});
                for (int i = 0; i < rows; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        // Cell (i,j) of the open grid means plant channel (i,j)
                        // and controller channel (j,i) — the same physical
                        // path, read from each end.
                        int out = i, in = j;
                        if (s == 1) { out = j; in = i; }
                        state.bode_grid[s][i * cols + j] = caliburn::computeBode(
                            state.systems[s], out, in,
                            state.freq_min_hz, state.freq_max_hz,
                            state.num_freq_points);
                    }
                }
            }
        }
```

`systems[2]` stays excluded from Bode and pole-zero, as it is today — so the loop grid holds the closed loop alone. The controller's cells are mostly blank by construction (one nonzero per loop); that is informative, not broken. The open grid becomes a picture of which channels are actually closed.

- [ ] **Step 4: Rewrite the time-response loop**

Replace `state.time_input_j` in lines 126-141 with the per-system index:

```cpp
        for (int s : {0, 3}) {
            if (!state.system_valid[s]) continue;
            // Plant is p x m, closed loop is p x p — the same shape mismatch,
            // in the time domain.  No parallel time-domain index pair: nothing
            // in the panel ever suggested decoupling it from the Bode channel
            // was deliberate.  See issue #9.
            const int u_idx = (s == 0) ? state.input_j : state.ref_r;
            switch (state.input_type) {
                case caliburn::InputType::Step:
                    state.time_resp[s] = caliburn::computeStepResponse(
                        state.systems[s], u_idx,
                        state.amplitude, state.duration, state.dt_sim);
                    break;
                case caliburn::InputType::Impulse:
                    state.time_resp[s] = caliburn::computeImpulseResponse(
                        state.systems[s], u_idx,
                        state.amplitude, state.duration, state.dt_sim);
                    break;
                case caliburn::InputType::Ramp:
                    state.time_resp[s] = caliburn::computeRampResponse(
                        state.systems[s], u_idx,
                        state.slope, state.duration, state.dt_sim);
                    break;
            }
        }
```

- [ ] **Step 5: Retire the time panel's own channel slider**

In `src/panels/time_response_panel.cpp`, delete the "Input ch" table column (lines 57-62) and reduce the surrounding `BeginTable` column count by one. Add a one-line readout in its place so the panel still says what it is exciting:

```cpp
        ImGui::TableNextColumn();
        ImGui::TextDisabled("plant u%d  /  closed-loop r%d",
                            state.input_j, state.ref_r);
```

- [ ] **Step 6: Add the `ref_r` slider to the channel selector**

In `src/panels/model_panel.cpp`, in the `// --- Channel selector ---` block (lines 360-374), add after the `Input j` slider:

```cpp
    if (ImGui::SliderInt("Reference r", &state.ref_r, 0, std::max(p - 1, 0))) {
        state.needs_recompute = true;
    }
```

Keep the `p <= 3 && m <= 3` gate on Show All Channels as it is.

- [ ] **Step 7: Rebuild the Bode grid panel around the two tables**

In `src/panels/bode_panel.cpp`, replace the `if (state.show_all_channels)` block (lines 281-316). The plant-dimensioned read at 282-283 and the `idx = i * m + j` at 292 are the bug; each table must be sized by the shape it holds.

```cpp
    if (state.show_all_channels) {
        const int p = state.plant.outputs();
        const int m = state.plant.inputs();

        auto drawGrid = [&](const char* id, const char* title, int rows,
                            int cols, const int* systems, int n_systems) {
            bool any = false;
            for (int k = 0; k < n_systems; ++k) {
                const int s = systems[k];
                if (state.trace_visible[s] && state.system_valid[s] &&
                    !state.bode_grid[s].empty()) any = true;
            }
            if (!any) return;
            ImGui::SeparatorText(title);
            if (!ImGui::BeginTable(id, cols, ImGuiTableFlags_Borders)) return;
            for (int i = 0; i < rows; ++i) {
                ImGui::TableNextRow();
                for (int j = 0; j < cols; ++j) {
                    ImGui::TableNextColumn();
                    ImGui::Text("(%d,%d)", i, j);
                    for (int k = 0; k < n_systems; ++k) {
                        const int s = systems[k];
                        if (!state.trace_visible[s] || !state.system_valid[s])
                            continue;
                        const int idx = i * cols + j;
                        if (idx >= (int)state.bode_grid[s].size()) continue;
                        const auto& bode = state.bode_grid[s][idx];
                        if (bode.points.empty()) continue;
                        const int nn = (int)bode.points.size();
                        std::vector<double> ff(nn), mm(nn);
                        for (int q = 0; q < nn; ++q) {
                            ff[q] = bode.points[q].freq_hz;
                            mm[q] = bode.points[q].magnitude_db;
                        }
                        char pid[64];
                        std::snprintf(pid, sizeof(pid), "##bode_%s_%d_%d_%d",
                                      id, i, j, s);
                        if (ImPlot::BeginPlot(pid, ImVec2(-1, 120))) {
                            ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                            ImPlot::PlotLine(system_names[s], ff.data(),
                                             mm.data(), nn,
                                             ImPlotSpec(ImPlotProp_LineColor,
                                                        system_colors[s]));
                            ImPlot::EndPlot();
                        }
                    }
                }
            }
            ImGui::EndTable();
        };

        // Open grid: plant (i,j) and controller (j,i) — the same channel read
        // from each end.  Loop grid: the p x p loop systems at (i,r).
        const int open_systems[2] = {0, 1};
        const int loop_systems[1] = {3};
        drawGrid("open", "Open channels  (y_i \xe2\x86\x90 u_j)", p, m,
                 open_systems, 2);
        drawGrid("loop", "Loop channels  (y_i \xe2\x86\x90 r_j)", p, p,
                 loop_systems, 1);
    } else {
```

- [ ] **Step 8: Render the stated reasons**

Still in `bode_panel.cpp`, just before `ImGui::End()`, add:

```cpp
    for (int s = 0; s < NUM_SYSTEMS; ++s) {
        if (s == 2 || state.channel_reason[s].empty()) continue;
        ImGui::TextDisabled("%s: %s", system_names[s],
                            state.channel_reason[s].c_str());
    }
```

Add the same three lines to `src/panels/pole_zero_panel.cpp` before its `ImGui::End()`. A suppressed trace must be distinguishable from breakage — the stated reason is what turns a dead channel into the diagnostic.

- [ ] **Step 9: Build, test, commit**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
git add src/app_state.h src/visualizer.cpp src/panels/bode_panel.cpp \
        src/panels/time_response_panel.cpp src/panels/model_panel.cpp
git commit -m "feat: loop list in AppState, and rule-C channel indexing

output_i/input_j were plant-shaped indices applied to four systems of three
shapes, so the controller read out of bounds and the closed loop and time
panel carried the same latent bug. Rule C adds ref_r and maps each system to
the channel it actually has; the all-channels view splits into two role-grids,
each sized by the shape it holds. Dead channels are suppressed with a stated
reason instead of a fabricated flat 0 phase.

Refs #2, #9"
```

---

### Task 6: Pairing diagnostics — RGA and Channel Share

**Files:**
- Create: `src/analysis/loop_diagnostics.h`
- Create: `src/analysis/loop_diagnostics.cpp`
- Create: `tests/test_loop_diagnostics.cpp`
- Modify: `CMakeLists.txt`

**[plan-level]** No ticket said where this computation lives; the prototype had it inside the panel. It goes in `analysis_lib` for the same reason the builders did — it is analysis, it needs no ImGui, and the `.transpose()`-not-`.adjoint()` trap is exactly the kind of bug only a test catches.

- [ ] **Step 1: Create `src/analysis/loop_diagnostics.h`**

```cpp
// src/analysis/loop_diagnostics.h
#pragma once

#include "../linear_system.h"
#include "controller_builders.h"

#include <Eigen/Core>
#include <string>
#include <vector>

namespace caliburn {

// The two readouts are mutually exclusive: a square paired sub-matrix of size
// >= 2 gives the RGA, everything else gives Channel Share.  One panel slot.
struct LoopDiagnostics {
    // --- RGA (pairing diagnostic) ---
    bool has_rga = false;
    Eigen::MatrixXd lambda;            // q x q, the real part of Λ
    std::vector<int> sub_out, sub_in;  // sub-matrix index maps into (p, m)
    double rga_number = 0.0;           // ‖Λ − Π‖_sum against the actual pairing

    // --- Channel Share (NOT a pairing measure) ---
    bool has_share = false;
    int share_in = 0;                  // the paired input the shares are of
    std::vector<double> share;         // per plant output, sums to 1
    std::vector<double> mag_db;        // |g_i(jω)| in dB, UNSCALED, per output

    // --- per loop ---
    std::vector<double> loop_lambda;   // NaN where there is no RGA
    std::vector<char> loop_dead;       // structurally dead paired channel

    bool scales_at_default = true;
    std::string headline;
};

// G(jω) as a complex matrix, one evalTransferFunction per channel.
Eigen::MatrixXcd evalPlantMatrix(const LinearSystem& sys, double freq_hz);

// |g_ij(jω)| below tol at EVERY ω on the Bode grid.  Positive rescaling cannot
// turn a zero into a nonzero at any frequency, so this test is immune to the
// scale-dependence of Channel Share — unlike any threshold on the share itself,
// or a ratio against the column maximum.
bool isStructurallyDeadChannel(const LinearSystem& sys, int output_i, int input_j,
                               double freq_min_hz, double freq_max_hz,
                               int num_points, double tol = 1e-9);

LoopDiagnostics computeLoopDiagnostics(const LinearSystem& plant,
                                       const std::vector<Loop>& loops,
                                       const std::vector<float>& output_scales,
                                       double freq_hz,
                                       double freq_min_hz, double freq_max_hz,
                                       int num_points);

// The RGA number alone at one ω, for the sweep strip.  NaN where no RGA exists.
double rgaNumberAt(const LinearSystem& plant, const std::vector<Loop>& loops,
                   double freq_hz);

}  // namespace caliburn
```

- [ ] **Step 2: Create `src/analysis/loop_diagnostics.cpp`**

```cpp
// src/analysis/loop_diagnostics.cpp
#include "loop_diagnostics.h"
#include "frequency_response.h"

#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace caliburn {

namespace {
constexpr double kPi = 3.14159265358979323846;

void pairedSets(const std::vector<Loop>& loops,
                std::vector<int>& outs, std::vector<int>& ins) {
    for (const auto& l : loops) {
        if (std::find(outs.begin(), outs.end(), l.out) == outs.end())
            outs.push_back(l.out);
        if (std::find(ins.begin(), ins.end(), l.in) == ins.end())
            ins.push_back(l.in);
    }
    std::sort(outs.begin(), outs.end());
    std::sort(ins.begin(), ins.end());
}
}  // anonymous namespace

Eigen::MatrixXcd evalPlantMatrix(const LinearSystem& sys, double freq_hz) {
    const std::complex<double> s(0.0, 2.0 * kPi * freq_hz);
    Eigen::MatrixXcd G(sys.outputs(), sys.inputs());
    for (int i = 0; i < sys.outputs(); ++i)
        for (int j = 0; j < sys.inputs(); ++j)
            G(i, j) = evalTransferFunction(sys, i, j, s);
    return G;
}

bool isStructurallyDeadChannel(const LinearSystem& sys, int output_i, int input_j,
                               double freq_min_hz, double freq_max_hz,
                               int num_points, double tol) {
    const double log_min = std::log10(freq_min_hz);
    const double log_max = std::log10(freq_max_hz);
    const int n = std::max(num_points, 2);
    for (int k = 0; k < n; ++k) {
        const double f = std::pow(10.0, log_min + (log_max - log_min) * k / (n - 1));
        const std::complex<double> s(0.0, 2.0 * kPi * f);
        if (std::abs(evalTransferFunction(sys, output_i, input_j, s)) > tol)
            return false;
    }
    return true;
}

LoopDiagnostics computeLoopDiagnostics(const LinearSystem& plant,
                                       const std::vector<Loop>& loops,
                                       const std::vector<float>& output_scales,
                                       double freq_hz,
                                       double freq_min_hz, double freq_max_hz,
                                       int num_points) {
    LoopDiagnostics d;
    const int p = plant.outputs();

    for (float v : output_scales)
        if (std::abs(v - 1.0f) > 1e-6f) { d.scales_at_default = false; break; }

    if (loops.empty()) {
        // Nothing paired: there is no input to take a column from.
        d.headline = "No loops \xe2\x80\x94 controller is None";
        return d;
    }

    const Eigen::MatrixXcd G = evalPlantMatrix(plant, freq_hz);

    d.loop_lambda.assign(loops.size(), std::nan(""));
    d.loop_dead.assign(loops.size(), 0);
    for (std::size_t k = 0; k < loops.size(); ++k) {
        d.loop_dead[k] = isStructurallyDeadChannel(
            plant, loops[k].out, loops[k].in,
            freq_min_hz, freq_max_hz, num_points) ? 1 : 0;
    }

    std::vector<int> outs, ins;
    pairedSets(loops, outs, ins);

    // The RGA exists only on a SQUARE paired sub-matrix of size >= 2.  A 1x1
    // sub-matrix yields Λ = [1] for every plant at every ω, so a green 1.00
    // would be actively misleading.  Note this is a property of the LOOP LIST,
    // not the plant shape: it also excludes a multi-input plant with one loop,
    // and the MISO fan-out whose sub-matrix is 1x2 — the case that passes a
    // naive `loops.size() >= 2` check.  See issues #4, #8.
    if (outs.size() == ins.size() && outs.size() >= 2) {
        const int q = static_cast<int>(outs.size());
        Eigen::MatrixXcd sub(q, q);
        for (int r = 0; r < q; ++r)
            for (int c = 0; c < q; ++c) sub(r, c) = G(outs[r], ins[c]);

        // Λ = G ∘ (G⁻¹)ᵀ on the COMPLEX matrix.  .transpose(), NEVER
        // .adjoint() — the two give wildly different answers and no
        // real-matrix test can tell them apart (S&P footnote 5, p. 83).
        const Eigen::MatrixXcd L =
            sub.array() * sub.inverse().transpose().array();

        d.has_rga = true;
        d.lambda = L.real();
        d.sub_out = outs;
        d.sub_in = ins;

        Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(q, q);
        for (const auto& l : loops) {
            const int r = (int)(std::find(outs.begin(), outs.end(), l.out) - outs.begin());
            const int c = (int)(std::find(ins.begin(), ins.end(), l.in) - ins.begin());
            Pi(r, c) = 1.0;
        }
        // Severity comes from the RGA number, NEVER from |λ_ii|: for
        // G = [[1,2],[1,1]], |λ_11| = 1.0000 looks perfect while λ_11 = −1 and
        // the RGA number is 8.
        d.rga_number = (d.lambda - Pi).cwiseAbs().sum();

        for (std::size_t k = 0; k < loops.size(); ++k) {
            const int r = (int)(std::find(outs.begin(), outs.end(), loops[k].out) - outs.begin());
            const int c = (int)(std::find(ins.begin(), ins.end(), loops[k].in) - ins.begin());
            d.loop_lambda[k] = d.lambda(r, c);
        }

        char buf[192];
        std::snprintf(buf, sizeof(buf),
                      "RGA @ %.3g Hz \xe2\x80\x94 RGA number %.2f",
                      freq_hz, d.rga_number);
        d.headline = buf;
        return d;
    }

    // Channel Share.  Not a pairing measure; never coloured or thresholded as
    // one.  λ_i = |g_i|²/Σ|g_k|² against the user's per-output scales.
    d.has_share = true;
    d.share_in = loops.front().in;
    d.share.assign(p, 0.0);
    d.mag_db.assign(p, 0.0);
    double total = 0.0;
    for (int i = 0; i < p; ++i) {
        const double scale = (i < (int)output_scales.size() &&
                              output_scales[i] > 1e-12f)
                                 ? output_scales[i] : 1.0;
        const double gi = std::abs(G(i, d.share_in)) / scale;
        d.share[i] = gi * gi;
        total += d.share[i];
        // The dB column is load-bearing, not decoration: the share is
        // normalised, so it cannot distinguish "drives both outputs strongly"
        // from "drives neither" — both read 50/50.  S&P's C3 has two halves,
        // and for a p x 1 matrix the column sum is exactly 1 by construction,
        // so the input-deletion half is vacuous.  Only the output half
        // survives, and it needs the magnitudes.
        d.mag_db[i] = 20.0 * std::log10(std::max(std::abs(G(i, d.share_in)), 1e-300));
    }
    for (int i = 0; i < p; ++i) d.share[i] = (total > 0.0) ? d.share[i] / total : 0.0;

    char buf[192];
    const char* why =
        (plant.inputs() == 1)              ? "single-input plant"
      : (outs.size() != ins.size())        ? "pairing is non-square"
                                           : "one loop of several inputs \xe2\x80\x94 RGA needs two";
    std::snprintf(buf, sizeof(buf),
                  "Channel Share of u%d @ %.3g Hz \xe2\x80\x94 %s; scales %s",
                  d.share_in, freq_hz, why,
                  d.scales_at_default ? "at default (all 1)" : "user-set");
    d.headline = buf;
    return d;
}

double rgaNumberAt(const LinearSystem& plant, const std::vector<Loop>& loops,
                   double freq_hz) {
    // Cheap path: the sweep strip only ever wants the number.
    std::vector<int> outs, ins;
    pairedSets(loops, outs, ins);
    if (outs.size() != ins.size() || outs.size() < 2) return std::nan("");

    const Eigen::MatrixXcd G = evalPlantMatrix(plant, freq_hz);
    const int q = static_cast<int>(outs.size());
    Eigen::MatrixXcd sub(q, q);
    for (int r = 0; r < q; ++r)
        for (int c = 0; c < q; ++c) sub(r, c) = G(outs[r], ins[c]);
    const Eigen::MatrixXd lambda =
        (sub.array() * sub.inverse().transpose().array()).real();

    Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(q, q);
    for (const auto& l : loops) {
        const int r = (int)(std::find(outs.begin(), outs.end(), l.out) - outs.begin());
        const int c = (int)(std::find(ins.begin(), ins.end(), l.in) - ins.begin());
        Pi(r, c) = 1.0;
    }
    return (lambda - Pi).cwiseAbs().sum();
}

}  // namespace caliburn
```

**Do not evaluate at ω = 0.** `colPivHouseholderQr().solve` on a singular `sI − A` fails *silently*, and pairing is done at the intended closed-loop bandwidth anyway — S&P give an example where `Λ(0) = Λ(j∞) = I` yet the correct pairing is the reverse. The ω slider's range (Task 8) starts at 0.01 Hz.

- [ ] **Step 3: Create `tests/test_loop_diagnostics.cpp`**

```cpp
// tests/test_loop_diagnostics.cpp
#include "analysis/loop_diagnostics.h"
#include "test_helpers.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace caliburn;

namespace {

// The Ball-Balancer shape: G = (k)/s^2 * [[0,1],[1,0]].  B(2,1) and B(3,0) mean
// input 0 drives output 1's chain, so the plant is exactly anti-diagonal.
LinearSystem makeAntiDiagonalIntegrator(double k) {
    LinearSystem g;
    g.A = Eigen::MatrixXd::Zero(4, 4);
    g.A(0, 2) = 1.0;
    g.A(1, 3) = 1.0;
    g.B = Eigen::MatrixXd::Zero(4, 2);
    g.B(2, 1) = k;
    g.B(3, 0) = k;
    g.C = Eigen::MatrixXd::Zero(2, 4);
    g.C(0, 0) = 1.0;
    g.C(1, 1) = 1.0;
    g.D = Eigen::MatrixXd::Zero(2, 2);
    return g;
}

// Coupled, genuinely complex on the jω axis, invertible.
LinearSystem makeCoupled() {
    LinearSystem g;
    g.A = (Eigen::MatrixXd(2, 2) << -1.0, 0.0, 0.0, -2.0).finished();
    g.B = (Eigen::MatrixXd(2, 2) << 1.0, 1.0, 1.0, 3.0).finished();
    g.C = Eigen::MatrixXd::Identity(2, 2);
    g.D = Eigen::MatrixXd::Zero(2, 2);
    return g;
}

LinearSystem makeSimo() {  // 1 input, 2 outputs
    LinearSystem g;
    g.A = (Eigen::MatrixXd(2, 2) << -1.0, 0.0, 0.0, -5.0).finished();
    g.B = (Eigen::MatrixXd(2, 1) << 1.0, 4.0).finished();
    g.C = Eigen::MatrixXd::Identity(2, 2);
    g.D = Eigen::MatrixXd::Zero(2, 1);
    return g;
}

const std::vector<float> kUnitScales2 = {1.0f, 1.0f};

// An integrator plant is fine, provably: by scaling invariance Λ is
// [[0,1],[1,0]] at EVERY ω > 0, with cond(G) = 1 throughout.  The infinite DC
// gain that made the spec's −C·A⁻¹·B heuristic undefined cancels exactly here.
void test_antidiagonal_rga_is_exact_at_every_omega() {
    const LinearSystem g = makeAntiDiagonalIntegrator(5.0 * 9.81 / 7.0);
    const std::vector<Loop> identity = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};

    for (double f : {1e-4, 1e-2, 1.0, 1e2, 1e4}) {
        const auto d = computeLoopDiagnostics(g, identity, kUnitScales2,
                                              f, 0.01, 100.0, 60);
        ASSERT_TRUE(d.has_rga);
        ASSERT_TRUE(!d.has_share);
        ASSERT_NEAR(d.lambda(0, 0), 0.0, 1e-9);
        ASSERT_NEAR(d.lambda(0, 1), 1.0, 1e-9);
        ASSERT_NEAR(d.lambda(1, 0), 1.0, 1e-9);
        ASSERT_NEAR(d.lambda(1, 1), 0.0, 1e-9);
        // 4.00 for the naive diagonal pairing — the loudest possible signal.
        ASSERT_NEAR(d.rga_number, 4.0, 1e-9);
    }

    const std::vector<Loop> swapped = {Loop{0, 1, {}, {}}, Loop{1, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, swapped, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(d.rga_number, 0.0, 1e-9);
}

// The transpose/adjoint trap.  Row and column sums of the true RGA are exactly
// 1 (S&P A.4.2); the adjoint form does not preserve that for complex G.  This
// is the only test that can fail if .adjoint() creeps in.
void test_rga_row_and_column_sums_are_one() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> loops = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    const auto d = computeLoopDiagnostics(g, loops, kUnitScales2,
                                          0.5, 0.01, 100.0, 60);
    ASSERT_TRUE(d.has_rga);
    for (int r = 0; r < 2; ++r) ASSERT_NEAR(d.lambda.row(r).sum(), 1.0, 1e-9);
    for (int c = 0; c < 2; ++c) ASSERT_NEAR(d.lambda.col(c).sum(), 1.0, 1e-9);
}

// A 1x1 sub-matrix gives Λ = [1] for every plant at every ω.  Trigger on the
// loop list, not on plant.inputs() == 1.
void test_single_loop_on_square_plant_falls_back_to_share() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, one, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(d.has_share);
}

// MISO fan-out: distinct outputs 2, distinct inputs 1 -> non-square 2x1.
// Two loops, so it passes a naive loops.size() >= 2 check.
void test_miso_fan_out_falls_back_to_share() {
    const LinearSystem g = makeSimo();
    const std::vector<Loop> two = {Loop{0, 0, {}, {}}, Loop{1, 0, {}, {}}};
    const auto d = computeLoopDiagnostics(g, two, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(d.has_share);
    ASSERT_EQ(d.share_in, 0);
}

void test_share_sums_to_one_and_respects_scales() {
    const LinearSystem g = makeSimo();
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};

    const auto flat = computeLoopDiagnostics(g, one, kUnitScales2,
                                             1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(flat.share[0] + flat.share[1], 1.0, 1e-12);
    ASSERT_TRUE(flat.scales_at_default);

    const std::vector<float> scaled = {1.0f, 10.0f};
    const auto s = computeLoopDiagnostics(g, one, scaled, 1.0, 0.01, 100.0, 60);
    ASSERT_NEAR(s.share[0] + s.share[1], 1.0, 1e-12);
    ASSERT_TRUE(!s.scales_at_default);
    // Scaling output 1 down by 10 must raise output 0's share.
    ASSERT_TRUE(s.share[0] > flat.share[0]);
    // The dB column is unscaled, so it does not move.
    ASSERT_NEAR(s.mag_db[1], flat.mag_db[1], 1e-12);
}

void test_structurally_dead_channel() {
    const LinearSystem g = makeAntiDiagonalIntegrator(1.0);
    // The paired identity channels are dead at every ω; the off-diagonal ones
    // are not.
    ASSERT_TRUE(isStructurallyDeadChannel(g, 0, 0, 0.01, 100.0, 60));
    ASSERT_TRUE(isStructurallyDeadChannel(g, 1, 1, 0.01, 100.0, 60));
    ASSERT_TRUE(!isStructurallyDeadChannel(g, 0, 1, 0.01, 100.0, 60));
    ASSERT_TRUE(!isStructurallyDeadChannel(g, 1, 0, 0.01, 100.0, 60));

    const std::vector<Loop> identity = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    const auto d = computeLoopDiagnostics(g, identity, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_EQ((int)d.loop_dead[0], 1);
    ASSERT_EQ((int)d.loop_dead[1], 1);
}

void test_empty_loop_list_shows_nothing() {
    const LinearSystem g = makeCoupled();
    const auto d = computeLoopDiagnostics(g, {}, kUnitScales2,
                                          1.0, 0.01, 100.0, 60);
    ASSERT_TRUE(!d.has_rga);
    ASSERT_TRUE(!d.has_share);
}

void test_rga_number_sweep_matches_full_computation() {
    const LinearSystem g = makeCoupled();
    const std::vector<Loop> loops = {Loop{0, 0, {}, {}}, Loop{1, 1, {}, {}}};
    for (double f : {0.05, 0.5, 5.0, 50.0}) {
        const auto d = computeLoopDiagnostics(g, loops, kUnitScales2,
                                              f, 0.01, 100.0, 60);
        ASSERT_NEAR(rgaNumberAt(g, loops, f), d.rga_number, 1e-12);
    }
    const std::vector<Loop> one = {Loop{0, 0, {}, {}}};
    ASSERT_TRUE(std::isnan(rgaNumberAt(g, one, 1.0)));
}

}  // anonymous namespace

int main() {
    test_antidiagonal_rga_is_exact_at_every_omega();
    test_rga_row_and_column_sums_are_one();
    test_single_loop_on_square_plant_falls_back_to_share();
    test_miso_fan_out_falls_back_to_share();
    test_share_sums_to_one_and_respects_scales();
    test_structurally_dead_channel();
    test_empty_loop_list_shows_nothing();
    test_rga_number_sweep_matches_full_computation();
    std::printf("All loop_diagnostics tests passed.\n");
    return 0;
}
```

- [ ] **Step 4: Wire into CMake, build, test, commit**

Add `src/analysis/loop_diagnostics.cpp` to `analysis_lib`, and a `test_loop_diagnostics` target cloned from the others.

```bash
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build --output-on-failure
git add src/analysis/loop_diagnostics.h src/analysis/loop_diagnostics.cpp \
        tests/test_loop_diagnostics.cpp CMakeLists.txt
git commit -m "feat: RGA and Channel Share pairing diagnostics

RGA on the complex G(jw) at a user-set frequency, defined only on a square
paired sub-matrix of size >= 2; Channel Share everywhere else, against
user-set per-output scales and never presented as a pairing measure. The
row-sum test is what catches .adjoint() creeping in where .transpose()
belongs — no real-matrix test can.

Refs #4, #8"
```

---

### Task 7: Loop Locus and the four-way PZMode

**Files:**
- Create: `src/analysis/loop_locus.h`
- Create: `src/analysis/loop_locus.cpp`
- Modify: `src/visualizer.cpp`
- Modify: `src/panels/pole_zero_panel.cpp`
- Modify: `CMakeLists.txt`

Today neither locus mode has anything to do with the dynamic controller: `UnityFB` sweeps a proportional gain on the **bare plant channel** and ignores `ctrl_ss` entirely. `LoopLocus` is the first locus that responds to Tier 1 tuning.

- [ ] **Step 1: Create `src/analysis/loop_locus.h`**

```cpp
// src/analysis/loop_locus.h
#pragma once

#include "../linear_system.h"
#include "controller_builders.h"
#include "pole_zero.h"

#include <vector>

namespace caliburn {

// Root locus of 1 + kappa*L(s) = 0, where kappa scales the WHOLE compensator on
// rule C's channel (input_j, output_i) — every loop with out == output_i and
// in == input_j, so duplicate pairings sum naturally and no loop selector is
// needed.  Other loops are held at their current tuning.
//
// This module sits ABOVE pole_zero.h so that module stays what it is today: a
// general routine over any LinearSystem with no knowledge of controllers.
//
// kappa is log-spaced over [kappa_min, kappa_max], kappa_min > 0.  Returns
// empty when the channel carries no loop.
std::vector<RootLocusPoint> computeLoopLocus(
    const LinearSystem& plant,
    const std::vector<Loop>& loops,
    LoopKind kind,
    int output_i, int input_j,
    double kappa_min, double kappa_max, int num_points);

// The selected loop's gain margin with every other loop held at its tuning:
// the first sign change in max Re(pole), linearly interpolated.  Returns a
// negative value when there is no crossing in range — never a fabricated number.
double loopGainMargin(const std::vector<RootLocusPoint>& locus);

}  // namespace caliburn
```

- [ ] **Step 2: Create `src/analysis/loop_locus.cpp`**

```cpp
// src/analysis/loop_locus.cpp
#include "loop_locus.h"
#include "system_connect.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <limits>

namespace caliburn {

std::vector<RootLocusPoint> computeLoopLocus(
    const LinearSystem& plant,
    const std::vector<Loop>& loops,
    LoopKind kind,
    int output_i, int input_j,
    double kappa_min, double kappa_max, int num_points) {
    std::vector<RootLocusPoint> result;
    if (loops.empty() || kappa_min <= 0.0 || kappa_max <= kappa_min ||
        num_points < 2)
        return result;

    // Same liveness predicate the Bode panel uses — one rule, two panels.
    // While inspecting an off-loop plant channel (the Ball-Balancer pairing
    // demo) the locus is correctly blank.
    bool live = false;
    for (const auto& l : loops)
        if (l.out == output_i && l.in == input_j) { live = true; break; }
    if (!live) return result;

    const double log_min = std::log10(kappa_min);
    const double log_max = std::log10(kappa_max);
    result.reserve(num_points);

    for (int step = 0; step < num_points; ++step) {
        const double kappa = std::pow(
            10.0, log_min + (log_max - log_min) * step / (num_points - 1));

        // kappa * C_k(s).  Scaling all three PID gains together is exactly
        // that, since kappa*(Kp + Ki/s + Kd*s/(tau_f*s+1)) leaves tau_f alone.
        // kappa != 0 never changes WHICH gains are zero, so the minimal
        // realization emits the identical state set at every point — the
        // closed-loop dimension is constant by construction, which is what
        // makes matchPoles and the panel's fixed-width drawing loop safe.
        std::vector<Loop> scaled = loops;
        for (auto& l : scaled) {
            if (l.out != output_i || l.in != input_j) continue;
            l.pid.Kp = static_cast<float>(l.pid.Kp * kappa);
            l.pid.Ki = static_cast<float>(l.pid.Ki * kappa);
            l.pid.Kd = static_cast<float>(l.pid.Kd * kappa);
            l.leadlag.Kc = static_cast<float>(l.leadlag.Kc * kappa);
        }

        // Literally the production path, so the locus cannot disagree with the
        // closed-loop pole plot drawn beside it.  buildLoopController's
        // contract is untouched.
        const LinearSystem ctrl = buildLoopController(scaled, kind, plant);
        const LinearSystem closed =
            feedbackConnect(seriesConnect(ctrl, plant));

        const int n = closed.states();
        std::vector<std::complex<double>> poles(n);
        if (n > 0) {
            Eigen::EigenSolver<Eigen::MatrixXd> es(closed.A, false);
            for (int i = 0; i < n; ++i) poles[i] = es.eigenvalues()(i);
        }

        if (!result.empty()) matchPoles(poles, result.back().poles);
        result.push_back({kappa, poles});
    }
    return result;
}

double loopGainMargin(const std::vector<RootLocusPoint>& locus) {
    auto maxRe = [](const RootLocusPoint& pt) {
        double m = -std::numeric_limits<double>::infinity();
        for (const auto& p : pt.poles) m = std::max(m, p.real());
        return m;
    };
    for (std::size_t k = 1; k < locus.size(); ++k) {
        const double a = maxRe(locus[k - 1]);
        const double b = maxRe(locus[k]);
        if (a < 0.0 && b >= 0.0) {
            const double t = -a / (b - a);
            return locus[k - 1].gain + t * (locus[k].gain - locus[k - 1].gain);
        }
    }
    return -1.0;
}

}  // namespace caliburn
```

- [ ] **Step 3: Rewrite the locus dispatch in `visualizer.cpp`**

Replace lines 149-164 (the `if (state.pz_mode == PZMode::UnityFB) ... else if ... StateFB`) with:

```cpp
        // Clear on EVERY non-computing path.  Today, when pz_mode == StateFB
        // and ctrl_type != GainMatrix, neither branch runs and root_locus is
        // never cleared, so the panel keeps drawing the PREVIOUS mode's locus
        // labelled as State FB, with its marker read off rl_current_alpha.
        // Clearing here kills that by construction rather than by a guard.
        state.root_locus.clear();
        state.rl_loop_gain_margin = -1.0;
        switch (state.pz_mode) {
            case caliburn::PZMode::PlantLocus:
                state.root_locus = caliburn::computeRootLocus(
                    state.plant, state.output_i, state.input_j,
                    state.rl_k_min, state.rl_k_max, state.rl_num_points);
                break;
            case caliburn::PZMode::LoopLocus:
                if (loop_backed && !state.loops.empty()) {
                    state.root_locus = caliburn::computeLoopLocus(
                        state.plant, state.loops, loop_kind,
                        state.output_i, state.input_j,
                        state.rl_kappa_min, state.rl_kappa_max,
                        state.rl_num_points);
                    state.rl_loop_gain_margin =
                        caliburn::loopGainMargin(state.root_locus);
                }
                break;
            case caliburn::PZMode::StateFB:
                if (state.ctrl_type == caliburn::ControllerType::GainMatrix &&
                    state.ctrl_K.size() > 0) {
                    state.root_locus = caliburn::computeStateFeedbackLocus(
                        state.plant, state.ctrl_K,
                        state.rl_alpha_min, state.rl_alpha_max,
                        state.rl_num_points);
                }
                break;
            case caliburn::PZMode::PoleZero:
                break;
        }
```

Add `#include "analysis/loop_locus.h"` — or rely on `app_state.h`, which does not pull it; add the include to `visualizer.cpp` explicitly.

- [ ] **Step 4: Four radios, κ controls and the margin readout**

In `src/panels/pole_zero_panel.cpp`, replace the mode selector (lines 28-36) and the control tables (lines 38-73):

```cpp
    int mode = static_cast<int>(state.pz_mode);
    ImGui::RadioButton("Pole-Zero", &mode, 0); ImGui::SameLine();
    ImGui::RadioButton("Plant Locus", &mode, 1); ImGui::SameLine();
    ImGui::RadioButton("Loop Locus", &mode, 2); ImGui::SameLine();
    ImGui::RadioButton("State FB", &mode, 3);
    if (static_cast<PZMode>(mode) != state.pz_mode) {
        state.pz_mode = static_cast<PZMode>(mode);
        state.needs_recompute = true;
    }

    // All four stay selectable.  Disabling them plus an auto-fallback would
    // silently yank the user out of their chosen mode on a controller-type
    // switch, and hiding them would change the control row's width as the
    // controller is edited with nothing on screen saying the mode exists.
    if (state.pz_mode == PZMode::PlantLocus) {
        if (ImGui::BeginTable("##rl_ctrls", 3)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##K", &state.rl_current_k,
                               state.rl_k_min, state.rl_k_max, "%.2f",
                               ImGuiSliderFlags_Logarithmic);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K min"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmin", &state.rl_k_min, 0.0f, 1.0f))
                state.needs_recompute = true;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("K max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##Kmax", &state.rl_k_max, 1.0f, 1000.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    } else if (state.pz_mode == PZMode::LoopLocus) {
        if (ImGui::BeginTable("##ll_ctrls", 3)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##kappa", &state.rl_kappa_current,
                               state.rl_kappa_min, state.rl_kappa_max, "%.3g",
                               ImGuiSliderFlags_Logarithmic);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba min"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            // Clamped away from zero: kappa = 0 would collapse the loop's
            // states and make the closed-loop dimension ragged mid-sweep.
            if (ImGui::SliderFloat("##kmin", &state.rl_kappa_min, 0.001f, 1.0f,
                                   "%.3g", ImGuiSliderFlags_Logarithmic))
                state.needs_recompute = true;
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xba max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##kmax", &state.rl_kappa_max, 1.0f, 1000.0f,
                                   "%.3g", ImGuiSliderFlags_Logarithmic))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
        if (state.rl_loop_gain_margin > 0.0) {
            ImGui::TextColored(ImVec4(0.3f, 1, 0.3f, 1),
                "\xce\xba* = %.3g   margin %.2f\xc3\x97 / %.1f dB",
                state.rl_loop_gain_margin,
                state.rl_loop_gain_margin / std::max(state.rl_kappa_current, 1e-6f),
                20.0 * std::log10(state.rl_loop_gain_margin /
                                  std::max(state.rl_kappa_current, 1e-6f)));
        } else if (!state.root_locus.empty()) {
            ImGui::TextDisabled("no crossing in [%.3g, %.3g]",
                                state.rl_kappa_min, state.rl_kappa_max);
        }
    } else if (state.pz_mode == PZMode::StateFB) {
        if (ImGui::BeginTable("##sfb_ctrls", 2)) {
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            ImGui::SliderFloat("##alpha", &state.rl_current_alpha,
                               state.rl_alpha_min, state.rl_alpha_max);
            ImGui::TableNextColumn();
            ImGui::TextUnformatted("\xce\xb1 max"); ImGui::SameLine();
            ImGui::SetNextItemWidth(-FLT_MIN);
            if (ImGui::SliderFloat("##alphamax", &state.rl_alpha_max, 0.1f, 10.0f))
                state.needs_recompute = true;
            ImGui::EndTable();
        }
    }

    // A locus mode with an unmet precondition draws the s-plane furniture and
    // states why, rather than the previous mode's data under the wrong label.
    if (state.pz_mode != PZMode::PoleZero && state.root_locus.empty()) {
        const char* why =
            (state.pz_mode == PZMode::LoopLocus)
                ? "Loop Locus needs a PID or Lead/Lag controller with at least "
                  "one loop on the selected channel"
          : (state.pz_mode == PZMode::StateFB)
                ? "State FB needs a Gain Matrix K controller"
                : "no locus data";
        ImGui::TextDisabled("%s", why);
    }
```

- [ ] **Step 5: Fix the current-value marker and the grey reference**

In the drawing branch, replace the marker selection (lines 195-197):

```cpp
                double current_gain = state.rl_current_k;
                if (state.pz_mode == PZMode::LoopLocus)
                    current_gain = state.rl_kappa_current;
                else if (state.pz_mode == PZMode::StateFB)
                    current_gain = state.rl_current_alpha;
```

and replace the grey `OL Poles` overlay (lines 217-228) — which reads `state.pole_zero[0]`, a plant-shaped assumption of exactly the kind rule C exists to remove — with the sweep's own first point:

```cpp
            // One rule for both locus modes: the reference is where the sweep
            // actually starts.  For Plant Locus with k_min = 0 that IS the
            // plant's poles, so nothing changes there.  For Loop Locus the
            // plant's poles are simply the wrong reference — the branches start
            // near the poles of plant-plus-other-loops-plus-this-compensator.
            if (!state.root_locus.empty()) {
                std::vector<double> pre, pim;
                for (const auto& pole : state.root_locus.front().poles) {
                    pre.push_back(pole.real());
                    pim.push_back(pole.imag());
                }
                char lbl[48];
                std::snprintf(lbl, sizeof(lbl), "Start (%s = %.3g)",
                              state.pz_mode == PZMode::LoopLocus ? "\xce\xba"
                              : state.pz_mode == PZMode::StateFB ? "\xce\xb1"
                                                                 : "K",
                              state.root_locus.front().gain);
                ImPlot::PlotScatter(lbl, pre.data(), pim.data(),
                                   static_cast<int>(pre.size()),
                                   ImPlotSpec(ImPlotProp_Marker, ImPlotMarker_Cross,
                                              ImPlotProp_MarkerSize, 6.0f,
                                              ImPlotProp_MarkerLineColor,
                                              ImVec4(0.5f, 0.5f, 0.5f, 0.7f),
                                              ImPlotProp_LineWeight, 1.5f));
            }
```

Also update the axis-limit pass at lines 89-94 to use `root_locus.front()` rather than `pole_zero[0]` for the same reason.

- [ ] **Step 6: Correct the help text per mode**

The current help promises *"Start at OL poles (K=0), end at OL zeros (K→∞)"* while **no locus mode draws zeros at all** — a pre-existing lie that two new modes would otherwise inherit. Replace the last three lines of the `drawHelpMarker` string with:

```
        "Pole-Zero: poles and zeros of each visible system.\n"
        "Plant Locus: proportional gain K on the bare plant channel (i,j).\n"
        "  The controller is NOT in the path — this is the pre-tuning tool.\n"
        "Loop Locus: kappa scaling the whole compensator on channel (i,j),\n"
        "  a true locus of 1 + kappa*L(s) = 0. kappa = 1 is the live design.\n"
        "State FB: alpha scaling the gain matrix K.\n\n"
        "No mode draws zeros. Unstable where a branch crosses into Re > 0."
```

- [ ] **Step 7: Wire into CMake, build, test, commit**

Add `src/analysis/loop_locus.cpp` to `analysis_lib`.

```bash
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build --output-on-failure
git add src/analysis/loop_locus.h src/analysis/loop_locus.cpp \
        src/visualizer.cpp src/panels/pole_zero_panel.cpp CMakeLists.txt
git commit -m "feat: Loop Locus — kappa on the compensator — and a 4-way PZMode

The root locus responds to Tier 1 tuning for the first time: kappa scales the
whole compensator on rule C's channel, re-running the production assembly at
each point so the locus cannot disagree with the closed-loop pole plot beside
it. Reports the loop gain margin from the sweep it already has. Clearing
root_locus on every non-computing path fixes a live bug where State FB
without a K redrew the previous mode's locus under the wrong label.

Refs #10"
```

---

### Task 8: The pairing grid

**Files:**
- Modify: `src/app_state.h`
- Modify: `src/visualizer.cpp`
- Modify: `src/panels/model_panel.cpp`

The p×m matrix **is** the editor: the same widget that grades the pairing is the one that creates it. Click a cell to pair, click again to unpair.

**354 px is the binding constraint** — the width of the docked Model Configuration panel in the checked-in `imgui.ini`. It eliminates the ledger layout (a table row per loop with gains as columns): at 354 px the `Kp`/`Ki`/`Kd`/`τf` headers degrade to `…` and the fields are one character wide. Adopting it would mean changing the app's default dock layout to buy a table. The reference implementation is `src/panels/loop_ui_prototype.cpp` on branch `prototype/loop-list-ui` (`de70a09`), variant C, with screenshots in `docs/prototypes/loop-ui/`.

- [ ] **Step 1: Cache the diagnostics in `AppState`** **[plan-level]**

The dead-channel test sweeps the whole Bode grid per loop, so it must not run per frame. Add to `AppState` beside the analysis results:

```cpp
    LoopDiagnostics diagnostics;
```

and `#include "analysis/loop_diagnostics.h"` to `app_state.h`. In `visualizer.cpp`, after the liveness block of Task 5:

```cpp
        state.diagnostics = caliburn::computeLoopDiagnostics(
            state.plant, state.loops, state.output_scales,
            state.diag_omega_hz,
            state.freq_min_hz, state.freq_max_hz, state.num_freq_points);
```

Every control that feeds it — the ω slider, the scales, a cell click — sets `needs_recompute`.

- [ ] **Step 2: Add the piecewise gain slider helper**

In the anonymous namespace of `src/panels/model_panel.cpp`:

```cpp
// Ki and Kd must reach exactly 0 and still resolve small values, and the spec
// asks for "linear 0—1, then log 1—hi".  ImGuiSliderFlags_Logarithmic DOES work
// with v_min = 0 — verified in the prototype's bench, no assert, no crash,
// ImGui's zero epsilon handles it — but it is rejected on emphasis: a value of
// 0.5 sits at ~62% of travel, compressing the 1—hi decades into the last third,
// the inverse of the split the spec asks for.  The piecewise scale below puts
// 0.5 at exactly 25%.
bool gainSliderPiecewise(const char* id, float* v, float hi) {
    const float lin_hi = 1.0f;
    float t = (*v <= lin_hi)
                  ? 0.5f * (*v / lin_hi)
                  : 0.5f + 0.5f * std::log(*v / lin_hi) / std::log(hi / lin_hi);
    char fmt[32];
    std::snprintf(fmt, sizeof(fmt), "%.4g", *v);   // display the true value
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::SliderFloat(id, &t, 0.0f, 1.0f, fmt)) {
        *v = (t <= 0.5f) ? 2.0f * t * lin_hi
                         : lin_hi * std::pow(hi / lin_hi, (t - 0.5f) * 2.0f);
        return true;
    }
    return false;
}

// Severity is graded from the RGA number, never from |λ_ii|, and applies to RGA
// cells only — Channel Share is never coloured.
ImVec4 severityFill(double rga_number) {
    if (rga_number < 0.5) return ImVec4(0.12f, 0.32f, 0.16f, 1.0f);
    if (rga_number < 2.0) return ImVec4(0.36f, 0.30f, 0.10f, 1.0f);
    return ImVec4(0.40f, 0.14f, 0.14f, 1.0f);
}
```

- [ ] **Step 3: Replace the controller section**

Replace `src/panels/model_panel.cpp` lines 282-358 (from `ImGui::SeparatorText("Controller");` to the close of the `GainMatrix` branch):

```cpp
    // --- Controller section ---
    ImGui::SeparatorText("Controller");
    const char* ctrl_types[] = {"None", "PID", "Lead/Lag",
                                "State-Space", "Gain Matrix K"};
    int ctrl_idx = static_cast<int>(state.ctrl_type);
    if (ImGui::Combo("Type", &ctrl_idx, ctrl_types, 5)) {
        state.ctrl_type = static_cast<ControllerType>(ctrl_idx);
        state.needs_recompute = true;
        const bool loop_backed = state.ctrl_type == ControllerType::PID ||
                                 state.ctrl_type == ControllerType::LeadLag;
        state.trace_visible[3] = state.ctrl_type != ControllerType::None;
        state.trace_visible[1] =
            loop_backed || state.ctrl_type == ControllerType::StateSpace;
    }

    const bool loop_backed = state.ctrl_type == ControllerType::PID ||
                             state.ctrl_type == ControllerType::LeadLag;

    if (loop_backed) {
        const int gp = state.plant.outputs();
        const int gm = state.plant.inputs();
        const auto& diag = state.diagnostics;

        ImGui::TextDisabled("click a cell to pair / unpair");

        // Columns: per-output scale, row label, then one cell per plant input.
        // The scale column supersedes a collapsed "scales: default (1, 1)"
        // status line — the values are on screen, so nothing has to state
        // whether they are still at default, and the row header is where
        // per-output metadata would naturally go if that fog ever clears.
        if (ImGui::BeginTable("##pairgrid", gm + 2,
                              ImGuiTableFlags_Borders |
                              ImGuiTableFlags_SizingFixedFit)) {
            ImGui::TableSetupColumn("scale", ImGuiTableColumnFlags_WidthFixed, 56);
            ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthFixed, 24);
            for (int j = 0; j < gm; ++j) {
                char h[8];
                std::snprintf(h, sizeof(h), "u%d", j);
                ImGui::TableSetupColumn(h, ImGuiTableColumnFlags_WidthFixed, 74);
            }
            ImGui::TableHeadersRow();

            for (int i = 0; i < gp; ++i) {
                ImGui::PushID(i);
                ImGui::TableNextRow();

                ImGui::TableNextColumn();
                if (i < (int)state.output_scales.size()) {
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    if (ImGui::DragFloat("##sc", &state.output_scales[i], 0.01f,
                                         0.001f, 1000.0f, "%.3g"))
                        state.needs_recompute = true;
                }

                ImGui::TableNextColumn();
                ImGui::Text("y%d", i);

                for (int j = 0; j < gm; ++j) {
                    ImGui::TableNextColumn();
                    ImGui::PushID(j);

                    int found = -1;
                    for (std::size_t k = 0; k < state.loops.size(); ++k)
                        if (state.loops[k].out == i && state.loops[k].in == j) {
                            found = static_cast<int>(k);
                            break;
                        }
                    const bool paired = found >= 0;

                    // Two lines: the diagnostic over |g_ij| in dB, dimmed.  The
                    // dB line is REQUIRED under Channel Share, which cannot
                    // tell "drives both strongly" from "drives neither"; under
                    // the RGA it is free and still worth reading.
                    char top[24] = "-";
                    if (diag.has_rga) {
                        const auto ro = std::find(diag.sub_out.begin(),
                                                  diag.sub_out.end(), i);
                        const auto co = std::find(diag.sub_in.begin(),
                                                  diag.sub_in.end(), j);
                        if (ro != diag.sub_out.end() && co != diag.sub_in.end())
                            std::snprintf(top, sizeof(top), "%.2f",
                                diag.lambda(ro - diag.sub_out.begin(),
                                            co - diag.sub_in.begin()));
                    } else if (diag.has_share && j == diag.share_in &&
                               i < (int)diag.share.size()) {
                        std::snprintf(top, sizeof(top), "%.0f%%",
                                      100.0 * diag.share[i]);
                    }

                    char cell[64];
                    if (diag.has_share && i < (int)diag.mag_db.size() &&
                        j == diag.share_in) {
                        std::snprintf(cell, sizeof(cell), "%s\n%.1f dB##c",
                                      top, diag.mag_db[i]);
                    } else {
                        std::snprintf(cell, sizeof(cell), "%s\n##c", top);
                    }

                    // Paired-ness is the BORDER; severity is the FILL.  One
                    // encoding cannot carry both — the prototype had only fill
                    // and so could not show severity at all.
                    const bool colour_cell = paired && diag.has_rga;
                    ImGui::PushStyleColor(ImGuiCol_Button,
                        colour_cell ? severityFill(diag.rga_number)
                                    : ImVec4(0.16f, 0.16f, 0.18f, 1.0f));
                    ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize,
                                        paired ? 2.0f : 0.0f);
                    ImGui::PushStyleColor(ImGuiCol_Border,
                                          ImVec4(0.35f, 0.72f, 1.0f, 1.0f));

                    if (ImGui::Button(cell, ImVec2(-FLT_MIN, 34))) {
                        if (paired) {
                            state.loops.erase(state.loops.begin() + found);
                            state.selected_loop =
                                state.loops.empty() ? -1 : 0;
                        } else {
                            // Clicking an unpaired cell creates the loop and
                            // selects it.  Selection is the last cell clicked.
                            state.loops.push_back(Loop{i, j, {}, {}});
                            state.selected_loop =
                                static_cast<int>(state.loops.size()) - 1;
                        }
                        state.needs_recompute = true;
                    } else if (paired && ImGui::IsItemClicked()) {
                        state.selected_loop = found;
                    }

                    ImGui::PopStyleColor(2);
                    ImGui::PopStyleVar();

                    if (paired && found < (int)diag.loop_dead.size() &&
                        diag.loop_dead[found]) {
                        ImGui::SameLine(0, 2);
                        ImGui::TextColored(ImVec4(1, 0.35f, 0.35f, 1), "!");
                        if (ImGui::IsItemHovered())
                            ImGui::SetTooltip(
                                "u%d does not affect y%d at any frequency", j, i);
                    }
                    ImGui::PopID();
                }
                ImGui::PopID();
            }
            ImGui::EndTable();
        }

        // One omega, plus a collapsed sweep.  No separate steady-state
        // frequency: S&P rule 2 is a DIC theorem requiring a stable plant, and
        // neither Ball-Balancer nor Inverted Pendulum qualifies, so a
        // permanently-inapplicable readout would earn its space on no preset.
        ImGui::TextUnformatted("\xcf\x89 [Hz]"); ImGui::SameLine();
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::SliderFloat("##diagw", &state.diag_omega_hz, 0.01f, 100.0f,
                               "%.3g", ImGuiSliderFlags_Logarithmic))
            state.needs_recompute = true;

        ImGui::Checkbox("RGA number vs \xcf\x89", &state.diag_sweep);
        if (state.diag_sweep) {
            float vals[40];
            for (int k = 0; k < 40; ++k) {
                const double f = state.freq_min_hz *
                    std::pow(state.freq_max_hz / state.freq_min_hz, k / 39.0);
                const double v = rgaNumberAt(state.plant, state.loops, f);
                vals[k] = std::isnan(v) ? 0.0f : static_cast<float>(v);
            }
            ImGui::PlotLines("##sweep", vals, 40, 0, nullptr, 0.0f, FLT_MAX,
                             ImVec2(-FLT_MIN, 46));
        }

        // The subtitle carries the definition AND the explicit negative.  The
        // failure mode being guarded against is a user reading a number in the
        // pairing area and assuming it grades their pairing.
        ImGui::TextWrapped("%s", state.diagnostics.headline.c_str());
        if (state.diagnostics.has_share) {
            ImGui::TextDisabled(
                "Channel Share is not the RGA and does not grade the pairing.");
        }

        // Gain block for the selected loop, full width, below the grid.
        if (state.selected_loop >= 0 &&
            state.selected_loop < (int)state.loops.size()) {
            Loop& l = state.loops[state.selected_loop];
            char hdr[64];
            std::snprintf(hdr, sizeof(hdr), "Gains  y%d \xe2\x86\x90 u%d",
                          l.out, l.in);
            ImGui::SeparatorText(hdr);

            bool g_changed = false;
            if (state.ctrl_type == ControllerType::PID) {
                ImGui::TextUnformatted("Kp"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##Kp", &l.pid.Kp, 0.01f, 100.0f,
                                                "%.4g", ImGuiSliderFlags_Logarithmic);
                ImGui::TextUnformatted("Ki"); ImGui::SameLine();
                g_changed |= gainSliderPiecewise("##Ki", &l.pid.Ki, 50.0f);
                ImGui::TextUnformatted("Kd"); ImGui::SameLine();
                g_changed |= gainSliderPiecewise("##Kd", &l.pid.Kd, 20.0f);
                ImGui::TextUnformatted("\xcf\x84f"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##tauf", &l.pid.tau_f,
                                                0.001f, 1.0f, "%.4g",
                                                ImGuiSliderFlags_Logarithmic);
            } else {
                // Kind is one top-level combo for the whole list; Lead/Lag MODE
                // is per loop, inside the gain block.
                const char* modes[] = {"Lead", "Lag", "Lead-Lag"};
                int mi = static_cast<int>(l.leadlag.mode);
                if (ImGui::Combo("Mode", &mi, modes, 3)) {
                    l.leadlag.mode = static_cast<CompensatorMode>(mi);
                    g_changed = true;
                }
                ImGui::TextUnformatted("Kc"); ImGui::SameLine();
                ImGui::SetNextItemWidth(-FLT_MIN);
                g_changed |= ImGui::SliderFloat("##Kc", &l.leadlag.Kc, 0.1f, 100.0f,
                                                "%.4g", ImGuiSliderFlags_Logarithmic);
                if (l.leadlag.mode != CompensatorMode::Lag) {
                    ImGui::TextUnformatted("\xce\xb1 lead"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    g_changed |= ImGui::SliderFloat("##al", &l.leadlag.alpha_lead,
                                                    0.01f, 0.99f, "%.4g",
                                                    ImGuiSliderFlags_Logarithmic);
                    ImGui::TextUnformatted("T lead"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    g_changed |= ImGui::SliderFloat("##Tl", &l.leadlag.T_lead,
                                                    0.001f, 100.0f, "%.4g",
                                                    ImGuiSliderFlags_Logarithmic);
                }
                if (l.leadlag.mode != CompensatorMode::Lead) {
                    ImGui::TextUnformatted("\xce\xb1 lag"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    g_changed |= ImGui::SliderFloat("##ag", &l.leadlag.alpha_lag,
                                                    1.01f, 100.0f, "%.4g",
                                                    ImGuiSliderFlags_Logarithmic);
                    ImGui::TextUnformatted("T lag"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(-FLT_MIN);
                    g_changed |= ImGui::SliderFloat("##Tg", &l.leadlag.T_lag,
                                                    0.001f, 100.0f, "%.4g",
                                                    ImGuiSliderFlags_Logarithmic);
                }
            }
            if (g_changed) state.needs_recompute = true;
        } else {
            ImGui::TextDisabled("select a paired cell to edit its gains");
        }
    } else if (state.ctrl_type == ControllerType::StateSpace) {
        // ... unchanged from today ...
    } else if (state.ctrl_type == ControllerType::GainMatrix) {
        // ... unchanged from today ...
    }
```

Keep the existing `StateSpace` and `GainMatrix` bodies verbatim — only the enum comparison changes.

**No add / remove / reorder controls exist**, and none should be added. Cells are the affordance, grid order is loop order, so reorder is meaningless. One consequence: the grid cannot express an exact duplicate pairing, so the duplicate-pairing caution needs no UI surface — the builder still sums duplicates, which is correct for a loop list built programmatically and costs nothing. **All fan-out shapes remain expressible**: many-to-one is several cells in one row, one-to-many several cells in one column.

Add `#include "../analysis/loop_diagnostics.h"` and `<cmath>` where needed.

- [ ] **Step 4: Build and commit**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
git add src/app_state.h src/visualizer.cpp src/panels/model_panel.cpp
git commit -m "feat: pairing grid — the diagnostic matrix is the loop editor

The p x m grid creates and destroys pairings by click, so there are no add,
remove or reorder controls and grid order is loop order. Ball-Balancer's dead
identity pairing becomes a picture (paired cells 0.00, unpaired 1.00) and
Quarter-Car degrades to a single Channel Share column for free — one panel
slot, two mutually exclusive readouts. Ki/Kd ship a hand-rolled piecewise
scale; ImGuiSliderFlags_Logarithmic works at v_min = 0 but inverts the
emphasis the spec asks for.

Refs #6, #8"
```

---

### Task 9: Correct the design spec

**Files:**
- Modify: `docs/specs/2026-07-23-controllers-tier1-pid-leadlag-design.md`

- [ ] **Step 1: Add a superseded banner**

Insert after the `# Tier 1: ...` heading:

```markdown
> **Superseded in part.** The MIMO half of this spec — `std::vector<PIDParams>`
> sized by `plant.inputs()`, `buildDiagonalPID`, and the DC-gain coupling
> heuristic — assumed square plants and does not hold: only 3 of the 6 presets
> are square, and the one interesting multivariable preset has a singular `A`,
> so `−C·A⁻¹·B + D` is undefined for it. The loop-pairing model replaces it; see
> `docs/plans/2026-08-30-controllers-tier1-loop-pairing.md`.
>
> **Still correct and reused verbatim:** the state-space realization tables for
> all four PID forms and both Lead/Lag forms.
```

- [ ] **Step 2: Replace the Verification section (lines 305-315)**

Step 8 of the current list — *"Select a MIMO preset (Quarter-Car)"* — is unachievable: Quarter-Car is 1-input / 2-output. Replace the whole section with:

```markdown
## Verification

1. `cmake --build build && ctest --test-dir build` — clean compile, all suites
   green including `test_controller_builders` and `test_loop_diagnostics`.
2. **First-Order** (1×1): pick PID. The grid is a single cell, already paired by
   the identity seed. Sliding Kp shifts the Bode magnitude and moves the
   closed-loop pole.
3. Raise Ki — the step response's steady-state error goes to zero, and the
   closed-loop pole count rises by one as the integrator state appears.
4. Raise Kd — overshoot reduces and the phase margin improves.
5. **The n = 0 path.** Set Ki = Kd = 0 exactly. The controller contributes no
   states: the closed-loop pole count equals the plant's, the Controller trace
   is flat at 20·log₁₀(Kp) with no poles plotted, and neither build config
   faults.
6. **Lead/Lag** on First-Order: in Lead mode, sweeping α and T shows phase lead;
   switching to Lag shows the low-frequency gain boost.
7. **Ball-Balancer** (2×2 — the only genuinely multivariable preset): the seeded
   identity pairing reads λ = 0.00 on both paired cells and 1.00 on both
   unpaired ones, with an RGA number of 4.00 at any ω > 0 and a dead-channel
   warning on each paired cell. Clicking the two off-diagonal cells takes the
   RGA number to 0.00 and clears the warnings.
8. **Quarter-Car** (1-input / 2-output — SIMO, *not* MIMO): the grid is a single
   column and the readout is **Channel Share**, never labelled or coloured as
   RGA, each row showing the share over |gᵢ| in dB. Sweeping ω from 1 Hz to
   10 Hz moves the share from `zs` toward `zu` — body bounce to wheel hop.
   Changing an output scale changes the shares and the headline stops reading
   "at default".
9. **Inverted Pendulum** (SIMO): add a second loop into the same input, so both
   `y0 → u0` and `y1 → u0` are paired — the MISO fan-in. The controller is 1×2,
   the closed loop is still 2×2, and both compensator channels are live.
10. **Rule C.** On Quarter-Car with only `y0` paired, selecting Output i = 1
    suppresses the Controller trace with *"no loop pairs y1 → u0"* rather than
    drawing a flat 0° phase; setting Reference r = 1 suppresses the closed-loop
    trace for the matching reason.
11. **Loop Locus.** On Second-Order with a PID loop, the branches start at
    κ = 0.01 with the grey reference on the sweep's first point, the marker sits
    at κ = 1, and the readout gives κ* or *"no crossing in [0.01, 100]"*.
    Selecting State FB with no K applied clears the plot and states why, rather
    than redrawing the previous mode's locus.
12. **Grid mode.** Ball-Balancer with Show All Channels draws two tables — open
    p×m and loop p×p — and no cell shows another channel's data.
13. Existing State-Space and Gain Matrix controllers still work unchanged.
```

- [ ] **Step 3: Commit**

```bash
git add docs/specs/2026-07-23-controllers-tier1-pid-leadlag-design.md
git commit -m "docs: correct the Tier 1 spec's verification list and mark the MIMO half superseded

Step 8 named Quarter-Car as a MIMO preset; it is 1-input / 2-output. The new
list is written against the presets that actually exist, and covers the n=0
path, the MISO fan-in and rule C's suppression behaviour.

Refs #7"
```

---

### Task 10: Visual verification

Run `./build/visualizer` and work steps 2-13 of the corrected Verification list above, in order. Everything there is reachable from the shipped presets — that is the point of rewriting it.

- [ ] **Step 1: Launch and confirm the build is optimized**

```bash
grep CMAKE_BUILD_TYPE build/CMakeCache.txt   # expect Release
./build/visualizer
```

Drag a PID gain slider on Ball-Balancer with Loop Locus selected and confirm the UI stays responsive — the 200-point sweep is ~2.4 ms at ten closed-loop states, so no throttling or debouncing is needed anywhere.

- [ ] **Step 2: Work the Verification list**

Steps 2-13. The three that are most likely to expose a mistake, and what a mistake looks like:

| Check | A wrong implementation shows |
|---|---|
| Ball-Balancer RGA number | Anything other than 4.00 for identity / 0.00 for the swap means `.adjoint()` crept in, or the sub-matrix is being taken from the wrong rows |
| Quarter-Car readout | A green `1.00` labelled RGA means the trigger was written as `plant.inputs() == 1` or `loops.size() >= 2` instead of the square-and-≥2 sub-matrix test |
| P-only loop, Controller trace | A pole at the origin or at −1/τf means the minimal-realization rule was not applied — a phantom mode |

- [ ] **Step 3: Confirm the panel is legible at its real width**

Undock nothing. The Model Configuration panel is 354 px in the checked-in `imgui.ini`; if the grid needs more, the grid is wrong, not the layout.

- [ ] **Step 4: Check the web build still configures**

```bash
emcmake cmake -S . -B build-web && cmake --build build-web -j8
```

`analysis_lib` is not gated on `if(NOT EMSCRIPTEN)`, so `controller_builders.cpp`, `loop_diagnostics.cpp` and `loop_locus.cpp` join the web build automatically. The panel/UI side is the part that was never verified — if this fails, it fails there, and that is the open fog item below, not a regression from this plan.

---

## Deliberately not in this plan

Carried from the map's **Out of scope**, so they are not re-proposed:

- **Derived `ctrl_ss` display and a "bake to State-Space" escape hatch** — read-only Ac/Bc/Cc/Dc for the built controller. Ruled out; an implementer may hide the sliders by default instead.
- **A Tier 2/3 forward-compatibility seam** in `ControllerType` / `AppState` for pole placement, LQR and eigenstructure assignment. Deferred to those tiers, accepting a possible refactor.
- **Anti-windup and actuator saturation** — nonlinear, beyond what a linear analyzer models.

And from the ticket record:

- **Filtering the controller's pole plot to the selected loop's block.** More truthful, but it needs `buildLoopController` to report per-loop state ranges — a new return contract — and would make the controller and the plant behave differently for the same plot. A cheap follow-up if it grates once built, not a blocker.
- **A per-loop enable/bypass flag.** Muting a loop is already expressible as zero gains, and the field would need its own semantics in the assembler.
- **An independent loop selector for the locus.** The swept loop is derived from rule C's channel, so there is no second answer to "which loop is dead" that must stay in step with the Bode panel's.

## Still open

- **Per-output metadata** (name, unit, allowed deviation) on `ModelEntry`. Outputs are anonymous indices today and the user-set scale vector covers the immediate need, but this would buy output *names* app-wide — plot legends, channel selectors, the pairing grid's row headers, where the scale `DragFloat` already sits beside the `y0` label. Still blocked on how the free-text `C`-editing path, where `p` changes under you, would carry it.
- **Emscripten/web parity for the panel/UI sources.** The analysis sources are covered by `analysis_lib` being ungated; the visualizer/panel side is unverified. Task 10 Step 4 is a check, not a fix.

## Files touched

| File | Change |
|---|---|
| `src/analysis/pole_zero.h` / `.cpp` | Two `n == 0` guards; `matchPoles` promoted out of the anonymous namespace |
| `src/analysis/frequency_response.cpp` | `n == 0` early return in `evalTransferFunction` |
| `src/analysis/controller_builders.h` / `.cpp` | **New.** `PIDParams`, `LeadLagParams`, `CompensatorMode`, `Loop`, `LoopKind`; `buildPID`, `buildLeadLag`, `buildLoopController` |
| `src/analysis/loop_diagnostics.h` / `.cpp` | **New.** RGA on the square paired sub-matrix, Channel Share elsewhere, structurally-dead-channel test |
| `src/analysis/loop_locus.h` / `.cpp` | **New.** κ-on-compensator locus and the loop gain margin |
| `src/app_state.h` | `loops`, `cached_p/m`, `output_scales`, `selected_loop`, `ref_r`, `channel_live`/`channel_reason`, `diag_*`, `rl_kappa_*`, `diagnostics`; `ControllerType` and `PZMode` extended; `time_input_j` deleted; `seedIdentityLoops` |
| `src/visualizer.cpp` | Dimension guard and seeding; loop-backed controller assembly; rule-C channel map and liveness; two-role channel grid; time indices; 4-way locus dispatch with a clear on every non-computing path |
| `src/panels/model_panel.cpp` | Controller section becomes the pairing grid; `ref_r` slider; piecewise Ki/Kd scale |
| `src/panels/bode_panel.cpp` | Two role-grids replacing the plant-dimensioned single grid; stated reasons for suppressed channels |
| `src/panels/pole_zero_panel.cpp` | Four radios, κ controls, gain-margin readout, `root_locus.front()` reference, per-mode help text |
| `src/panels/time_response_panel.cpp` | `time_input_j` slider retired in favour of `input_j` / `ref_r` |
| `CMakeLists.txt` | Release by default; three new `analysis_lib` sources; two new test targets |
| `tests/test_controller_builders.cpp` | **New.** |
| `tests/test_loop_diagnostics.cpp` | **New.** |
| `docs/specs/2026-07-23-controllers-tier1-pid-leadlag-design.md` | Superseded banner; Verification section rewritten |
