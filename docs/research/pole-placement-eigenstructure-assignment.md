# Pole Placement, Eigenvalue Assignment, and Eigenstructure Assignment

Research for the Caliburn Linear System Analyzer. All algorithms target C++17 with Eigen 3.4, no external control libraries.

---

## 1. Controllability Check (Prerequisite)

Already implemented in `src/analysis/system_properties.cpp` as `checkControllability()`.

**What exists:**
```
Mc = [B, AB, A^2 B, ..., A^(n-1) B]    // n x (n*m) matrix
rank(Mc) == n  =>  controllable
```

**What needs to change for pole placement:** Nothing. The existing `PropertyResult` returns `rank`, `required_rank`, and `pass`. Before any placement attempt, call `checkControllability(plant)` and abort if `!pass`.

For the observer design (Section 4), the dual check is `checkObservability()`, also already implemented.

---

## 2. Pole Placement via Ackermann's Formula (SISO)

### When It Applies

- Single-input systems only (m = 1, `plant.inputs() == 1`)
- Controllable system (rank(Mc) = n)

### User Inputs

- Desired closed-loop pole locations: vector of n complex numbers (conjugate pairs for real K)
- The plant (A, B, C, D)

### Algorithm

Given desired poles {p_1, p_2, ..., p_n}, form the desired characteristic polynomial:

```
alpha(s) = (s - p_1)(s - p_2)...(s - p_n) = s^n + a_{n-1} s^{n-1} + ... + a_1 s + a_0
```

**Step 1:** Compute polynomial coefficients from desired poles.

Given poles as complex numbers, expand the product. Since the gain K must be real, the poles must come in conjugate pairs (or be real). The expansion produces real coefficients {a_0, a_1, ..., a_{n-1}}.

**Step 2:** Evaluate the characteristic polynomial at the matrix A:

```
alpha(A) = A^n + a_{n-1} A^{n-1} + ... + a_1 A + a_0 I
```

This is a matrix polynomial evaluation. Compute iteratively:

```cpp
MatrixXd phi = MatrixXd::Identity(n, n) * coeffs[0];  // a_0 * I
MatrixXd Ak = MatrixXd::Identity(n, n);
for (int k = 1; k <= n; ++k) {
    Ak = Ak * A;  // A^k
    phi += coeffs[k] * Ak;
}
// coeffs[n] = 1.0 (monic polynomial)
```

**Step 3:** Compute controllability matrix Mc = [B, AB, ..., A^{n-1}B].

Already available from `checkControllability()`.

**Step 4:** Compute the last row of Mc^{-1}.

For a single-input system, Mc is n x n (square). The gain vector is:

```
K = e_n^T * Mc^{-1} * alpha(A)
```

where e_n^T = [0, 0, ..., 0, 1] is the last standard basis row vector.

Rather than inverting Mc explicitly, solve `e_n^T = x^T * Mc` for x, then `K = x^T * alpha(A)`:

```cpp
// Solve Mc^T * x = e_n  (column vector form)
VectorXd e_n = VectorXd::Zero(n);
e_n(n - 1) = 1.0;
VectorXd x = Mc.transpose().colPivHouseholderQr().solve(e_n);
// K = x^T * phi(A)  => K is a 1 x n row vector
RowVectorXd K = x.transpose() * phi_A;
```

**Step 5:** Form closed-loop system: `A_cl = A - B * K`

This is already implemented as `stateFeedbackClose()` in `system_connect.cpp`.

### Output

- Gain matrix K (1 x n row vector)
- Closed-loop system via `stateFeedbackClose(plant, K)`

### Eigen Implementation

```cpp
#include <Eigen/Core>
#include <Eigen/QR>
#include <complex>
#include <vector>

namespace caliburn {

struct PolePlacementResult {
    Eigen::MatrixXd K;          // m x n gain matrix
    bool success;
    std::string error;          // Human-readable failure reason
};

// Expand (s - p_1)(s - p_2)...(s - p_n) into real polynomial coefficients.
// Returns coefficients [a_0, a_1, ..., a_{n-1}, 1] (length n+1, monic).
// Requires conjugate-paired complex roots for real coefficients.
std::vector<double> polyFromRoots(
    const std::vector<std::complex<double>>& roots)
{
    int n = static_cast<int>(roots.size());
    // Start with p(s) = 1
    std::vector<std::complex<double>> coeffs(n + 1, {0.0, 0.0});
    coeffs[0] = {1.0, 0.0};

    for (int i = 0; i < n; ++i) {
        // Multiply current polynomial by (s - roots[i])
        // New coeffs[k] = coeffs[k-1] - roots[i] * coeffs[k]
        for (int k = i + 1; k > 0; --k) {
            coeffs[k] = coeffs[k - 1] - roots[i] * coeffs[k];
        }
        coeffs[0] = -roots[i] * coeffs[0];
    }

    // coeffs is in descending order: coeffs[0]*s^n + ... + coeffs[n]
    // Reverse to ascending: [a_0, a_1, ..., a_{n-1}, 1]
    std::vector<double> real_coeffs(n + 1);
    for (int i = 0; i <= n; ++i) {
        real_coeffs[i] = coeffs[n - i].real();
        // Imaginary parts should be ~0 if roots are conjugate-paired
    }
    return real_coeffs;
}

// Evaluate matrix polynomial: alpha(A) = sum_k coeffs[k] * A^k
Eigen::MatrixXd polyEvalMatrix(
    const Eigen::MatrixXd& A,
    const std::vector<double>& coeffs)
{
    int n = A.rows();
    int degree = static_cast<int>(coeffs.size()) - 1;
    Eigen::MatrixXd result = coeffs[0] * Eigen::MatrixXd::Identity(n, n);
    Eigen::MatrixXd Ak = Eigen::MatrixXd::Identity(n, n);
    for (int k = 1; k <= degree; ++k) {
        Ak = Ak * A;
        result += coeffs[k] * Ak;
    }
    return result;
}

PolePlacementResult ackermannPlace(
    const LinearSystem& plant,
    const std::vector<std::complex<double>>& desired_poles)
{
    PolePlacementResult result;
    int n = plant.states();
    int m = plant.inputs();

    if (m != 1) {
        result.success = false;
        result.error = "Ackermann's formula requires single-input (m=1)";
        return result;
    }

    if (static_cast<int>(desired_poles.size()) != n) {
        result.success = false;
        result.error = "Number of desired poles must equal number of states";
        return result;
    }

    // Check controllability
    auto ctrb = checkControllability(plant);
    if (!ctrb.pass) {
        result.success = false;
        result.error = "System is not controllable";
        return result;
    }

    // Mc is n x n for SISO (since m=1)
    const Eigen::MatrixXd& Mc = ctrb.matrix;

    // Characteristic polynomial from desired poles
    auto coeffs = polyFromRoots(desired_poles);

    // Evaluate polynomial at A
    Eigen::MatrixXd phi_A = polyEvalMatrix(plant.A, coeffs);

    // Solve for K = e_n^T * Mc^{-1} * phi(A)
    Eigen::VectorXd e_n = Eigen::VectorXd::Zero(n);
    e_n(n - 1) = 1.0;
    Eigen::VectorXd x = Mc.transpose().colPivHouseholderQr().solve(e_n);

    result.K = x.transpose() * phi_A;  // 1 x n
    result.success = true;
    return result;
}

}  // namespace caliburn
```

### Numerical Considerations

- **Condition of Mc:** For high-order systems (n > ~10), the controllability matrix becomes very ill-conditioned. Matrix powers A^k grow or shrink exponentially. For the expected use case (n <= 8, educational models), this is not a problem.
- **Conjugate pairs:** The UI must enforce that complex poles come in conjugate pairs. A simple approach: let the user enter poles as either real values or (real, imag) pairs, and auto-add the conjugate.
- **Verification:** After computing K, verify by checking eigenvalues of (A - B*K) match the desired poles. This is already computable with `Eigen::EigenSolver`.

### Limitations

- SISO only (m = 1)
- Numerically fragile for n > 10
- Does not allow shaping eigenvectors (only eigenvalues)

---

## 3. Eigenvalue Assignment (MIMO Pole Placement)

### The Problem

For multi-input systems (m > 1), Ackermann's formula does not apply. We need algorithms that handle the additional degrees of freedom from multiple inputs.

Key difference from SISO: with m inputs, there are infinitely many K matrices that place the same poles. The extra freedom can be used to improve conditioning, robustness, or shape the eigenvectors.

### Method A: Direct Eigenstructure Assignment (Recommended)

This is the most practical approach for MIMO systems with Eigen. It places poles AND allows partial eigenvector specification. Even if you only care about poles (not eigenvectors), this method works and is numerically well-conditioned.

**Algorithm (Kautsky-Nichols-Van Dooren, 1985):**

Given plant (A, B) with n states and m inputs, desired eigenvalues {lambda_1, ..., lambda_n}:

For each desired eigenvalue lambda_i:
1. Compute the null space of [lambda_i * I - A, B]:
   ```
   N_i = null([lambda_i * I - A | B])    // dimension >= m
   ```
   This gives vectors [v_i; w_i] such that (lambda_i * I - A) * v_i = B * w_i.

2. The columns of N_i represent the achievable eigenvector/gain combinations. Any v_i from this space is a valid eigenvector for eigenvalue lambda_i.

3. If no partial eigenvector is specified, choose v_i to be as "robust" as possible (maximize the minimum angle between v_i and all other eigenvectors).

4. From each [v_i; w_i], extract the eigenvector v_i (first n components) and gain contribution w_i (last m components).

5. Assemble:
   ```
   V = [v_1, v_2, ..., v_n]    // n x n matrix of eigenvectors
   W = [w_1, w_2, ..., w_n]    // m x n matrix of gain contributions
   K = W * V^{-1}              // m x n gain matrix
   ```

**Implementation with Eigen:**

```cpp
PolePlacementResult mimoPlace(
    const LinearSystem& plant,
    const std::vector<std::complex<double>>& desired_poles)
{
    PolePlacementResult result;
    int n = plant.states();
    int m = plant.inputs();

    if (static_cast<int>(desired_poles.size()) != n) {
        result.success = false;
        result.error = "Number of desired poles must equal number of states";
        return result;
    }

    auto ctrb = checkControllability(plant);
    if (!ctrb.pass) {
        result.success = false;
        result.error = "System is not controllable";
        return result;
    }

    // Process poles, handling complex conjugate pairs together
    Eigen::MatrixXcd V(n, n);  // eigenvectors (complex)
    Eigen::MatrixXcd W(m, n);  // gain contributions (complex)

    int col = 0;
    while (col < n) {
        std::complex<double> lam = desired_poles[col];

        if (std::abs(lam.imag()) < 1e-12) {
            // Real pole: work with real arithmetic
            Eigen::MatrixXd M(n, n + m);
            M.leftCols(n) = lam.real() * Eigen::MatrixXd::Identity(n, n) - plant.A;
            M.rightCols(m) = plant.B;

            // Null space via SVD
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(M, Eigen::ComputeFullV);
            // Last column(s) of V_svd span the null space
            Eigen::VectorXd null_vec = svd.matrixV().col(n + m - 1);

            V.col(col) = null_vec.head(n).cast<std::complex<double>>();
            W.col(col) = null_vec.tail(m).cast<std::complex<double>>();
            col += 1;
        } else {
            // Complex pole: process with conjugate pair
            Eigen::MatrixXcd M(n, n + m);
            M.leftCols(n) = lam * Eigen::MatrixXcd::Identity(n, n)
                           - plant.A.cast<std::complex<double>>();
            M.rightCols(m) = plant.B.cast<std::complex<double>>();

            Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M, Eigen::ComputeFullV);
            Eigen::VectorXcd null_vec = svd.matrixV().col(n + m - 1);

            V.col(col) = null_vec.head(n);
            W.col(col) = null_vec.tail(m);

            // Conjugate for the paired pole
            V.col(col + 1) = null_vec.head(n).conjugate();
            W.col(col + 1) = null_vec.tail(m).conjugate();
            col += 2;
        }
    }

    // K = W * V^{-1}
    Eigen::MatrixXcd K_complex = W * V.inverse();

    // Result should be real (imaginary parts ~0)
    result.K = K_complex.real();
    result.success = true;

    // Verify: check that max imaginary part is small
    double max_imag = K_complex.imag().cwiseAbs().maxCoeff();
    if (max_imag > 1e-6) {
        result.error = "Warning: K has significant imaginary part ("
                     + std::to_string(max_imag) + "). "
                     "Check that poles are conjugate-paired.";
    }

    return result;
}
```

**Handling the null space dimension:**
- For a single-input system, the null space of [lambda*I - A | B] is 1-dimensional (no eigenvector choice).
- For m inputs, the null space is m-dimensional, giving freedom to choose among m possible eigenvector directions.
- The simple implementation above picks the last singular vector. A more sophisticated version (Kautsky-Nichols iteration) would iterate to maximize robustness.

### Method B: Bass-Gura Formula

Extension of Ackermann to canonical form. Not recommended for implementation because:
- Requires transformation to controller canonical form (numerically poor)
- Same numerical issues as Ackermann for high-order systems
- Does not leverage multi-input freedom

### Method C: Sylvester Equation Approach

Place poles by solving the Sylvester equation A*X - X*Lambda = B*W where Lambda = diag(desired poles). This is mathematically equivalent to the eigenstructure approach but uses a different computational path. Eigen does not have a built-in Sylvester solver, so this would require implementing one (Bartels-Stewart algorithm). Not worth the added complexity.

### Recommendation for This Project

Use **Method A** (direct eigenstructure assignment) for all cases:
- It handles both SISO and MIMO
- It naturally extends to eigenvector assignment (Section 3)
- It uses only SVD, which Eigen provides natively
- For SISO, fall back to Ackermann if desired (simpler code path, but not strictly necessary)

---

## 3. Eigenstructure Assignment

### Concept

Standard pole placement specifies WHERE the poles go (eigenvalues). Eigenstructure assignment additionally specifies HOW the modes couple to the states (eigenvectors).

**Why it matters:**
- Two systems with identical poles can have very different transient behaviors
- Eigenvectors determine which states participate in which mode
- In aircraft flight control: you might want the "dutch roll" mode to affect yaw but not roll
- In robotics: you might want decoupled joint responses

### User Inputs

For each desired pole lambda_i, the user optionally provides:
- **Desired eigenvector direction** d_i (partial specification): an n-vector indicating the preferred eigenvector direction. The algorithm finds the achievable eigenvector closest to this direction.
- If no direction is provided, the algorithm chooses automatically (same as Method A above).

### Algorithm (Extension of Method A)

For each desired eigenvalue lambda_i with desired direction d_i:

1. Compute null space N_i of [lambda_i * I - A | B] as before.

2. Extract the eigenvector portion: let N_v = first n rows of N_i (the "achievable eigenvector subspace").

3. Project the desired direction d_i onto the achievable subspace:
   ```
   v_i = N_v * (N_v^H * N_v)^{-1} * N_v^H * d_i
   ```
   This gives the closest achievable eigenvector to d_i.

4. Find the corresponding gain contribution w_i by solving:
   ```
   [v_i; w_i] = N_i * alpha_i
   ```
   where alpha_i is chosen so the eigenvector part matches v_i.

5. Assemble K = W * V^{-1} as before.

### Implementation Sketch

```cpp
struct EigenstructureSpec {
    std::complex<double> eigenvalue;
    std::optional<Eigen::VectorXcd> desired_direction;  // partial eigenvector
};

PolePlacementResult eigenstructureAssign(
    const LinearSystem& plant,
    const std::vector<EigenstructureSpec>& specs)
{
    int n = plant.states();
    int m = plant.inputs();
    PolePlacementResult result;

    Eigen::MatrixXcd V(n, n);
    Eigen::MatrixXcd W(m, n);

    int col = 0;
    for (const auto& spec : specs) {
        std::complex<double> lam = spec.eigenvalue;

        // Build [lambda*I - A | B]
        Eigen::MatrixXcd M(n, n + m);
        M.leftCols(n) = lam * Eigen::MatrixXcd::Identity(n, n)
                       - plant.A.cast<std::complex<double>>();
        M.rightCols(m) = plant.B.cast<std::complex<double>>();

        // Full SVD to get null space
        Eigen::JacobiSVD<Eigen::MatrixXcd> svd(
            M, Eigen::ComputeFullV);

        // Null space columns: those corresponding to zero singular values
        int null_dim = (n + m) - n;  // = m (for full-rank [lambda*I - A | B])
        Eigen::MatrixXcd null_space = svd.matrixV().rightCols(null_dim);

        Eigen::MatrixXcd N_v = null_space.topRows(n);    // achievable eigenvectors
        Eigen::MatrixXcd N_w = null_space.bottomRows(m); // corresponding gains

        Eigen::VectorXcd alpha;
        if (spec.desired_direction.has_value()) {
            // Project desired direction onto achievable subspace
            // Minimize || N_v * alpha - d ||
            alpha = N_v.colPivHouseholderQr().solve(
                spec.desired_direction.value());
        } else {
            // No preference: pick first null vector (simple choice)
            // Better: pick to maximize angle separation (Kautsky iteration)
            alpha = Eigen::VectorXcd::Zero(null_dim);
            alpha(0) = 1.0;
        }

        V.col(col) = N_v * alpha;
        W.col(col) = N_w * alpha;

        // Handle conjugate pair for complex eigenvalue
        if (std::abs(lam.imag()) > 1e-12 && col + 1 < n) {
            V.col(col + 1) = V.col(col).conjugate();
            W.col(col + 1) = W.col(col).conjugate();
            col += 2;
        } else {
            col += 1;
        }
    }

    Eigen::MatrixXcd K_complex = W * V.inverse();
    result.K = K_complex.real();
    result.success = true;
    return result;
}
```

### Practical Considerations for the UI

1. **Is eigenstructure assignment practical for a real-time UI tool?** Yes. The computation is dominated by m SVD decompositions of n x (n+m) matrices, plus one n x n matrix inversion. For n <= 8 (all current presets), this completes in microseconds.

2. **UI for specifying eigenvectors:** This is the hard part. Options:
   - **Simple mode (eigenvalue-only):** User just enters desired poles. Algorithm picks eigenvectors automatically. This covers 90% of use cases.
   - **Advanced mode:** For each pole, show a direction vector input. Default to the open-loop eigenvector as a starting point. Let the user modify individual components.
   - **Geometric mode:** On the time response plot, show which states participate in each mode. Let the user drag to adjust coupling. This is aspirational.

3. **Recommendation:** Start with eigenvalue-only mode (equivalent to basic pole placement). Add eigenvector specification as a future enhancement.

---

## 4. State Observer (Luenberger Observer)

### Why It Is Needed

State feedback u = -K*x requires measuring all states x. In practice, we only measure y = C*x (the outputs). An observer reconstructs the full state from the output measurements.

### Observer Design

The Luenberger observer is the dual of pole placement:

```
x_hat_dot = A * x_hat + B * u + L * (y - C * x_hat)
           = (A - L*C) * x_hat + B * u + L * y
```

where L is the observer gain matrix (n x p).

**The observer pole placement problem:** Find L such that eigenvalues of (A - L*C) are at desired observer pole locations.

**Duality:** This is identical to the state feedback problem for the dual system (A^T, C^T):
```
Find L^T such that eig(A^T - C^T * L^T) = desired observer poles
=> L = place(A^T, C^T, observer_poles)^T
```

So the same `mimoPlace()` function can be reused:

```cpp
PolePlacementResult designObserver(
    const LinearSystem& plant,
    const std::vector<std::complex<double>>& observer_poles)
{
    // Dual system: (A^T, C^T)
    LinearSystem dual;
    dual.A = plant.A.transpose();
    dual.B = plant.C.transpose();
    dual.C = Eigen::MatrixXd::Identity(plant.states(), plant.states());
    dual.D = Eigen::MatrixXd::Zero(plant.states(), plant.C.rows());

    auto result = mimoPlace(dual, observer_poles);
    if (result.success) {
        result.K = result.K.transpose();  // L = K_dual^T, now n x p
    }
    return result;
}
```

### Observer Pole Selection Rules

- Observer poles should be 2-10x faster than controller poles (further left in the s-plane)
- Common heuristic: multiply controller pole real parts by a factor of 3-5
- Too fast => high observer gains => noise amplification
- Too slow => state estimates lag behind reality

### Separation Principle

The controller and observer can be designed independently:
- Design K assuming full-state feedback
- Design L assuming known A, B, C

The combined system (plant + observer with state feedback) has eigenvalues that are the union of:
- eig(A - B*K) — controller poles
- eig(A - L*C) — observer poles

This is a fundamental result. The closed-loop eigenvalues are exactly the desired controller poles plus the desired observer poles.

### Observability Check

Before designing an observer, check `checkObservability(plant)`. If the system is not observable, the observer cannot reconstruct all states.

---

## 5. Closed-Loop Representations

### 5a. State Feedback Only (u = -K*x + r)

Already implemented as `stateFeedbackClose()`:

```
A_cl = A - B*K
B_cl = B
C_cl = C - D*K
D_cl = D
```

The existing `ControllerType::GainMatrix` path in `visualizer.cpp` (line 84-89) handles this.

### 5b. Observer-Based Output Feedback (Augmented System)

When using an observer, the combined plant + observer system has 2n states:

```
State vector: [x; x_hat]     (plant states; estimated states)

Plant:     x_dot = A*x + B*u
Observer:  x_hat_dot = (A - L*C)*x_hat + B*u + L*C*x
Control:   u = -K*x_hat + r

Substituting u into both:
  x_dot     = A*x - B*K*x_hat + B*r
  x_hat_dot = L*C*x + (A - B*K - L*C)*x_hat + B*r

Augmented state-space:
  A_aug = [  A        -B*K      ]
          [  L*C    A - B*K - L*C ]

  B_aug = [B]
          [B]

  C_aug = [C   -D*K]     (output from plant states)

  D_aug = D
```

**Alternative (error-based) formulation** using e = x - x_hat:

```
State vector: [x; e]

A_aug = [A - B*K    B*K   ]
        [   0     A - L*C  ]

B_aug = [B]
        [0]

C_aug = [C - D*K   D*K]

D_aug = D
```

The second form makes the separation principle visible: the eigenvalues of A_aug are eig(A - B*K) union eig(A - L*C), since A_aug is block-triangular.

### Implementation

```cpp
LinearSystem observerFeedbackClose(
    const LinearSystem& plant,
    const Eigen::MatrixXd& K,   // m x n state feedback gain
    const Eigen::MatrixXd& L)   // n x p observer gain
{
    int n = plant.states();
    int m = plant.inputs();
    int p = plant.outputs();

    LinearSystem result;

    // Augmented system: 2n states
    result.A = Eigen::MatrixXd::Zero(2 * n, 2 * n);
    result.A.topLeftCorner(n, n) = plant.A - plant.B * K;
    result.A.topRightCorner(n, n) = plant.B * K;
    result.A.bottomRightCorner(n, n) = plant.A - L * plant.C;
    // bottom-left is zero (error form)

    result.B = Eigen::MatrixXd::Zero(2 * n, m);
    result.B.topRows(n) = plant.B;
    // bottom rows are zero (error form: reference doesn't enter error dynamics)

    result.C = Eigen::MatrixXd::Zero(p, 2 * n);
    result.C.leftCols(n) = plant.C - plant.D * K;
    result.C.rightCols(n) = plant.D * K;

    result.D = plant.D;

    return result;
}
```

### 5c. How This Connects to the Existing Framework

The current `AppState` supports these controller types:

| `ControllerType` | Controller slot | Closed-loop formation |
|---|---|---|
| `None` | — | No closed-loop system |
| `StateSpace` | `ctrl_ss` (A,B,C,D) | `series(ctrl, plant)` then `feedbackConnect()` |
| `GainMatrix` | `ctrl_K` | `stateFeedbackClose(plant, K)` |

To add observer-based feedback, there are two approaches:

**Approach 1: New ControllerType (recommended)**

Add `ControllerType::ObserverFeedback` with fields for both K and L. The recompute section in `visualizer.cpp` would call `observerFeedbackClose(plant, K, L)`.

```cpp
enum class ControllerType { None, StateSpace, GainMatrix, ObserverFeedback };

// In AppState:
Eigen::MatrixXd observer_L;        // n x p observer gain
std::vector<std::complex<double>> desired_poles;     // for display
std::vector<std::complex<double>> observer_poles;     // for display
```

**Approach 2: Represent observer as StateSpace controller**

The observer + state feedback can be written as a dynamic output feedback controller:

```
Controller state-space:
  A_c = A - B*K - L*C + L*D*K
  B_c = L
  C_c = K
  D_c = 0
```

This is a proper dynamic controller that takes y as input and produces u as output. It could be entered as a `ControllerType::StateSpace` and connected via the existing `seriesConnect` + `feedbackConnect` path.

**Advantage of Approach 2:** No changes to the connection framework. The existing series/feedback pipeline handles it.

**Disadvantage:** The user loses visibility into K and L separately. The observer poles are hidden inside A_c.

**Recommendation:** Use Approach 1 for the UI (explicit K and L) but also provide a utility to convert to Approach 2 form for analysis.

---

## 6. Summary of Functions to Implement

### Core Algorithms (new file: `src/analysis/pole_placement.h/.cpp`)

| Function | Inputs | Output | Notes |
|---|---|---|---|
| `polyFromRoots()` | vector of complex poles | vector of real polynomial coefficients | Utility |
| `polyEvalMatrix()` | matrix A, polynomial coefficients | matrix alpha(A) | Utility |
| `ackermannPlace()` | plant, desired_poles | K (1 x n) | SISO only |
| `mimoPlace()` | plant, desired_poles | K (m x n) | General case, uses SVD |
| `eigenstructureAssign()` | plant, eigenvalue+direction specs | K (m x n) | Full eigenvector control |
| `designObserver()` | plant, observer_poles | L (n x p) | Uses duality + mimoPlace |

### System Formation (extend `src/analysis/system_connect.h/.cpp`)

| Function | Inputs | Output | Notes |
|---|---|---|---|
| `stateFeedbackClose()` | plant, K | LinearSystem | Already exists |
| `observerFeedbackClose()` | plant, K, L | LinearSystem (2n states) | New |
| `observerAsController()` | plant, K, L | LinearSystem (controller form) | New, utility |

### Verification Utilities

| Function | Purpose |
|---|---|
| `verifyPlacement()` | Compute eig(A - B*K), compare against desired poles |
| `placementError()` | Max distance between achieved and desired pole locations |

### UI (modify `src/panels/model_panel.cpp`)

- Desired pole entry: text fields or interactive pole-zero plot (click to place)
- Observer pole entry: similar, with multiplier heuristic
- Display: computed K matrix, computed L matrix, achieved poles, placement error

### AppState Additions

```cpp
// In AppState:
std::vector<std::complex<double>> desired_cl_poles;   // user-specified
std::vector<std::complex<double>> desired_obs_poles;   // user-specified
Eigen::MatrixXd observer_L;                            // observer gain
bool use_observer = false;
```

---

## 7. Eigen Library Capabilities Summary

Everything required is available in Eigen 3.4 without additional libraries:

| Need | Eigen Class/Function |
|---|---|
| Eigenvalue decomposition | `Eigen::EigenSolver<MatrixXd>` |
| Complex eigenvalues | `Eigen::EigenSolver::eigenvalues()` returns `VectorXcd` |
| SVD (for null space) | `Eigen::JacobiSVD<MatrixXd>` or `Eigen::BDCSVD<MatrixXd>` |
| Complex SVD | `Eigen::JacobiSVD<MatrixXcd>` |
| QR factorization | `ColPivHouseholderQR` (already used) |
| Matrix inverse | `MatrixXcd::inverse()` (fine for n <= 8) |
| Generalized eigenvalues | `Eigen::GeneralizedEigenSolver` (already used for zeros) |

No external control library (e.g., Slicot, control-toolbox) is needed.

---

## 8. Complexity Assessment

| Feature | Implementation Effort | Numerical Risk | UI Complexity |
|---|---|---|---|
| Ackermann (SISO) | Low | Low (n <= 8) | Low |
| MIMO eigenvalue placement | Medium | Low-Medium | Medium |
| Eigenstructure assignment | Medium | Low-Medium | High (eigenvector input) |
| Luenberger observer | Low (reuses placement) | Low | Medium |
| Observer-based closed-loop | Low | None | Medium |
| Interactive pole placement UI | High | None | High |

**Recommended build order:**
1. Core placement functions (Ackermann + MIMO place) + verification
2. Observer design (trivial given placement functions)
3. Closed-loop formation (observerFeedbackClose)
4. UI: pole entry as text fields, computed K display
5. UI: interactive pole drag on pole-zero plot
6. Eigenstructure assignment (advanced feature, later)
