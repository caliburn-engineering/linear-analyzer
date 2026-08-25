# Relative Gain Array at a Frequency as a Pairing Diagnostic

Research for the Caliburn Linear System Analyzer. Resolves issue #4 (part of map #1). All algorithms target C++17 with Eigen 3.4, no external control libraries.

**Primary sources.** The canonical treatment is Skogestad & Postlethwaite, *Multivariable Feedback Control: Analysis and Design*, 2nd edition (Wiley, 2005) — hereafter **S&P**. The RGA appears there in three places: §3.4 (introduction, pp. 82–91), §10.6 (decentralized control and the pairing rules, pp. 429–454), and Appendix A.4 (algebraic properties and the non-square extension, pp. 526–529). The author hosts the full text at [skoge.folk.ntnu.no/book/ps/bookall.pdf](https://skoge.folk.ntnu.no/book/ps/bookall.pdf) and the ch. 3 + ch. 10 extract at [skoge.folk.ntnu.no/vgprosessregulering/papers-pensum/ch3-ch10-2ndedition.pdf](https://skoge.folk.ntnu.no/vgprosessregulering/papers-pensum/ch3-ch10-2ndedition.pdf). Every S&P citation below was read from those PDFs and equation numbers are quoted as printed.

**Note on Bristol.** The original paper is E. H. Bristol, "On a new measure of interaction for multivariable process control", *IEEE Transactions on Automatic Control*, vol. 11, no. 1, pp. 133–134, Jan. 1966 ([DOI 10.1109/TAC.1966.1098266](https://doi.org/10.1109/TAC.1966.1098266)). It is a two-page correspondence behind the IEEE paywall and **was not read directly for this document**. Where Bristol's reasoning is used below it is taken from S&P §3.4.1, which opens "We follow Bristol (1966) here" and reproduces the derivation in full. Claims attributed to Bristol are therefore attributed *via* S&P, not verified against the 1966 text.

**Note on reference implementations.** There is no RGA function in `python-control` (checked `control/` source tree and the ReadTheDocs global index on `main` — no `relative_gain_array`, no `rga`), and no RGA function documented in MATLAB's Control System Toolbox or Robust Control Toolbox (checked the Robust Control Toolbox function index; the only MATLAB implementations found are File Exchange community submissions, which are not primary sources). **The de facto reference implementation is the MATLAB snippet printed by Skogestad himself in S&P Table 3.1, p. 87** — reproduced verbatim in §1.2 below. That snippet, not a library, is what an implementer should transcribe.

---

## 1. Definition of the RGA at a Frequency

### 1.1 The form is confirmed

S&P eq. (3.54), p. 82:

> "The RGA (Bristol, 1966) of a non-singular square **complex** matrix `G` is a square **complex** matrix defined as
> `RGA(G) = Λ(G) ≜ G × (G⁻¹)ᵀ`
> where `×` denotes element-by-element multiplication (the Hadamard or Schur product)."

Identically restated as eq. (A.77), p. 526. Elementwise (S&P eq. 3.61 / A.78):

```
λij = [G]ij · [G⁻¹]ji
```

Note the **index swap**: element `(i,j)` of Λ multiplies element `(i,j)` of `G` by element `(j,i)` of `G⁻¹`. That is what the transpose in `(G⁻¹)ᵀ` encodes. Getting this backwards produces `Λ = G ∘ G⁻¹` which is not the RGA and has none of its properties.

So the ticket's stated form `Λ(jω) = G(jω) ∘ (G(jω)⁻¹)ᵀ` is **correct as written**, with one caveat that matters enormously in code — see §1.3.

For the 2×2 case S&P give the closed form (eq. 3.55):

```
Λ = [ λ11      1 - λ11 ]                            1
    [ 1 - λ11  λ11     ]        λ11 = ---------------------------
                                       1 - (g12·g21)/(g11·g22)
```

Useful as a unit test but **not** as the implementation — the general form is needed for p×p.

### 1.2 Λ is formed from the complex matrix directly, not from magnitudes

This is the ticket's sharpest question and the literature is unambiguous.

**The definition itself says "complex".** S&P eq. (3.54) and eq. (A.77) both define Λ for a *complex* matrix and both state the result is a *complex* matrix. A.4.1 adds: "If `A` is real then `Λ(A)` is also real" — the converse implication being that when `A` is complex, Λ generally is not.

**Skogestad's own frequency-domain code confirms it.** S&P Table 3.1, p. 87, "Matlab program to calculate frequency-dependent RGA":

```matlab
omega = logspace(-5,2,61);
for i = 1:length(omega)
    Gf = freqresp(G,omega(i));                      % G(jω)
    RGAw(:,:,i) = Gf.*inv(Gf).';                    % RGA at frequency omega
    RGAno(i)    = sum(sum(abs(RGAw(:,:,i) - eye(2))));  % RGA number
end
RGA = frd(RGAw,omega);
```

`Gf` is the complex frequency response. `inv(Gf)` is the complex inverse. `.'` is applied to the complex inverse. **No magnitude is taken anywhere until the RGA number line**, where `abs()` is applied to `Λ − I`, i.e. to the already-complex Λ.

**The worked example is explicitly complex.** S&P Example 3.11, eq. (3.66), p. 86 — a 2×2 pressurized-vessel model evaluated at ω = 0.01 rad/s:

```
Λ = [ 0.2469 + 0.0193i    0.7531 - 0.0193i ]
    [ 0.7531 - 0.0193i    0.2469 + 0.0193i ]
```

A complex Λ printed with complex entries, and the rows and columns sum to exactly 1 as required by property A2 (§4.1 below). There is no ambiguity: **build `G(jω)` as complex, invert as complex, Hadamard as complex, then take magnitudes only for display and for the RGA number.**

**Doing it the other way round is not an approximation — it changes the answer's sign.** Reproducing eq. (3.65) with the coefficients as printed in the 2nd edition and evaluating at ω = 0.01 rad/s gives:

| method | λ₁₁ |
|---|---|
| `Λ(G(jω))` — correct, complex throughout | `+0.2416 + 0.0190i` |
| S&P's printed value, eq. (3.66) | `+0.2469 + 0.0193i` |
| `Λ(\|G(jω)\|)` — magnitudes taken first | **`−0.4696`** |
| `G ∘ (G⁻¹)ᴴ` — conjugate transpose instead of transpose | **`−0.2408 − 0.0280i`** |

(The ~2% gap between the first two rows is a coefficient-rounding difference in the printed model, not a methodological one — both agree on sign, magnitude and the row/column sums. It does not affect any conclusion. Verified numerically in pure-Python complex arithmetic; the discrepancy was not resolved against the 1st edition and is flagged rather than explained away.)

The magnitudes-first result is **negative**, which under pairing rule 2 (§4.2) would fire a "wrong pairing / integral instability" alarm on a plant whose correct λ₁₁ is a healthy positive 0.24. Magnitudes-first is not a cheap approximation of the RGA. It is a different, wrong quantity.

### 1.3 Transpose, **not** conjugate transpose — the single most likely implementation bug

S&P flag this explicitly in footnote 5 on p. 83:

> "The symbol `'` in Matlab gives the conjugate transpose (`Aᴴ`), and we must use `.'` to get the 'regular' transpose (`Aᵀ`)."

The definition requires the plain transpose. With Eigen:

| operation | Eigen | correct? |
|---|---|---|
| plain transpose | `.transpose()` | **yes** |
| conjugate transpose | `.adjoint()` | **no** |

`Eigen::MatrixXcd::adjoint()` is the conjugate transpose and is the natural-looking thing to reach for on a complex matrix. It is wrong here. As the table in §1.2 shows, on a genuinely complex `G` the two differ in sign, so this bug will not be caught by a real-matrix (DC) unit test — it only surfaces at ω > 0. **Test the RGA at a frequency where `G(jω)` has non-trivial phase, or the bug ships.**

A cheap runtime detector: for a real `G`, `.transpose()` and `.adjoint()` agree, so a DC test proves nothing; but for any `G`, the correct Λ has rows and columns summing to exactly 1 (property A2), and the `.adjoint()` version generally does not. See the self-check in §6.3.

---

## 2. Which Frequency to Evaluate At

### 2.1 What the literature says: the closed-loop bandwidth / crossover region

S&P state this as a named rule, twice, in identical terms.

**§3.4.1, p. 85 (Pairing rule 1, forward reference to p. 450):**

> "Prefer pairings such that the rearranged system, with the selected pairings along the diagonal, has an RGA matrix close to identity **at frequencies around the closed-loop bandwidth**."

**§10.6.9, p. 450 (the rule in its home section, titled "Pairing rule 1. RGA at crossover frequencies"):**

> "Prefer pairings such that the rearranged system, with the selected pairings along the diagonal, has an RGA matrix close to identity at frequencies around the closed-loop bandwidth."

Note the title says *crossover*, the body says *closed-loop bandwidth*. S&P use the two interchangeably here; for a well-damped loop they coincide to within a factor of ~1–2. The book does not distinguish them for this purpose.

The rationale, S&P §10.6.9 p. 450:

> "Pairing rule 1 is to ensure that we have diagonal dominance where interactions from other loops do not cause instability."

### 2.2 The literature actively warns against DC-only

S&P are unusually pointed about this. §3.4.1 Remark, p. 84:

> "The assumption of `yk = 0` ('perfect control of `yk`') in (3.57) is satisfied at steady-state (ω = 0) provided we have integral action in the loop, but it will generally not hold exactly at other frequencies. Unfortunately, this has led many authors to dismiss the RGA as being 'only useful at steady-state' or 'only useful if we use integral action'. On the contrary, **in most cases it is the value of the RGA at frequencies close to crossover which is most important**, and both the gain and the phase of the RGA elements are important."

And §10.6.9 Remark, p. 450:

> "Even if we have `λii(0) = 1` and `λii(∞) = 1` for all `i`, this does not necessarily mean that the diagonal pairing is the best, even for a 2×2 plant. The reason for this is that the behaviour at 'intermediate' bandwidth frequencies is more important."

That remark points at Example 3.11, which is the canonical demonstration: a plant with `Λ(0) = I` and `Λ(j∞) = I` — both endpoints saying "pair on the diagonal" — for which the correct answer at the actual bandwidth is the *reverse* pairing (S&P p. 86). Both DC and high-frequency evaluation give the wrong answer on that plant.

**Conclusion for the ticket: DC-with-fallback is not merely awkward for integrator plants, it is the evaluation point the canonical source specifically warns produces wrong pairings.** Retiring the DC heuristic is right for a second, independent reason beyond the Ball-Balancer's infinite DC gain.

### 2.3 The one thing that *is* anchored at DC

Pairing rule 2 is stated at steady-state and only at steady-state. S&P §10.6.9, p. 450:

> "**Pairing rule 2.** For a stable plant avoid pairings that correspond to negative steady-state RGA elements, `λij(0) < 0`."

This is not a stylistic choice. It descends from Theorem 10.6 (p. 443), the Grosdidier–Holt–Morari DIC result, which is a statement about integral action and therefore about `s = 0`. A negative `λii` at ω = 500 rad/s does not imply what a negative `λii(0)` implies.

**Design consequence:** a single user-chosen ω cannot serve both rules. The honest UI shows rule 1 at the chosen ω and, *separately and only when `G(0)` is finite*, rule 2 at DC. For plants with poles at the origin, rule 2 is simply unavailable — see §5.4.

### 2.4 Fixed default vs. user-chosen ω — trade-offs

| option | for | against |
|---|---|---|
| **Fixed ω = 0 (DC)** | one number, no user input; matches rule 2 exactly | undefined for integrator plants; §2.2 shows it gives wrong pairings on plants whose Λ varies with frequency; contradicts the source's own advice |
| **Plant gain crossover** | automatic; no user input; "crossover" is the word in rule 1's title | **ill-defined for a MIMO plant.** `computeBode` in `frequency_response.h` is per-channel — there are p·m crossovers, not one. Well-defined only via `σ̄(G(jω)) = 1`, which requires an SVD sweep and is scaling-dependent (a crossover in metres differs from one in millimetres). Also circular: crossover of the *open-loop plant* is not the *closed-loop bandwidth* rule 1 asks for — that depends on the controller you have not designed yet |
| **User-chosen ω, defaulted** | matches the rule as stated ("the intended closed-loop bandwidth" is a design decision, not a plant property); one slider; cheap to recompute | requires the user to have an intent; a bad default is silently wrong |
| **Sweep Λ over the Bode grid and plot it** | this is what S&P actually do — Figure 3.8(a) plots `\|λij\|` vs. frequency and 3.8(b) plots the RGA number vs. frequency; reveals frequency-dependent pairing flips, which a single ω hides by construction | more UI than a number; needs the phase warning of §4.4 attached to any magnitude plot |

**Recommendation.** Make ω a user input, defaulted, and *additionally* plot the RGA number over the existing Bode frequency grid. The plot is nearly free — the app already sweeps a log grid and already evaluates `evalTransferFunction` per point — and it is the presentation the canonical source uses. The single-ω readout is then a cursor on that plot rather than the whole diagnostic.

For the default value, the defensible choices in descending order of honesty:

1. **The dominant plant dynamics.** `ω_default = |λ_dominant|` where `λ_dominant` is the finite plant eigenvalue of largest magnitude, or the geometric mean of the finite eigenvalue magnitudes. Cheap (`Eigen::EigenSolver` is already used for pole-zero), never zero, and defensible as "the region where the plant is actually doing something."
2. **The geometric centre of the existing Bode grid**, `sqrt(f_min · f_max)`. Arbitrary but transparent and never degenerate.
3. **Not DC.** Whatever is chosen, `ω = 0` must not be the default, per §2.2.

I could not find a primary source prescribing an automatic default ω for RGA evaluation. S&P always either plot the sweep or state the bandwidth as a design assumption ("if we use decentralized control and the closed-loop bandwidth is around 0.01 rad/s", p. 86). **Treat any automatic default as an app convenience, not a literature-backed value, and label it as such in the UI.**

### 2.5 Units gotcha specific to this codebase

`FrequencyResponse` and `computeBode` work in **Hz** (`freq_hz`, `gain_crossover_hz`) and convert internally via `omega = 2·π·freq_hz`. Every S&P frequency quoted in this document is in **rad/s**. Whichever the RGA panel exposes, label it, and convert once at the boundary — `s = std::complex<double>(0.0, 2.0*M_PI*freq_hz)` to match the existing convention in `frequency_response.cpp:41-43`.

---

## 3. Non-Square and Subset Handling

The map's model pairs a subset of channels. Two distinct generalizations exist and they answer **different questions**. Conflating them is the main risk here.

### 3.1 The square sub-matrix approach is sound — and it is the right one for pairing

Take the loop list, and build `G_S(jω)` as the p_loops × p_loops matrix whose row `k` is loop `k`'s output index and whose column `k` is loop `k`'s input index. Compute `Λ = G_S ∘ (G_S⁻¹)ᵀ` on that. **This is correct**, and it follows directly from Bristol's derivation as S&P restate it (§3.4.1, pp. 83–84).

Bristol's `λij` is the ratio of two gains from `uj` to `yi`:

- **all other loops open** — S&P eq. (3.56): `(∂yi/∂uj)` with `uk = 0 ∀k ≠ j`, which equals `gij`.
- **all other loops closed with perfect control** — S&P eq. (3.57)–(3.58): `(∂yi/∂uj)` with `yk = 0 ∀k ≠ i`, which equals `1/[G⁻¹]ji`.

"All other loops" means the loops of *the decentralized scheme you are designing*. In the subset case:

- Inputs not in the loop list are **inert** — the map's standing decision is that unpaired reference channels are a zero column in `C`. They are never manipulated, so they are held at zero in both extreme cases. The "other loops open" gain is therefore still `gij`.
- Outputs not in the loop list are **not controlled**, so the "perfect control" condition `yk = 0` is imposed only over `k ∈ S`, and the inversion that realizes it is the inversion of `G_S`.

Both of Bristol's two extreme cases therefore land exactly on the square sub-matrix. **Λ(G_S) is Bristol's relative gain for the control structure actually being built.** This is a derivation from the primary definition, not a quoted theorem; S&P do not state the subset case explicitly. It is straightforward, but it is mine, not theirs.

### 3.2 Two caveats on the sub-matrix, one of them severe

**Caveat A — `Λ(G_S)` is not a sub-matrix of `Λ(G)`.** Every `λij` depends on the *whole* matrix through `det`: S&P eq. (A.78), p. 527, gives `λij = aij·det(Aij)/det(A)·(−1)^{i+j}`. Deleting rows and columns changes both determinants. Do **not** compute Λ on the full plant and then read off the paired entries — that is a different and meaningless number. Build `G_S` first, then invert.

**Caveat B — for a single loop, Λ is identically 1 and the diagnostic is vacuous.** If the loop list has one entry, `G_S = [g]` is 1×1, `G_S⁻¹ = [1/g]`, and `Λ = [g · (1/g)] = [1]` for every plant, every ω, always. This is not a "good pairing" signal; it is no signal at all.

**This directly contradicts an assumption in the ticket and the map.** The map records that three of the five presets — Inverted Pendulum, Quarter-Car and Double Mass-Spring — are 1-input / 2-output. A 1-input plant admits **at most one loop**, so the square sub-matrix is always 1×1 and the RGA pairing diagnostic says `1.0` regardless of which output you paired to. For 3 of 5 presets the sub-matrix RGA cannot distinguish a good pairing from a bad one. The pairing diagnostic is only informative when **two or more loops are paired**, which in the current preset library means the Ball-Balancer alone.

The UI must therefore say "not applicable — a single loop has no interaction to measure" for the 1-loop case rather than displaying a reassuring `Λ = 1.00`. A green `1.00` on a plant where the user paired the wrong output would be actively misleading.

### 3.3 The non-square RGA answers the *other* question — and it is the useful one for SIMO/MISO

S&P Appendix A.4.2, eq. (A.84), p. 528:

> "The RGA may be generalized to a non-square `l × m` matrix `A` by use of the pseudo-inverse `A†` defined in (A.62). We have `Λ(A) = A × (A†)ᵀ`."

The ticket's stated form `Λ = G ∘ (G†)ᵀ` is confirmed. But note **what S&P say it is for** — control property C3, §3.4.5, p. 89:

> "**Non-square plants.** The definition of the RGA may be generalized to non-square matrices by using the pseudo-inverse; see Appendix A.4.2. **Extra inputs:** If the sum of the elements in a column of RGA is small (≪ 1), then one may consider deleting the corresponding input. **Extra outputs:** If all elements in a row of RGA are small (≪ 1), then the corresponding output cannot be controlled."

The non-square RGA is a **channel-selection screen**, not a pairing measure. It tells you *which inputs are worth keeping* and *which outputs are reachable at all*. That is precisely the question the SIMO presets pose and precisely the question the 1×1 sub-matrix cannot answer.

**Concretely, for a p×1 plant** (one input, p outputs — the Quarter-Car / Inverted Pendulum / Double Mass-Spring case), the pseudo-inverse form collapses to something trivially implementable. With `g = G(jω)` a column vector, `G† = gᴴ/‖g‖²`, so `(G†)ᵀ = conj(g)/‖g‖²` and

```
λi = gi · conj(gi) / ‖g‖²  =  |gi|² / Σk |gk|²
```

which is **real, non-negative, and sums to exactly 1** — consistent with S&P property A.4.2-2(b) for a full-column-rank matrix ("the elements in each column of the RGA sum to 1"). It is the fraction of the input's total squared effect that lands on output `i`. For the Quarter-Car it answers "can this actuator move sprung mass position at all, or is its effect essentially all on the unsprung mass?", which is the SIMO analogue of a pairing decision. This is worth surfacing.

Two properties survive the generalization and the rest do not. S&P A.4.2, p. 528:

> "Properties 1 (transpose and inverse) and 2 (permutations) of the RGA also hold for non-square matrices, **but the remaining properties do not apply in the general case.**"

Full row rank keeps row sums = 1 and output-scaling invariance; full column rank keeps column sums = 1 and input-scaling invariance (A.4.2 items 1 and 2). In general only `Σᵢⱼ λij = rank(G)` survives (eq. A.87). **In particular, the two-sided scaling invariance that makes the square RGA unit-free is lost.**

### 3.4 The pseudo-inverse caveat that S&P do not mention

There is a live objection in the recent literature to using the Moore–Penrose pseudo-inverse for this at all.

Uhlmann, "On the Relative Gain Array (RGA) with Singular and Rectangular Matrices", arXiv:1805.10312v3 (2019), §III, proves that the MP inverse does **not** satisfy diagonal-scaling consistency:

> "`(DAE)⁻ᴾ ≠ E⁻¹ · A⁻ᴾ · D⁻¹`"

whereas the true inverse does, and shows the concrete consequence (his eqs. 19–22): for `G = ones(3,3)`, both the MP-RGA and his alternative give `(1/9)·ones(3,3)`; but after scaling the first row and column by 2, the MP-RGA changes to `(1/9)·[[4,1,1],[1,¼,¼],[1,¼,¼]]` while a unit-consistent generalized inverse leaves it unchanged. **The MP-RGA of a rectangular plant depends on your choice of units.**

The follow-up — Al Yousuf & Uhlmann, "On Use of the Moore-Penrose Pseudoinverse for Evaluating the RGA of Non-Square Systems", arXiv:2106.09766v2 (2022) — shows this is not academic. On a real non-square crude-distillation plant they demonstrate that changing the temperature unit by a factor of 10 changes the *set of pairings the MP-RGA recommends* (§3, eqs. 19–22). Their §4:

> "we have shown in a practical example that a simple decimation of the unit of temperature is sufficient to change the set of pairings indicated by the MP-RGA. This tends to strongly undermine confidence in the MP-RGA as a controller design tool when applied with non-square systems."

They also caution against over-reading their own alternative: "our results should not be interpreted as providing general evidence of the practical efficacy of the UC-RGA … future work is needed."

**Practical resolution for this project.** Do not implement the UC inverse — it is recent, not in Eigen, and its practical merit is explicitly unsettled by its own authors. Instead:

1. **Prefer the square sub-matrix wherever there are ≥ 2 loops.** The square RGA is exactly scale-invariant (S&P A.4.1 property 5: `Λ(D₁AD₂) = Λ(A)`) and none of this applies to it.
2. **When showing a non-square RGA, scale first and say so.** S&P §1.4, p. 5, prescribes the scaling: "dividing each variable by its maximum expected or allowed change", and S&P Remark 1, p. 6, warns "A number of the interpretations used in the book depend critically on a correct scaling." Since the p×1 formula above is literally `|gi|²/Σ|gk|²`, a plant whose outputs are in metres and radians will give an answer dominated by whichever unit happens to be numerically larger. Either require per-output scale factors or label the readout as unit-dependent. **Silently pinv-ing an unscaled plant is the failure mode both papers describe.**

---

## 4. Interpretation and UI Thresholds

### 4.1 The algebraic facts the UI can rely on

From S&P §3.4.4, p. 88 and Appendix A.4.1, p. 527:

| # | property | UI use |
|---|---|---|
| A1 / 5 | Scale-invariant: `Λ(D₁GD₂) = Λ(G)` for diagonal `D` | no units needed for the square RGA; also enables the equilibration guard in §5.3 |
| A2 / 3 | Every row and every column sums to 1 | **runtime self-check** — deviation from 1 measures how badly the inverse failed |
| A3 / 4 | `Λ = I` iff `G` is triangular (and for diagonal `G`) | `Λ − I` is a measure of *two-way* interaction; a triangular plant has one-way interaction and still gives `Λ = I` |
| 2 | `Λ(P₁GP₂) = P₁Λ(G)P₂` for permutations | **Λ is computed once, not per candidate pairing** — reordering the loops permutes Λ identically |
| A5 / 7 | Large RGA elements ⟹ large minimized condition number `γ*` (eq. A.80); converse false | large Λ implies ill-conditioning, but a large condition number alone does not imply large Λ (S&P Example 3.12: `diag(100,1)` has `γ = 100` but `Λ = I`) |

Property 2 is worth exploiting: when the user re-pairs, you do not need to re-invert. Only the comparison target changes.

### 4.2 What the values mean

Recall `λij = "open-loop gain (other loops open)" / "closed-loop gain (other loops closed)"` (S&P eq. 3.61 and the boxed derivation on p. 85). Everything below follows from that ratio.

| `λii` | meaning | source |
|---|---|---|
| **≈ 1** | Closing the other loops does not change this loop's gain. The loop can be tuned as if it were SISO. This is the target. | S&P p. 84: "we prefer to pair variables `uj` and `yi` so that `λij` is close to 1 at all frequencies, because this means that the gain from `uj` to `yi` is unaffected by closing the other loops" |
| **0 < λ < 1** | Closing the other loops *increases* the gain (`ĝ > g`). Interaction helps; loop will be more aggressive than its SISO tuning suggests. | direct from the ratio definition |
| **≈ 0.5** | Gain doubles when the other loops close. For a 2×2 plant, `λ11 = 0.5` forces `λ12 = 0.5` by property A2 — **the RGA cannot discriminate between the two pairings.** Severe two-way interaction and an ambiguous answer. | ratio definition + S&P property A2. The blending example (S&P Example 3.9, `Λ = [[0.2,0.8],[0.8,0.2]]`) is the near-miss case where 0.8 vs 0.2 still decides it |
| **λ > 1** | Closing the other loops *reduces* the gain — interaction fights the loop. Mild for λ ≈ 1.5–2. | ratio definition |
| **λ = 0** | Either `gij = 0` (the paired channel has no direct effect) or the cofactor vanishes. Loop `i` by itself does nothing; you would be relying entirely on interaction. | S&P §10.6.10, p. 451: "with sequential design one may choose to pair on an element with `gii = 0` (and `λii = 0`), which violates both pairing rules 1 and 3" |
| **λ < 0** | **The sign of the gain flips when the other loops close.** With integral action this is an integrity failure. | S&P Theorem 10.6, p. 443 — see below |
| **λ ≫ 1 (5–10 or more)** | Fundamentally difficult plant: strong interaction plus sensitivity to input uncertainty. Not just a bad pairing — a bad plant. | S&P control property C1, p. 89 |

**Negative λ — the loud case.** S&P Theorem 10.6, p. 443 (Grosdidier, Holt & Morari, 1985):

> "Consider a stable square plant `G` and a diagonal controller `K` with integral action in all elements, and assume that the loop transfer function `GK` is strictly proper. If a pairing of outputs and manipulated inputs corresponds to a negative steady-state relative gain, then the closed-loop system has at least one of the following properties:
> (a) The overall closed-loop system is unstable.
> (b) The loop with the negative relative gain is unstable by itself.
> (c) The closed-loop system is unstable if the loop with the negative relative gain is opened (broken)."

summarized as (eq. 10.79): "A stable (reordered) plant `G(s)` is DIC **only if** `λii(0) ≥ 0` for all `i`." For 2×2 this tightens to an equivalence (eq. 10.82): `DIC ⇔ λ11(0) > 0`.

Every one of (a), (b), (c) is a real failure. S&P p. 443: "situation (c) is also highly undesirable as it will imply instability if the loop with the negative relative gain somehow becomes inactive, e.g. due to input saturation." **This is the loudest warning the diagnostic can give — but note the theorem's hypotheses: stable plant, integral action, strictly proper `GK`, and steady state.** See §5.4.

**Large λ.** S&P C1, p. 89:

> "Large RGA elements (typically, **5 − 10 or larger**) at frequencies important for control indicate that the plant is fundamentally difficult to control due to strong interactions and sensitivity to uncertainty."

and specifically: "decouplers or other inverse-based controllers should not be used for plants with large RGA elements". The "5–10" band is the only explicit numeric threshold S&P give.

### 4.3 The scalar to show: the RGA number

S&P §3.4.3, eq. (3.67), p. 87:

> "**RGA number.** A simple measure for selecting pairings according to rule 1 is to prefer pairings with a small RGA number. For a diagonal pairing, `RGA number ≜ ‖Λ(G) − I‖sum` … where we have (somewhat arbitrarily) chosen the sum norm, `‖A‖sum = Σi,j |aij|`."

and for a non-diagonal pairing (p. 88):

> "The RGA number for other pairings is obtained by subtracting 1 for the selected pairings; for example, `‖Λ(G) − [[0,1],[1,0]]‖sum` for the off-diagonal pairing for a 2×2 plant."

So in general: **`RGAnum = ‖Λ − Π‖sum`**, where `Π` is the permutation matrix of the candidate pairing. If `G_S` is built with the paired channels already on the diagonal (§6.1), `Π = I` and the formula is just `(Λ - MatrixXcd::Identity()).cwiseAbs().sum()`.

The RGA number is the right single number for the UI because `|λ − 1|` correctly penalizes `λ = −1` (distance 2) and `λ = 5` (distance 4) and `λ = 0` (distance 1), all of which a magnitude readout mishandles. It also collapses a p×p matrix to one scalar that can be plotted against frequency — which is exactly S&P Figure 3.8(b).

One documented limitation, S&P p. 88 Remark:

> "There is a close relationship between a small RGA number and diagonal dominance, but unfortunately there are exceptions for plants of size 4×4 or larger, so a small RGA number does not always guarantee diagonal dominance."

Irrelevant at the current preset scale (max 2 loops) but worth a comment in the code.

### 4.4 The trap: never show magnitudes alone

S&P §3.4.3, p. 87, opening paragraph:

> "Note that in Figure 3.8(a) we plot only the magnitudes of `λij`, but this may be misleading when selecting pairings. **For example, a magnitude of 1 (seemingly a desirable pairing) may correspond to an RGA element of −1 (an undesirable pairing).** The phase of the RGA elements should therefore also be considered."

Verified numerically: for `G = [[1,2],[1,1]]`, `Λ = [[−1,2],[2,−1]]`. The diagonal element has `|λ11| = 1.0000` — which reads as a perfect pairing — while the actual value is `−1`, an integrity failure, and the RGA number is **8**. A UI that displays `|λii|` would render this plant green.

**Design rule: display `λii` as a complex number (or as real + imaginary, or magnitude + phase), and drive all colour/severity logic off the RGA number, never off `|λii|`.**

### 4.5 Proposed escalation ladder

The thresholds below are **engineering judgement, not literature values**, except where a source is named. S&P give exactly three anchors — "close to 1" (rule 1), "negative" (rule 2), "5–10 or larger" (C1) — and no UI bands. Label them in code as tunable.

| condition (on `Λ` at the chosen ω) | severity | text |
|---|---|---|
| fewer than 2 loops paired | **N/A** | "Pairing diagnostic needs at least two loops — a single loop has no interaction to measure." (§3.2 Caveat B) |
| `G_S` numerically singular / rank-deficient | **error** | "`G(jω)` is singular at this frequency — RGA undefined. Try a different ω." (§5.3) |
| any `Re(λii) < 0` | **critical** | "Pairing on a negative RGA element. Gain sign reverses when the other loops close; with integral action this is an integrity failure (S&P Thm 10.6)." Qualify if the plant is unstable — §5.4 |
| any `\|λii\| ≥ 5` | **critical** | "Large RGA element — plant is fundamentally difficult to control; avoid inverse-based controllers (S&P C1)." |
| any `\|λii\| < 0.1` | **critical** | "Near-zero relative gain — this loop's input has essentially no direct effect on its output. Check the pairing." |
| `RGAnum ≥ 1.0` | **warn** | "Strong two-way interaction; loops cannot be tuned independently." |
| `0.3 ≤ RGAnum < 1.0` | **caution** | "Moderate interaction; expect loop tuning to interact." |
| `RGAnum < 0.3` | **silent** | — |
| an alternative pairing has a strictly smaller RGA number | **suggest** | "Pairing `(y2←u1, y1←u2)` has RGA number 0.00 vs. 4.00 for the current pairing." |

The last row is the highest-value item in the table and is cheap: by property 2, Λ is computed once, and each candidate pairing is a permutation `Π` to subtract. For p loops there are `p!` permutations — 2 for the Ball-Balancer, 6 for a 3-loop plant, 24 for 4. Enumerate exhaustively up to p ≤ 6 (720) and skip above that.

---

## 5. Integrator Plants and Numerical Guards

### 5.1 Is the RGA well-conditioned at ω > 0 for a pure-integrator plant? Yes — provably, and exactly

The map records that the Ball-Balancer's `A` is pure integrators, that `−C·A⁻¹·B + D` is therefore undefined, and that "its DC gain is genuinely infinite, not merely numerically awkward". All true. **None of it affects the RGA at ω > 0**, and the reason is a theorem, not luck.

From `model_library.cpp:220-231`, the Ball-Balancer is

```
A(0,2)=A(1,3)=1,  B(2,1)=B(3,0)=5g/7,  C(0,0)=C(1,1)=1,  D=0
```

whose transfer matrix is

```
G(s) = (5g/7)/s² · [[0, 1],
                    [1, 0]]
```

— a scalar times a permutation matrix. Now apply two S&P properties:

- **A.4.1 property 5** (p. 527): "The RGA is scaling invariant. Therefore `Λ(D₁AD₂) = Λ(A)` where `D₁` and `D₂` are diagonal matrices." A scalar multiple is a diagonal `D₁ = c·I`.
- **A.4.1 property 2** (p. 527): "`Λ(P) = P` for any permutation matrix."

Therefore `Λ(G(jω)) = Λ(P) = P` for **every** ω > 0. The `1/(jω)²` factor — the entire source of the infinite DC gain — cancels identically.

Verified numerically over ω ∈ [10⁻⁶, 10⁶] rad/s:

| ω (rad/s) | max \|G(jω)\| | Λ | RGAnum (diagonal) | RGAnum (swapped) |
|---|---|---|---|---|
| 1e−6 | 7.007e+12 | `[[0,1],[1,0]]` | 4 | 0 |
| 1e−2 | 7.007e+04 | `[[0,1],[1,0]]` | 4 | 0 |
| 1e+0 | 7.007e+00 | `[[0,1],[1,0]]` | 4 | 0 |
| 1e+3 | 7.007e−06 | `[[0,1],[1,0]]` | 4 | 0 |
| 1e+6 | 7.007e−12 | `[[0,1],[1,0]]` | 4 | 2.2e−16 |

Exact to machine precision across eighteen decades of plant gain, and `cond(G(jω)) = 1` at every frequency (a scaled permutation is perfectly conditioned — it is *badly scaled*, not ill-conditioned; the two are different and only the latter hurts).

**And it detects exactly the failure the map said a magnitude comparison cannot.** The map notes the Ball-Balancer's coupling "is a channel *swap* (`B(2,1)`, `B(3,0)` — input 0 drives output 1's chain), exactly the pairing failure a magnitude comparison cannot detect." The RGA number for the naive diagonal pairing is **4** — the worst possible value for a 2×2 — and **0** for the swapped pairing. `λ11 = λ22 = 0` fires the near-zero-relative-gain alarm from §4.5. This is the strongest possible confirmation that RGA-at-a-frequency is the right replacement for the DC heuristic on this plant.

### 5.2 The general case

A plant whose `G(s)` is *exactly* a scalar transfer function times a constant matrix has a frequency-independent Λ. Most integrator plants are not that clean: they are `G(s) = M(s)/s^k` with `M(s)` non-constant, in which case Λ varies with ω but remains bounded, because the `1/s^k` common factor still cancels under scaling invariance. The general statement worth relying on is: **a common scalar factor in `G(s)` — including poles at the origin, and including a common time delay `e^{-θs}` — cancels out of the RGA exactly.** (This is why the delay in S&P eq. (3.65) does not appear in eq. (3.66).)

The RGA is therefore *structurally* immune to the thing that broke the DC heuristic. What it is not immune to is `G(jω)` being genuinely rank-deficient at the specific ω you chose — a transmission zero on the imaginary axis, or a plant with a genuinely uncontrollable direction. That is a real failure and needs a guard.

### 5.3 Numerical guards

**Guard 1 — reject `ω = 0` for plants with poles at the origin, loudly.** `evalTransferFunction` (`frequency_response.cpp:8-22`) forms `sI − A` and calls `colPivHouseholderQr().solve(b)`. At `s = 0` with a singular `A` this is a singular system, and `ColPivHouseholderQR::solve` **returns a least-squares-ish answer without signalling failure** — no exception, no status flag. You will get finite-looking garbage. Check `A` for eigenvalues at (or very near) the origin before allowing `ω = 0`, or check `sI − A` conditioning inside the RGA path.

**Guard 2 — equilibrate before inverting; it is free.** Scaling invariance (property 5) means you may replace `G` by `D₁GD₂` for any diagonal `D₁, D₂` without changing Λ at all. So normalize away the dynamic range before the inverse:

```cpp
// Λ(D1·G·D2) = Λ(G) exactly  [S&P A.4.1 property 5, p. 527]
// so we may equilibrate G to O(1) before inverting, for free.
Eigen::MatrixXcd equilibrate(const Eigen::MatrixXcd& G) {
    Eigen::MatrixXcd S = G;
    for (int i = 0; i < S.rows(); ++i) {
        double r = S.row(i).cwiseAbs().maxCoeff();
        if (r > 0.0) S.row(i) /= r;
    }
    for (int j = 0; j < S.cols(); ++j) {
        double c = S.col(j).cwiseAbs().maxCoeff();
        if (c > 0.0) S.col(j) /= c;
    }
    return S;
}
```

For the Ball-Balancer at ω = 10⁻⁶ this turns entries of magnitude 7×10¹² into entries of magnitude 1 with no change to the answer. Worth doing unconditionally.

**Guard 3 — test rank with an SVD, not with `det`.** `det(G)` scales as the product of the singular values and for the Ball-Balancer at ω = 10⁻⁶ is `4.9×10²⁵` — a number that tells you nothing about invertibility. Use the reciprocal condition number:

```cpp
Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Gs);          // BDCSVD for larger p
const auto& sv = svd.singularValues();
double rcond = sv(sv.size() - 1) / sv(0);
if (!(rcond > 1e-10)) {                              // also catches NaN
    return { .valid = false,
             .error = "G(jw) is singular or near-singular at this frequency "
                      "(rcond = " + std::to_string(rcond) + "); RGA undefined." };
}
```

`1e-10` is a suggestion, not a sourced value — it leaves ~6 decimal digits of headroom in double precision. Because the input is equilibrated first, `rcond` is a meaningful measure rather than an artefact of units.

**Guard 4 — validate with the row/column sum identity.** Property A2 says every row and column of Λ sums to exactly 1. That is a free, exact, post-hoc check on the whole pipeline — it catches a failed inverse, a bad `evalTransferFunction`, *and* the `.adjoint()`-instead-of-`.transpose()` bug of §1.3:

```cpp
double sum_err = 0.0;
for (int i = 0; i < p; ++i) {
    sum_err = std::max(sum_err, std::abs(L.row(i).sum() - std::complex<double>(1,0)));
    sum_err = std::max(sum_err, std::abs(L.col(i).sum() - std::complex<double>(1,0)));
}
if (sum_err > 1e-8) { /* mark result untrustworthy */ }
```

This should be an assertion in tests and a soft warning in the UI. Note it holds for the **square** RGA only; for the non-square case only the conditional versions of A.4.2 apply (§3.3).

**Guard 5 — do not use `MatrixXcd::inverse()` blindly.** For p ≤ 4 Eigen uses a direct cofactor formula and will happily return `inf`/`NaN` on a singular matrix. Gate it behind Guard 3, or use `FullPivLU` and check `isInvertible()`.

### 5.4 The caveat the ticket does not anticipate: neither integrator plant is DIC-eligible

Pairing rule 2 and Theorem 10.6 both require a **stable** plant. S&P Theorem 10.6, p. 443, opens "Consider a **stable** square plant `G` …", eq. (10.79) reads "A **stable** (reordered) plant `G(s)` is DIC only if `λii(0) ≥ 0`", and Remark 3 on p. 445 is blunt:

> "**Unstable plants are not DIC.** The reason for this is that with all `εi = 0` we are left with the uncontrolled plant `G`, and the system will be (internally) unstable if `G(s)` is unstable."

The Ball-Balancer has four poles at the origin — not asymptotically stable. The Inverted Pendulum is explicitly unstable. **The `λii < 0` ⇒ integrity-failure interpretation is formally inapplicable to both.** S&P's prescribed order of operations, §10.6.8, p. 449:

> "If the plant is unstable, then it is recommended that a lower-layer stabilizing controller is first implemented, at least for the 'fast' unstable modes."

i.e. stabilize first, then run the pairing analysis on the stabilized plant.

**Recommendation:** compute the plant's poles (already available via `computePoleZero`), and when any pole has `Re ≥ 0`, qualify the negative-λ warning: "plant is not asymptotically stable — the DIC integrity result (S&P Thm 10.6) does not apply; interpret this as an interaction measure only." Rule 1 (RGA close to identity at bandwidth) is a diagonal-dominance argument and does not carry the stability hypothesis, so it remains usable — which is fortunate, since rule 1 is the one that catches the Ball-Balancer's channel swap.

---

## 6. Reference Implementation (Eigen 3.4)

### 6.1 Building `G(jω)` on the paired sub-matrix

The key design decision: build `G_S` **in loop order**, so the paired channels land on the diagonal by construction and the comparison target is `I`. This makes the RGA number a one-liner and is exactly S&P's "the rearranged system, with the selected pairings along the diagonal".

```cpp
// src/analysis/rga.h
#pragma once

#include "../linear_system.h"
#include "frequency_response.h"     // evalTransferFunction
#include <Eigen/Core>
#include <complex>
#include <string>
#include <vector>

namespace caliburn {

// One entry of the loop-pairing list from map #1.
struct LoopPair {
    int output_i;   // row of G
    int input_j;    // column of G
};

struct RGAResult {
    Eigen::MatrixXcd lambda;      // p x p, complex  [S&P eq. (3.54)]
    double rga_number = 0.0;      // ||Lambda - I||_sum  [S&P eq. (3.67)]
    double rcond = 0.0;           // sigma_min / sigma_max of the equilibrated G_S
    double sum_error = 0.0;       // max |row/col sum - 1|; should be ~0  [property A2]
    bool valid = false;
    std::string error;
};

// Lambda(G_S(jw)) for the square sub-matrix of the paired channels.
// omega is in rad/s.
RGAResult computeRGA(const LinearSystem& plant,
                     const std::vector<LoopPair>& loops,
                     double omega);

// ||Lambda - Pi||_sum for an arbitrary permutation of the loop order.
// perm[k] == index of the loop whose input is assigned to output-slot k.
double rgaNumberForPermutation(const Eigen::MatrixXcd& lambda,
                               const std::vector<int>& perm);

}  // namespace caliburn
```

```cpp
// src/analysis/rga.cpp
#include "rga.h"
#include <Eigen/LU>
#include <Eigen/SVD>
#include <algorithm>
#include <cmath>

namespace caliburn {
namespace {

// Lambda(D1 * G * D2) == Lambda(G) exactly for diagonal D1, D2.
// [S&P Appendix A.4.1 property 5, p. 527]  -> equilibration is free.
Eigen::MatrixXcd equilibrate(const Eigen::MatrixXcd& G) {
    Eigen::MatrixXcd S = G;
    for (int i = 0; i < S.rows(); ++i) {
        double r = S.row(i).cwiseAbs().maxCoeff();
        if (r > 0.0) S.row(i) /= r;
    }
    for (int j = 0; j < S.cols(); ++j) {
        double c = S.col(j).cwiseAbs().maxCoeff();
        if (c > 0.0) S.col(j) /= c;
    }
    return S;
}

}  // namespace

RGAResult computeRGA(const LinearSystem& plant,
                     const std::vector<LoopPair>& loops,
                     double omega) {
    RGAResult out;
    const int p = static_cast<int>(loops.size());

    if (p < 2) {
        out.error = "RGA needs at least two loops - a single loop has no "
                    "interaction to measure (Lambda is identically 1).";
        return out;                                  // see section 3.2 caveat B
    }
    if (!(omega > 0.0)) {
        out.error = "RGA must be evaluated at omega > 0; omega = 0 is undefined "
                    "for plants with poles at the origin.";
        return out;                                  // see section 5.3 guard 1
    }

    // ---- Build G_S(jw) in LOOP ORDER, so paired channels sit on the diagonal.
    const std::complex<double> s(0.0, omega);
    Eigen::MatrixXcd Gs(p, p);
    for (int r = 0; r < p; ++r) {
        for (int c = 0; c < p; ++c) {
            // row r  <- output of loop r ;  col c <- input of loop c
            Gs(r, c) = evalTransferFunction(plant, loops[r].output_i,
                                            loops[c].input_j, s);
        }
    }
    if (!Gs.allFinite()) {
        out.error = "G(jw) is not finite at this frequency.";
        return out;
    }

    // ---- Guard: equilibrate, then test rank with an SVD (not with det).
    const Eigen::MatrixXcd Ge = equilibrate(Gs);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(Ge);
    const auto& sv = svd.singularValues();
    out.rcond = (sv(0) > 0.0) ? sv(sv.size() - 1) / sv(0) : 0.0;
    if (!(out.rcond > 1e-10)) {
        out.error = "G(jw) is singular or near-singular at this frequency "
                    "(rcond = " + std::to_string(out.rcond) + "); RGA undefined.";
        return out;                                  // see section 5.3 guard 3
    }

    // ---- Lambda = G .* (G^-1)^T
    // NOTE: .transpose(), NEVER .adjoint().  [S&P footnote 5, p. 83]
    // .adjoint() is the conjugate transpose and gives a different, wrong answer
    // whenever G(jw) is genuinely complex - which is always, for w > 0.
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(Ge);
    const Eigen::MatrixXcd Ginv = lu.inverse();
    out.lambda = Ge.cwiseProduct(Ginv.transpose());

    // ---- Self-check: every row and column of Lambda sums to exactly 1.
    // [S&P Appendix A.4.1 property 3, p. 527]
    const std::complex<double> one(1.0, 0.0);
    for (int i = 0; i < p; ++i) {
        out.sum_error = std::max(out.sum_error,
                                 std::abs(out.lambda.row(i).sum() - one));
        out.sum_error = std::max(out.sum_error,
                                 std::abs(out.lambda.col(i).sum() - one));
    }

    // ---- RGA number for the current (diagonal-by-construction) pairing.
    // [S&P eq. (3.67), p. 87]  ||A||_sum = sum_ij |a_ij|
    out.rga_number =
        (out.lambda - Eigen::MatrixXcd::Identity(p, p)).cwiseAbs().sum();

    out.valid = true;
    return out;
}

double rgaNumberForPermutation(const Eigen::MatrixXcd& lambda,
                               const std::vector<int>& perm) {
    // "The RGA number for other pairings is obtained by subtracting 1 for the
    //  selected pairings"  [S&P p. 88] - i.e. ||Lambda - Pi||_sum.
    // Lambda is computed ONCE; permuting G permutes Lambda identically.
    // [S&P Appendix A.4.1 property 2, p. 527]
    const int p = static_cast<int>(lambda.rows());
    Eigen::MatrixXcd Pi = Eigen::MatrixXcd::Zero(p, p);
    for (int k = 0; k < p; ++k) Pi(k, perm[k]) = 1.0;
    return (lambda - Pi).cwiseAbs().sum();
}

}  // namespace caliburn
```

### 6.2 The p×1 (SIMO/MISO) channel screen

Separate function, separate meaning, separate UI panel — see §3.3. Only report it with output scaling applied, per §3.4.

```cpp
// Non-square RGA for a single-input plant: Lambda = G .* (G^dagger)^T
// [S&P eq. (A.84), p. 528].  For an l x 1 column this reduces exactly to
// lambda_i = |g_i|^2 / sum_k |g_k|^2  -- real, non-negative, sums to 1.
// This is NOT a pairing measure: it answers S&P control property C3,
// "which outputs can this input actually reach?"  (p. 89)
//
// WARNING: unlike the square RGA, this is NOT scale-invariant.  Scale each
// output by its allowed deviation first  [S&P section 1.4, p. 5; Remark 1, p. 6];
// see Uhlmann arXiv:1805.10312 and Al Yousuf & Uhlmann arXiv:2106.09766.
Eigen::VectorXd simoChannelShare(const LinearSystem& plant, int input_j,
                                 double omega,
                                 const Eigen::VectorXd& output_scales) {
    const int l = plant.outputs();
    const std::complex<double> s(0.0, omega);
    Eigen::VectorXd share(l);
    double total = 0.0;
    for (int i = 0; i < l; ++i) {
        const double g =
            std::abs(evalTransferFunction(plant, i, input_j, s)) / output_scales(i);
        share(i) = g * g;
        total += share(i);
    }
    if (total > 0.0) share /= total;
    return share;
}
```

### 6.3 Tests worth writing

| test | expectation | source |
|---|---|---|
| triangular `G` | `Λ = I` exactly | S&P property A3 / A.4.1-4, p. 527 |
| diagonal `G` | `Λ = I` exactly | same |
| permutation `G = [[0,1],[1,0]]` | `Λ = G` | S&P A.4.1-2, p. 527 |
| `G = [[1,1],[0.4,−0.1]]` (blending, Example 3.9) | `Λ = [[0.2,0.8],[0.8,0.2]]` | S&P p. 85 |
| FCC steady state, S&P eq. (3.64) | `Λ ≈ [[1.50,0.99,−1.48],[−0.41,0.97,0.45],[−0.08,−0.95,2.03]]` | S&P p. 85 |
| distillation, S&P eq. (3.71) | `Λ ≈ [[35.1,−34.1],[−34.1,35.1]]` | S&P p. 89 |
| **complex `G(jω)`, Example 3.11 at ω = 0.01** | `Λ ≈ [[0.24+0.02i, 0.76−0.02i],[0.76−0.02i, 0.24+0.02i]]` | S&P eq. (3.66), p. 86 — **the only test that catches `.adjoint()`** |
| Ball-Balancer at any ω > 0 | `Λ = [[0,1],[1,0]]`; RGAnum(diag) = 4, RGAnum(swap) = 0 | §5.1 above |
| Ball-Balancer over ω ∈ [1e−6, 1e6] | Λ invariant to machine precision | §5.1 above |
| `G = [[1,2],[1,1]]` | `λ11 = −1`, `\|λ11\| = 1`, RGAnum = 8 | §4.4 — the magnitudes trap |
| any random non-singular `G` | all rows and columns of Λ sum to 1 | S&P A.4.1-3, p. 527 |
| any `G`, scaled `D₁GD₂` | Λ unchanged | S&P A.4.1-5, p. 527 |
| 1-loop pairing | returns `valid = false` with the "needs two loops" message | §3.2 caveat B |

The Example 3.11 row is the load-bearing one. Every other test in the table uses a real `G`, where `.transpose()` and `.adjoint()` agree.

---

## 7. Summary Against the Ticket's Five Questions

| # | question | answer |
|---|---|---|
| 1 | Hadamard / inverse-transpose form; complex handling | **Confirmed** as `Λ = G ∘ (G⁻¹)ᵀ` (S&P 3.54, A.77). Λ is formed **from the complex matrix directly**; magnitudes are taken only for display and for the RGA number (S&P Table 3.1, eq. 3.66). Magnitudes-first is not an approximation — on S&P's own worked example it flips λ₁₁ from +0.24 to **−0.47**, inverting the diagnostic. Use `.transpose()`, never `.adjoint()` (S&P footnote 5, p. 83) |
| 2 | which ω | **Closed-loop bandwidth / crossover** (S&P pairing rule 1, pp. 85 and 450). S&P explicitly warn that DC is misleading and give an example where `Λ(0) = Λ(j∞) = I` yet the correct pairing is the reverse one. But **pairing rule 2 is anchored at DC and only at DC** (S&P p. 450, Thm 10.6) — one ω cannot serve both rules. Recommend: user-chosen ω for rule 1, plus a sweep plot of the RGA number over the existing Bode grid (S&P Fig. 3.8b). No primary source prescribes an automatic default ω |
| 3 | subset vs. pseudo-inverse | **Square sub-matrix is sound** for pairing — it is exactly Bristol's two extreme cases restricted to the loops actually closed (derived from S&P 3.56–3.58; S&P do not state the subset case, so this derivation is ours). Caveats: `Λ(G_S)` is **not** a sub-matrix of `Λ(G)`; and for a **single loop it is identically 1 and says nothing**. The non-square pinv RGA answers a *different* question — channel selection, S&P C3 p. 89 — and loses scale invariance, with Uhlmann (arXiv:1805.10312) and Al Yousuf & Uhlmann (arXiv:2106.09766) demonstrating that a unit change alters the pairings it recommends |
| 4 | thresholds | λ≈1 ideal; 0<λ<1 gain rises when loops close (λ=0.5 is the 2×2 ambiguity point); λ>1 gain falls; λ=0 no direct effect; **λ<0 integrity failure** (S&P Thm 10.6); **λ ≥ 5–10 fundamentally difficult plant** (S&P C1). Drive the UI off the **RGA number** `‖Λ − Π‖sum` (S&P 3.67), never off `\|λii\|` — S&P p. 87 warn `\|λ\| = 1` can be `λ = −1`, verified here as a plant with RGA number 8 that a magnitude readout renders green |
| 5 | integrators & guards | **Well-conditioned, provably.** `G(s) = (5g/7)/s² · [[0,1],[1,0]]` for the Ball-Balancer, so by scaling invariance and `Λ(P) = P`, `Λ = [[0,1],[1,0]]` **exactly at every ω > 0** — verified to machine precision over 18 decades. It catches the channel swap: RGAnum 4 (diagonal) vs 0 (swapped). Guards: reject ω=0 for origin poles (`colPivHouseholderQr` fails silently); equilibrate before inverting (free, by scaling invariance); SVD-based rcond, not `det`; row/column-sum self-check |

### Things that contradict or complicate the ticket's assumptions

1. **The sub-matrix RGA is vacuous for 3 of the 5 presets.** Inverted Pendulum, Quarter-Car and Double Mass-Spring are 1-input, so at most one loop, so `Λ = [1]` always. The map treats "RGA on the paired sub-matrix" as the general pairing diagnostic; it is only informative for ≥ 2 loops, which today means the Ball-Balancer alone. The 1-loop case needs the non-square channel screen (§3.3, §6.2) instead, and the UI must say "not applicable" rather than showing a green `1.00`.
2. **One ω cannot serve both pairing rules.** Rule 1 wants the bandwidth; rule 2 is a steady-state DIC theorem. The ticket's "which ω" framing assumes a single choice.
3. **Neither integrator preset is DIC-eligible.** Theorem 10.6 requires a *stable* plant; S&P state flatly that "unstable plants are not DIC" (p. 445). The Ball-Balancer (poles at origin) and Inverted Pendulum (unstable) both fall outside it, so the `λii < 0` ⇒ integrity-failure message must be qualified for them. Rule 1 is unaffected and still works.
4. **The pseudo-inverse generalization is contested in the current literature.** S&P endorse it (A.4.2) but two recent papers show it is not unit-invariant and that this changes recommended pairings in practice. This does not affect the square case at all, but it means any non-square readout must be preceded by explicit output scaling.
5. **The map's open question "whether the preset library needs a square MIMO example with finite DC gain" is now answerable: no.** The Ball-Balancer is an excellent RGA demo precisely because its infinite DC gain cancels exactly and its coupling is a pure channel swap that the RGA number resolves as 4-vs-0. The map speculated this "may become moot once RGA-at-frequency lands" — it does.

---

## 8. Sources

**Primary — textbook**

- S. Skogestad and I. Postlethwaite, *Multivariable Feedback Control: Analysis and Design*, 2nd ed., Wiley, 2005. §3.4 (pp. 82–91), §10.6 (pp. 429–454), Appendix A.4 (pp. 526–529), §1.4 (pp. 5–7). Read from the author-hosted PDFs: [full book](https://skoge.folk.ntnu.no/book/ps/bookall.pdf), [ch. 3 + ch. 10 extract](https://skoge.folk.ntnu.no/vgprosessregulering/papers-pensum/ch3-ch10-2ndedition.pdf).

**Primary — papers**

- E. H. Bristol, "On a new measure of interaction for multivariable process control", *IEEE Trans. Automatic Control*, vol. 11, no. 1, pp. 133–134, 1966. [DOI 10.1109/TAC.1966.1098266](https://doi.org/10.1109/TAC.1966.1098266). **Paywalled; not read directly.** Used only via S&P §3.4.1, which reproduces the derivation.
- J. Uhlmann, "On the Relative Gain Array (RGA) with Singular and Rectangular Matrices", [arXiv:1805.10312v3](https://arxiv.org/abs/1805.10312), 2019. Read in full.
- R. Q. Al Yousuf and J. Uhlmann, "On Use of the Moore-Penrose Pseudoinverse for Evaluating the RGA of Non-Square Systems", [arXiv:2106.09766v2](https://arxiv.org/abs/2106.09766), 2022. Read in full.

**Cited via S&P, not read directly**

- P. Grosdidier, M. Morari and B. R. Holt (1985) — the DIC / negative-RGA theorem, S&P Theorem 10.6, p. 443.
- C. Chang and C. Yu (1990) — the non-square RGA extension, S&P Remark, p. 529.
- S. Skogestad and M. Morari (1988) — DIC, S&P Remark 1, p. 444; the 2×2 equivalence, eq. (10.82).
- W. Yu and M. Fan (1990) — the 3×3 DIC condition, S&P eq. (10.83).
- C. Johnson and H. Shapiro (1986) and E. Wolff (1994) — iterative RGA convergence, S&P p. 88.

**Checked and found absent**

- `python-control` — no `relative_gain_array` / `rga` in the `main` source tree or the [ReadTheDocs global index](https://python-control.readthedocs.io/en/latest/genindex.html).
- MATLAB Control System Toolbox / Robust Control Toolbox — no documented RGA function. Only MATLAB Central File Exchange community submissions exist, which are not primary sources. The reference implementation is S&P's own Table 3.1 snippet.

**Numerical verification**

All numbers in §1.2, §4.4 and §5.1 were computed from scratch in pure-Python complex arithmetic (no numpy available in the environment) using the state-space data in `src/analysis/model_library.cpp` and the models printed in S&P eqs. (3.65) and (3.71). The one unresolved discrepancy — S&P's printed eq. (3.66) versus the value recomputed from eq. (3.65), differing in the third decimal — is flagged in §1.2 rather than explained.
