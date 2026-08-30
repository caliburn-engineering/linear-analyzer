# Ball-Balancer — Context

One application over two halves: a 3-RRS table balancing a ball under real
rolling dynamics, and the analyzer that models a plant, closes a loop around it,
and shows the consequences across Bode, Nyquist, pole-zero and step-response
views.  The repository is still named `linear-analyzer` — the merge (issue #13)
renamed no directory.

## Glossary

Terms below are the project's ubiquitous language. Where a term has a plausible
synonym that is *wrong* here, the synonym is recorded with the reason — drifting
to it silently reintroduces a resolved confusion.

### Loop

A pairing of one plant output to one plant input, carrying its own gains. The
controller is a **list** of loops (`std::vector<Loop>`), not a per-channel
diagonal block — this is what lets square, SIMO and MISO plants share one data
model. A square plant with `m` loops reproduces the diagonal case.

Fan-out is legal in both directions; many-to-one is the primary MISO structure.
An empty loop list is equivalent to `ControllerType::None`.

Decided in [#2](https://github.com/caliburn-engineering/caliburn/issues/2).

### RGA (Relative Gain Array)

`Λ = G ∘ (G⁻¹)ᵀ` on the **complex** `G(jω)`, evaluated at the intended
closed-loop bandwidth — not at DC. A *pairing* diagnostic: it grades whether the
loops are matched to the right channels.

Defined only on the **square sub-matrix** of the paired channels, and therefore
only when the loop list's distinct paired outputs and distinct paired inputs are
equal and at least 2. Outside that, see **Channel Share**.

Severity comes from the RGA number `‖Λ − Π‖ₛᵤₘ`, never from `|λᵢᵢ|`.

Decided in [#4](https://github.com/caliburn-engineering/caliburn/issues/4).

### Channel Share

The fraction of one input's squared effect landing on each output, at a given ω:

```
λᵢ = |gᵢ(jω)|² / Σₖ |gₖ(jω)|²
```

Real, non-negative, sums to exactly 1. Shown wherever the RGA is vacuous or
undefined — a single-input plant, a multi-input plant with only one loop, or a
MISO fan-out whose paired sub-matrix is non-square.

**Channel Share is not a pairing measure and must never be presented as one.**
It answers "which outputs does this input reach", not "are these the right
pairings". It is also **not scale-invariant**, unlike the square RGA: it is
reported against user-supplied per-output scales (maximum allowed deviation),
and the readout states when those scales are still at their defaults.

Because it is normalised it cannot distinguish "reaches both outputs strongly"
from "reaches neither", so it is always displayed alongside `|gᵢ(jω)|` in dB.

Decided in [#8](https://github.com/caliburn-engineering/caliburn/issues/8).

> **Not "reach", "reachability", or "control authority".**
> *Reachability* is the textbook synonym for controllability, which this app
> already computes (`checkControllability`) and displays in the System Properties
> panel — the collision would be direct. *Control authority* asserts the input is
> commanded, which is false for the Quarter-Car preset, whose single input is the
> road profile disturbance.

### Structurally dead channel

A plant channel whose `|gᵢ(jω)|` is below tolerance at **every** ω on the Bode
grid — a genuine structural disconnect, as opposed to a numerical dip at one
frequency. Positive output rescaling cannot turn a zero into a nonzero at any
frequency, so this test is unaffected by the scale-dependence of Channel Share;
a threshold on the share itself, or a ratio against the column maximum, would
not be.

A loop pairing an output to an input that cannot affect it is the Ball-Balancer
failure case: `B(2,1)` / `B(3,0)` mean input 0 drives output 1's chain, so the
seeded identity pairing is dead.

Dead channels are suppressed with a stated reason rather than drawn with
misleading fallbacks — `arg(0)` renders a flat 0° phase, and the pole-zero pencil
returns `eig(A)` as coincident zeros, which reads as full cancellation.

Decided in [#8](https://github.com/caliburn-engineering/caliburn/issues/8) and
[#9](https://github.com/caliburn-engineering/caliburn/issues/9).

### Pairing grid

The controller section of the model panel: a `p × m` grid of cells, one per
plant channel, that is **simultaneously the pairing editor and the diagnostic
readout**. Clicking a cell creates or destroys the loop on that channel; each
cell shows the RGA `λᵢⱼ` (or the Channel Share, where the RGA is undefined) over
`|gᵢⱼ|` in dB. A leading column carries the per-output scale.

It is the only editor, so the loop list has no add, remove or reorder controls,
and grid order *is* loop order. One consequence: an exact duplicate pairing is
unreachable through the UI, so **#2's duplicate warning needs no surface** —
though `buildLoopController` still sums duplicates, which is correct for a loop
list built programmatically.

**Paired-ness is the cell border; severity is the cell fill.** They are separate
visual channels because a single one cannot carry both, and severity applies to
RGA cells only — Channel Share is never coloured (see above).

Decided in [#6](https://github.com/caliburn-engineering/caliburn/issues/6).

> **Not the "coupling matrix" or the "RGA table".**
> The grid outlives the RGA: on a single-input plant the same widget shows
> Channel Share, and on a 1×1 plant it is one cell with no diagnostic at all.
> Naming it after one of its readouts hides that it is first an editor.

### Plate view

The ball-balancer half of the application: the 3-RRS table, the ball on it, the
3D scene, and the Plate Control / Plate Plots panels.  It owns no window, no
ImGui context and no main loop — the analyzer's entry point owns all three, and
reconciling that was the whole of the merge.

Drawn into the **central dock node**, which is deliberately left empty so the
passthru dockspace shows the scene through it.  A window docked there would
paint over the plate.

Desktop only.  The Emscripten build compiles the analyzer half alone until the
web port (issue #15).

### Tilt convention (phi/theta vs alpha/beta)

The two halves name the plate's tilt differently, and the mapping is neither
identity nor a simple swap:

```
alpha = -phi        beta = theta
```

`TableKinematics` uses `R = Ry(theta) * Rx(phi)`, so the plate's local +Y edge
*rises* with phi while its local +X edge *falls* with theta.
`RollingBallDynamics` accelerates the ball by `+g*sin(alpha)` along y and
`+g*sin(beta)` along x.  Theta therefore carries straight through and phi has to
flip.

A wrong mapping compiles, and still moves the ball — it just rolls uphill.  That
is why the convention is pinned by `tests/test_ball_sim.cpp` rather than left to
the eye, and why it lives in one function, `ballTiltFromPose`.

### Plate frame

The ball's state `[x, y, vx, vy]` is expressed in the **plate's own frame**, not
the world.  Rendering lifts it into the world with the same rotation that places
the leg attachment points: `p = c + R * (x, y, r_ball)`.

The plate is a **disc**, so the ball leaves it at a radius of
`R_table - r_ball`.  `RollingBallDynamics::on_plate` tests the inscribed
*square* and is wrong here; `ballOnPlate` is the test that matches the geometry
being drawn.
