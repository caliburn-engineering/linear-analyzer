# Linear Analyzer — Context

Interactive analysis of linear time-invariant systems: state-space plants, Tier 1
compensators, and the frequency/time/pole-zero views over them.

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
