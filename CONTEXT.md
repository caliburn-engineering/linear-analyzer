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

Both targets since [#15](https://github.com/caliburn-engineering/caliburn/issues/15).
The renderer draws through the GL subset shared by OpenGL 3.3 core and WebGL2 —
vertex arrays, buffer objects, array draws — so the port was the shaders and
nothing else: `#version 330 core` and `#version 300 es` plus a precision
qualifier, supplied as a separate source string ahead of one shared body.

`glLineWidth` above 1.0 is honoured or ignored at the driver's discretion under
WebGL2, so line weight is decorative here and carries no information: every
element of the scene is told apart by colour.

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

### Balance loop

The closed loop that the LQR design surface exists to produce: `u = -K(x - x_ref)`
evaluated once per frame against the **nonlinear** plate, in `auto_balance.h`.
`x` is the cascade plant's state read straight off the simulator — three leg
deviations from the linearisation point, then the ball's `[x, y, x', y']` in the
plate frame — and `u` is a leg-angle command, so nothing adapts tilt to legs
anywhere.

Engaged from Plate Control, and only while the model panel is offering a gain
that was solved against **this** plate: the right controller type, a successful
solve, a `3 x 7` gain, a ball actually being simulated, and a plant whose
geometry *and gravity* match the simulated one.  Gravity is a cascade parameter
and the plate's is fixed, so a gain designed on the moon is exactly as wrong as
one designed for longer legs and is refused the same way.  Losing any of those
drops the loop rather than freezing the last command, because a stale gain is
not a controller — and `setDesign` is the single place that drop happens, so a
draw pass cannot disagree with the step about who is driving.

The plant a cascade parameter list describes is built by `cascadeMechanism` /
`cascadeGravity` / `cascadeServoTau` / `cascadeHomeLegAngle`, and the model
builder, the application and the tests all read those same functions.  A second
transcription is how the check comes to reject two identical plates — it did,
once, over a `float` slider widening to `0.3000000119`.

Decided in [#16](https://github.com/caliburn-engineering/caliburn/issues/16).

> **The sliders and the loop cannot both be live.**
> Plate Control's three servo sliders and `u = -Kx` write the same array, so
> while the loop is engaged the sliders go read-only and keep displaying what
> the controller is asking for — the most direct view of the loop working.
> Steering moves to a **ball-position setpoint**, which is a quantity the
> design actually regulates.  A ball at rest anywhere on a flat plate is an
> equilibrium of this plant, so that setpoint enters as a reference state and
> needs no feedforward at all: the regulator is already a tracker.

### Leg command vs leg angle

Two arrays, since the loop was closed.  `alpha_cmd_deg_` is what the sliders,
the animation or the controller **ask for**; `alpha_deg_` is where the legs
actually **are**, one first-order lag behind at the plant model's own `tau`.

The lag applies to **every** writer, not just the loop.  It is a property of
the plate rather than of the controller, and making it conditional would mean
the manual plate and the modelled plate are two different machines.  Home/Low/
High therefore command rather than teleport, and the Animation amplitude is the
amplitude *asked for* — attenuated about 1% at its default 0.5 Hz, which is
what a real servo does.  `snapServos` is the one path that writes both arrays,
because a reset is not driving anything.

Before this the slider *was* the leg angle, and closing `u = -Kx` around that
is an algebraic loop on the leg states — the gain would have been reacting to
its own command one frame earlier.  The servo dynamics the cascade model claims
had to become real in the simulator for the design to mean anything.

The lag is integrated exactly, `alpha <- cmd + (alpha - cmd)*exp(-dt/tau)`, not
by forward Euler: the plate runs at a fixed 60 Hz and tau defaults to 0.05 s,
three steps per time constant, and the tau slider goes lower still.

### Assembly mode

Which of the two solutions of the 3-RRS constraint equations a pose is.  The
mechanism is **built upward** — table above its knees — and there is a second,
**folded** root with the table lying flat on the base.  For this plate the
folded one is not an approximation or a numerical artefact: `R_ground ==
R_table` with `L1 == L2` makes `(phi, theta, z_c) = (0, 0, 0)` satisfy every
leg's length constraint exactly, at every servo angle.  It is a single pose,
the same one whatever the legs are doing.

A residual norm cannot tell the two apart, so **a converged forward-kinematics
solve is not by itself an answer about a mechanism**.  `forward_kinematics` is
a Newton solver and will return whichever root its seed is nearest;
`solve_pose` is the one that stays on the built assembly, re-solving from the
analytic level pose when the warm start lands on the folded one.

This is why issue #22 was invisible for so long.  A hard manoeuvre carried the
pose down through the knee plane, Newton settled on the folded root, and every
frame after that was seeded from a pose that was *already* a root — so it
converged in one iteration with a residual of 1e-11 and never left.  The plate
then had no tilt authority at all and every ball put on it rolled off.  Sixty
seconds in, and permanent.

The floor separating them is a tenth of the mean knee height.  Measured over
300k samples against assemblies verified by IK round-trip, the built population
runs from 0.163 to 2.00 mean knee heights and the folded root sits at zero, so
the floor clears the lowest real pose by 1.6x and the folded one by everything.
Both ends are pinned by `tests/test_assembly_mode.cpp`.

Decided in [#22](https://github.com/caliburn-engineering/caliburn/issues/22).

> **Not "the FK failed" and not "a singularity".**
> Nothing failed: the solve converged, to a real root, of the right equations.
> And a singularity is a place where the Jacobian loses rank — this is a
> perfectly well-conditioned second solution. Calling it either invites the fix
> that was already there and did not work: a tighter tolerance, or a retry
> keyed on convergence.

### Workspace vs. servo box

The three servo travel limits form a **box**; the set of leg triples the plate
can actually be assembled at is not one.  Barely half the box has any assembly,
and a *single* leg taken to either of its shipped stops with the other two at
home already has none — the mechanism's real range about the home pose is far
narrower than its servos'.

So clamping each leg to its own travel is not enough.  A per-leg clamp can hand
back a triple the plate cannot make, and a plate that cannot make its command
has no pose: it freezes at whatever tilt it last held, and a frozen tilted plate
rolls the ball off.  That was the other half of #22 — nearly two thousand frames
of a ten-minute run.

**The scatter was the folded assembly, not the controller.**  Before the fix
the loss counts jumped about with no pattern — 0.22 m/s of kick lost the ball
where 0.25 did not, and servo travel of ±35 degrees lost 61 balls in five
minutes where ±45 lost none and ±60 lost 65.  That is not something a
stabilising loop does, and it is what said the cause was discrete rather than
dynamic: whether a particular trajectory happened to carry the pose down
through the knee plane on some frame is a yes-or-no event, and everything
downstream of a yes was already broken.  Travel limits and kick magnitude
changed *which* trajectories crossed, not how much authority the loop had.

With the assembly pinned, the response is monotonic in the disturbance the way
it should always have been: at 360 swept directions the shipped tuning loses
none at 0.35 m/s, 2 at 0.42, 14 at 0.50 and 259 at 0.60.  A curve, not a
scatter — and one that says something true about the plant.

`retreatToWorkspace` is the one answer, used at both seams: pull the target back
toward a triple known to be assemblable until it is too.  Giving up magnitude
and keeping direction is what saturation ought to do, and the returned point is
always one that was actually tested, so the workspace does not have to be convex
for it to be right.

Both seams are needed, which is the part worth remembering.  `legCommand` clips
the command against the level pose; `stepServosOnPlate` clips the *step*
against where the legs already are.  Clipping only the command leaves the
straight line the first-order lag travels along, and that line can leave the
workspace even when both of its ends are inside it.  Two kick directions in
every 360 were exactly that: one frame, mid-flight, with no assembly.

The clip is not a rare corner.  Under the shipped tuning a single attract kick
saturates the legs from every direction tested and is workspace-clipped from
177 of 180 of them; under a deliberately aggressive tuning, from all 180.  It
runs on essentially every disturbed frame.

Decided in [#22](https://github.com/caliburn-engineering/caliburn/issues/22).

> **The plate has less authority than the linear design believes.**
> The cascade model knows about servo travel and nothing about the workspace,
> so the gain will ask for tilts the mechanism cannot make, and against a
> disturbance large enough it will lose the ball rather than recover it.
>
> There is no clean ceiling: the slivers of direction where it cannot come
> back get narrower as the disturbance shrinks rather than stopping at a
> threshold.  At 5760 swept directions, 0.29 m/s of ball velocity loses six,
> 0.28 loses two, and 0.26 loses none.  Attract mode's kick sits at 0.26 for
> that reason, and nobody has proved there is not a narrower sliver still.
> Designing a gain that respects the workspace is not done, and is the thing
> that would replace this argument with a guarantee.

### Attract mode

The application driving itself for a visitor who has not arrived yet: the
balance loop engaged as soon as a gain exists, and the ball already tracing a
120 mm circle at a ten-second lap.  Named after the arcade cabinet's demo reel,
which solves exactly this problem.

A motionless canvas is indistinguishable from a broken build, and a still
balanced ball reads as a photograph of one.  A visitor gives the page about ten
seconds and will not go looking for a play button.

**This used to be a disturbance schedule, and the change is worth recording
because the original reasoning was sound.**  The demo opened with the ball 72 mm
off centre and kicked it every four seconds, on the argument that what
distinguishes a control system from an animation is *recovery*: something
disturbs the ball and the loop puts it back.  In practice it read as a fault.

Two reasons, and the first is the one that matters:

- **A kick is legible only against a still baseline.**  A ball parked at the
  centre, nudged, returning — that is a story.  But the still moment between
  kicks is also the moment the page looks like a static image, so the demo
  spent most of its time looking broken in order to make the other part legible.
  A circle needs no baseline: at every instant the ball is somewhere it was
  told to be and the plate is visibly working to keep it there.  The loop is
  not demonstrating that it *can* respond, it is demonstrating that it *is*
  responding, continuously — which is also the more honest claim, since it is
  the one running every frame.
- **The opening displacement was a step input.**  Handing a state feedback
  72 mm of error at `t = 0`, with the legs exactly at home and no servo history
  to smear it, is a step: measured, the legs swung **20.7 degrees apart within
  three frames** and the table dropped 4.5 mm before settling inside 250 ms.
  Correct, and it looked like the mechanism glitching.  Starting the ball ON
  the path at the path's own velocity makes both the position and the velocity
  error zero at `t = 0`; the same measurement then gives a **3.0 degree**
  opening swing with the table height not moving at all.  Pinned by
  `test_the_opening_does_not_slam_the_legs`.

What is given up is stated plainly, because it was a real acceptance criterion:
**the demo no longer shows disturbance rejection unprompted.**  The visitor can
still see it — the Nudge buttons are right there — and the property is still
pinned by `tests/test_attract_mode.cpp`, which sweeps a 0.26 m/s disturbance
over all 720 directions and checks the ball comes home from every one.  It is
no longer *performed*.  That sweep also moved: the kick used to be production
code the demo delivered, and is now a test fixture, because what it tests is a
property of the plant and the gain rather than of the demo.

Two decisions carry the opening:

- **The opening state is a state, not a performance.**  It is set once in the
  constructor and then the visitor owns it.  The old mode had to watch for a
  visitor arriving so it could stand its kicks down, and never resume; a circle
  has nothing to stand down, so the watching went with it.
- **The circle rather than a cornered shape.**  Its tracking error is smooth,
  so the opening reads as competence rather than as the ball stumbling at every
  corner.  The corners are the better demonstration and are one dropdown away —
  but they are an argument the visitor should choose to hear, not the first
  thing they see.  120 mm at a ten-second lap is 75 mm/s, a third of the speed
  cap, tracked to 3.7 mm.

`setDesign` is where the loop first closes, and it has to be: on the opening
frame the LQR solve has not run, so a `balance_engaged_` set true in the
constructor would be cleared by the stale-gain drop and never set again.
Engaging on the first *usable* design is what "already stabilising at load"
actually amounts to.  It fires **once** — without the latch it would re-engage
a loop the visitor had deliberately dropped, on the very next frame, which is
the demo arguing with the person using it.

The application also opens with `ControllerType::LQR` selected and the
closed-loop trace visible, since a page whose whole claim is that the loop is
closed should not open its pole-zero map on the open-loop poles.

Decided in [#17](https://github.com/caliburn-engineering/caliburn/issues/17),
revised for the opening path in
[#24](https://github.com/caliburn-engineering/caliburn/issues/24).

> **Not "demo mode", "idle mode" or "screensaver".**
> All three name a state the application is *stuck in* until something releases
> it, and two of them imply a substitute for the real thing — a canned loop
> playing where the product would be.  Nothing here is canned: the plate, the
> gain and the rolling ball are the same ones the visitor gets, and the only
> difference is who is choosing the setpoint.  "Attract" also says what the
> mode is *for*, which is the test any addition to it has to pass.
>
> Prose may still call the running page "the demo" — that is the artefact, not
> the mode.

### Contact, and the two phases

The ball is no longer glued to the plate.  `RollingBallDynamics` still models
the rolling half — it is a good model of a ball in contact — but whether the
ball IS in contact is a question it could not ask, because it had no normal
force and no vertical state.

`ball_contact.h` asks it.  `N/m = g(n·z) + n·c̈ + n·[ω̇×(Rs) + ω×(ω×(Rs))] +
2n·(ω×Rṡ)` — gravity's share along the normal, the plate being driven up or
down under the ball, the plate's rotation swinging the contact point, and
Coriolis.  **`N/m ≤ 0` is separation**: a surface can push a ball and never
pull it.

This is not a corner case.  Measured under the shipped tuning and a 0.26 m/s
disturbance, the ball separates in 30 of 72 directions — briefly, a few frames,
hopping 6.7 mm.  The plate heaves hard because the legs do (`z_c = 0.30·sin α`,
so a 20° leg swing is 74 mm of table in a tenth of a second), and a loop
rejecting a disturbance slams the legs.

**Distinguish the tuning from the demo here, since they parted company.**  The
shipped *tuning* hops, as above, and a visitor pressing Nudge will see it.  The
shipped *demo* no longer does, because it no longer kicks the ball: tracing a
gentle circle never separates it, and `test_ten_minutes_unattended_never_loses_the_ball`
asserts exactly that — zero airborne frames in ten minutes.  Both facts are
pinned, and neither implies the other.  See *Attract mode*.

**Two phases, two frames, and the frame is not a detail.**  Rolling is natural
in the plate frame — that is where `RollingBallDynamics` integrates and where
the cascade plant's state vector is written, and neither should change because
the ball can leave.  Flight is natural in the world frame, because "only
gravity acts" is a statement about an inertial frame.  Each phase keeps its own
truth and `plateFrame` converts on demand.

A launch takes the ball's **full** world velocity, including the plate's motion
at the contact point — which is most of it when the legs are slamming.  Landing
is inelastic: a real ball bounces, but a restitution coefficient is a number
nobody here has measured, and the behaviour this exists to show does not depend
on it.

**Rolling resistance is scaled by that same normal force, not by `g`.**  It is a
normal-force effect — the contact patch deforms in proportion to how hard the
surface is pressed — so the retarding acceleration is `c_rr·(N/m)`, and
`c_rr·g` was only ever the level-plate-at-rest special case.  `derivatives` and
`stepBall` take it as an argument; `quasiStaticNormalAccel` supplies `g(n·z)`
for the one caller that cannot have the real thing, the frame where the rates
are not trustworthy, and for linearising about the flat equilibrium.

The correction is one of principle rather than of magnitude, and it is worth
being plain about that.  It moves the separation count not at all — 30 of 72
directions either way, and the peak hop by 0.0001 mm — because separation is
decided by `normalAccel`, which friction does not enter, and because near
separation `N` is small and so is anything scaled by it.  What it does change
is the dead band's derivation, from `asin(c_rr)` to `atan(c_rr)`; see below.

Decided in [#23](https://github.com/caliburn-engineering/caliburn/issues/23).

> **A separation is a claim about the plate's velocity, so it is only as good
> as that velocity.**
> The rates come through the velocity Jacobian, `-J_pose⁻¹ J_alpha`, which
> amplifies without bound near a singularity: an over-aggressive gain drives
> the plate to a Jacobian condition number of 2464, and the rates that come
> back would launch the ball a metre into the air off a 0.26 m/s nudge.  That
> hop is arithmetic, not physics.  `PlateMotion::rates_trustworthy` carries the
> credibility and the contact test declines to act without it, using the
> threshold the application already shows the user — condition 20, where its
> own readout turns red and says "Poor".  This is #22's lesson again: do not
> act on a solve you cannot validate.

> **An over-aggressive tuning does not merely overshoot — it throws the ball
> off the plate.**
> Q on ball position at 2000 loses the ball in 54 of 720 kick directions, with
> launches reaching 689 mm.  Before this the same tuning looked merely fast,
> because a glued ball cannot be thrown.  That is the honest ceiling on #19's
> aggressive preset.

### Trajectory tracking

A moving setpoint: circle, square, triangle, or the fixed point the loop has
always had.  Balancing at the centre proves stability; tracing a shape proves
*tracking*, and it is far more legible to someone who does not read a
pole-zero map.

The path writes `sp_x_mm_` / `sp_y_mm_` rather than going round them, so the
control law is untouched and the sliders become the readout of where the ball
is being sent — the same arrangement the servo sliders have under the loop.
Polygons are traversed at constant **speed**, not constant angle: a corner is
where the interesting behaviour is, and sweeping an angle would crawl through
it.

**Velocity feedforward, decided by measurement.**  The reference state has
always claimed the ball should be at the setpoint *and stationary*, which is
false the moment the setpoint moves — so the loop spent its effort fighting the
motion it was asked for.  Subtracting the path's velocity from the measured one
is the same arithmetic as putting it in `x_ref`, since only the difference
enters `u = -K(x - x_ref)`.  Measured on a 120 mm circle at a ten-second lap:
**3.66 mm of mean error with it, 35.52 mm without** — nearly ten times, for two
lines.  Without it the ball does not follow the circle so much as sit inside it.

The corner survives that: a square is tracked to 14.4 mm at its corners against
a circle's 5.1 mm on the same size and lap.  That is a bandwidth limit made
visible, and it is the point of offering cornered shapes at all.

Decided in [#24](https://github.com/caliburn-engineering/caliburn/issues/24).

> **The sliders are bounded, and there are two bounds because there are two
> ways to ask for the impossible.**
> `kMaxSetpointSpeed` (250 mm/s) binds the large paths — beyond about 400 mm/s
> the plate loses the ball outright — and `kMinLapSeconds` (2 s) binds the
> small ones, where 250 mm/s round a 20 mm circle is a half-second lap and an
> angular rate the plate cannot follow however short the distance.  Neither
> implies the other.  A demo whose controls include a setting that breaks it is
> not offering a choice, it is offering a trap.

### Rolling-friction dead band

`RollingBallDynamics` opposes motion with `c_rr * (N/m) * sign(v)`, so a plate
tilted by less than `atan(c_rr) = 0.573 deg` **cannot start the ball moving at
all**: the tangential pull is `g sin(t)` and the resistance `c_rr g cos(t)`, so
the band closes where `tan(t) = c_rr`.

The derivation moved with #23 and the number did not.  Before the normal force
was computed, the resistance was scaled by `g` outright, which put the edge at
`asin(c_rr)` instead — 0.5729673 deg against 0.5729387 deg, a difference of
2.9e-5 deg, or one twenty-thousandth of the tolerance the test pins it to.  The
two agree to four decimals for any `c_rr` this small, so no test number changed;
what changed is which of them is *true*.  A state feedback has no integral term and therefore no way out: it
parks the ball wherever the tilt it is asking for falls inside the band.

This is not modelled in the linear plant at all, and it is why the LQR designer
does not open on unit weights.  Under `Q = R = I` the residual is 37 mm off
centre with the plate sitting at the dead-band edge — a loop that looks broken
while behaving exactly as designed.  The residual scales inversely with the
position gain, so the shipped default weights the ball position and lands
inside 2 mm.  Both the good case and the dead band itself are pinned by
`tests/test_auto_balance.cpp`; the number is not a memory.

> **Not "stiction" and not "the loop has steady-state error".**
> Stiction is a break-away force distinct from the sliding one, which this
> model does not have.  And a steady-state error suggests a gain that could be
> raised until it goes away — the band is a *region* of equilibria, and inside
> it the loop is not converging slowly, it is not moving.
