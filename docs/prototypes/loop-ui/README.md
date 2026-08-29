# PROTOTYPE — controller UI for a variable-length loop list

Throwaway. Wayfinder ticket #6. This branch exists only so the prototype
survives as a primary source; nothing here belongs on `master`.

## Run it

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --target visualizer -j8
    ./build/visualizer

Flip variants with the yellow bar at the bottom of the screen, or the ← / →
arrow keys. Env overrides for headless capture: `PROTO_VARIANT=0|1|2`,
`PROTO_PRESET=<index>`, `PROTO_BENCH=1` (opens the slider bench).

Presets worth trying: **2 = Ball-Balancer** (2×2, seeded identity pairing is
dead — the swap demo), **4 = Quarter-Car** (1 input / 2 outputs — Channel Share
instead of RGA).

## The three variants

| | structure | primary affordance |
|---|---|---|
| **A — Ledger** | one table row per loop, gains inline as columns | edit any gain of any loop without selecting |
| **B — Master / detail** | compact loop list + full-width gain block for the selected loop | the spec's slider block, one loop at a time |
| **C — Pairing grid** | the p×m diagnostic matrix *is* the pairing editor | click a cell to pair / unpair |

## What is real

The diagnostics are computed from the live plant via `evalTransferFunction`:
RGA on the complex `G(jω)` of the paired square sub-matrix (#4), Channel Share
where that is undefined (#8), and the structurally-dead test swept across the
whole Bode grid (#8). The loop list is seeded identity and reseeded by a `(p,m)`
guard (#2).

Nothing is wired into the recompute path — no controller is built, no closed
loop is formed. This draws a controller; it does not compute one.

## Screenshots

Both panel widths matter: the docked Model Configuration panel is **354 px** in
the checked-in `imgui.ini`. `*-354px.png` is that real width; `*-660px.png` is
the panel widened to see whether the layout is width-bound.
