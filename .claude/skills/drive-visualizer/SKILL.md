---
name: drive-visualizer
description: Launch, focus, click and screenshot the linear-analyzer visualizer GUI on Hyprland/Wayland. Use when verifying a UI change in the app, taking screenshots of panels, or driving the app to a specific state.
---

# Driving the linear-analyzer visualizer

Desktop GLFW/OpenGL app. On this machine (Hyprland + Wayland) it is driven with
`hyprctl` for focus, `ydotool` for input, and `grim` for screenshots.

## Launch

Run from the repo root — the app reads `imgui.ini` and `vendor/fonts/` relative to cwd.

```bash
cd projects/linear-analyzer
setsid ./build/visualizer >/tmp/viz.log 2>&1 &
sleep 4
```

## Find the window

```bash
A=$(hyprctl clients -j | python3 -c "
import json,sys
for c in json.load(sys.stdin):
    if 'Linear System' in c.get('title',''):
        print(c['address'])")
```

Geometry (`at` and `size`) comes from the same JSON — you need it to convert
screenshot pixel coordinates to screen coordinates.

## Focus it — REQUIRED before any click

**This Hyprland build uses a Lua dispatch API.** The old string form
(`hyprctl dispatch focuswindow address:0x...`) fails with a Lua parse error, and
so does `hl.dsp.window.focus` — there is no `focus` under `hl.dsp.window`.
The dispatcher is `hl.dsp.focus`, and it takes a **table**:

```bash
hyprctl dispatch "hl.dsp.focus({window='address:$A'})"
```

Without this, `ydotool` clicks go to whatever window currently has focus —
usually the terminal — and the app silently ignores them.

Discover any other dispatcher the same way, rather than guessing:

```bash
hyprctl repl 'local t={} for k,v in pairs(hl.dsp) do t[#t+1]=k end table.sort(t) return table.concat(t,", ")'
hyprctl dispatch "hl.dsp.focus({})"   # error message lists the accepted keys
```

`hyprctl repl <lua>` evaluates and prints; `hyprctl dispatch <lua>` wraps the
argument in `hl.dispatch(...)`.

## Click

`ydotoold` must be running (`pgrep ydotoold`). Coordinates are **screen**
coordinates: window `at` + pixel offset within the screenshot.

```bash
ydotool mousemove -a -x <X> -y <Y>; sleep 0.3; ydotool click 0xC0
```

`0xC0` is press+release of the left button. Re-focus the window before each
click burst if you ran anything that may have stolen focus.

## Screenshot

```bash
grim -g "<x>,<y> <w>x<h>" /tmp/shot.png     # region
```

Crop tightly to the panel you care about — a full 1896x1030 frame is ~1.9 MB and
mostly plots you are not checking. Then read the PNG.

## Gotchas

- **A blank or unchanged screenshot after a click usually means focus, not a
  code bug.** Re-check `hyprctl activewindow -j` before concluding anything.
- The window may be tiled small. It resizes when it is the only tiled client.
- The bundled `NotoSans-Regular.ttf` has **no arrow or math-operator glyphs** —
  U+2190 (left arrow), U+2220, U+221E all render as a hollow box, even though
  `visualizer.cpp` lists those ranges in `glyph_ranges`. Greek renders fine.
  Use ASCII (`<-`, `->`) in UI strings; do not add arrows expecting them to draw.
- Watch for **C++ hex-escape greediness** in UI strings: `"\xcf\x84f"` is parsed
  as `\xcf` then `\x84f` (out of range), not `tau` + `f`. Split the literal:
  `"\xcf\x84" "f"`.

## Click discipline

**Screenshot the region first, locate the widget, then click.** The model panel
scrolls, so a coordinate that was a combo box a moment ago can be a matrix
slider now — and a stray click on a `SliderFloat` silently edits the plant.
If a value looks wrong afterwards, restart the app rather than trying to undo.

Scrolling: `ydotool mousemove -w -x 0 -y -12` scrolls **down**, positive `-y`
scrolls up. `ydotool mousemove -a` does **not** work on this machine — the
cursor does not move — which is why focus uses `hl.dsp.cursor.move`.

Combos near the bottom of the window flip their popup upward; move the panel so
the widget sits mid-height before clicking, or the entry you want will not be
where you predict.

Killing the app: `pkill -x visualizer`. **Not** `pkill -f build/visualizer` —
`-f` matches the agent's own shell command line and kills the session.

## Workspaces and scroll traps

- **The window may be on another workspace.** `grim` captures the *visible*
  output and `ydotool` clicks land on the *visible* workspace, so a stale
  window elsewhere means you screenshot someone else's app and click into it.
  `hl.dsp.focus({window=...})` switches workspace as a side effect, so focus
  first, then confirm with `hyprctl activewindow -j` before clicking.
- **ImPlot eats the scroll wheel to zoom.** Scrolling with the cursor over any
  plot zooms that plot instead of scrolling its panel, and there is no undo —
  restart the app. Scroll panels from a non-plot column (the model panel is
  safe at x ~170).
