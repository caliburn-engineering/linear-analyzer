// src/panels/loop_ui_prototype.h
//
// PROTOTYPE — THROWAWAY. Wayfinder ticket #6.
//
// Question: what does the controller section of the model panel look like once
// the loop list is variable-length?
//
// Three structurally different variants of the controller section, rendered in
// place of the real one inside the existing "Model Configuration" panel, and
// switched from a floating bar at the bottom of the screen (← / → or the arrow
// keys). Not wired into the recompute path — nothing here computes a
// controller; it only draws one.
//
// Delete this file and its two call sites in model_panel.cpp when #6 resolves.
#pragma once

#include "../app_state.h"

namespace caliburn::proto {

// Draws the prototype controller section (call in place of the real one).
void drawControllerSection(AppState& state);

// Draws the floating variant switcher. Call once per frame, after the panel.
void drawSwitcher();

}  // namespace caliburn::proto
