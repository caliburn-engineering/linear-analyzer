# Physical Parameter Tuning & Derivation Display

## Summary

Add physical parameter sliders and ODE-to-state-space derivation steps to each preset model in the linear analyzer. Users can tune physical quantities (mass, damping, stiffness, etc.) and see how the governing equations transform from ODEs through Laplace to state-space form.

## Data Model

### New Structs (model_library.h)

```cpp
struct PhysicalParam {
    std::string name;    // Display name ("Mass")
    std::string symbol;  // Short symbol ("m")
    std::string unit;    // Unit string ("kg")
    float value;         // Default/current value
    float min_val;       // Slider minimum
    float max_val;       // Slider maximum
    bool logarithmic;    // Use logarithmic slider scale
};

struct DerivationStep {
    std::string title;   // "Equations of Motion (ODE)"
    std::string content; // Multi-line explanation with equations
};
```

### ModelEntry Extension

Add three fields to `ModelEntry`:
- `std::vector<PhysicalParam> params` — physical parameters with slider ranges
- `std::vector<DerivationStep> derivation` — ordered derivation steps
- `std::function<LinearSystem(const std::vector<PhysicalParam>&)> builder` — reconstructs SS from params

### AppState Addition

- `std::vector<PhysicalParam> current_params` — mutable copy of active preset's params

## Physical Parameters Per Model

| Model | Parameters | Default Values |
|---|---|---|
| First-Order | K (gain), tau (time constant) | K=1, tau=1 s |
| Mass-Spring-Damper | m (mass), b (damping), k (stiffness) | m=1 kg, b=1.4 Ns/m, k=1 N/m |
| Ball-Balancer | g (gravity) | g=9.81 m/s^2 |
| Inverted Pendulum | M (cart mass), m (pend. mass), l (length), g (gravity) | M=1 kg, m=0.1 kg, l=0.5 m, g=9.81 m/s^2 |
| Quarter-Car | ms (sprung), mu (unsprung), ks (spring), kt (tire), bs (damper) | ms=300 kg, mu=50 kg, ks=20kN/m, kt=200kN/m, bs=1kNs/m |
| Double Mass-Spring | m1, m2, k1, k2, b1, b2 | m1=m2=1 kg, k1=k2=1 N/m, b1=b2=0.1 Ns/m |

## Derivation Steps Per Model

Each model has 3-4 steps following the progression: Physical System -> ODE -> Laplace/TF -> State-Space.

Content uses UTF-8 Greek letters (matching existing codebase: tau, omega, zeta, theta) and plain-text equation formatting suitable for ImGui::TextWrapped rendering.

## Model Panel Layout

New order within the Plant section:

1. Preset dropdown + description
2. **Physical Parameters** — CollapsingHeader (default open). One slider per param, labeled "name [unit]". Calls builder on change, updates SS matrices and TF params.
3. **Derivation** — CollapsingHeader (default closed). Each step is a TreeNode that expands independently.
4. Transfer Function — existing auto-detected section for 1st/2nd order SISO
5. State-Space Matrices — existing text fields + slider grids
6. Controller, Channel, Traces, Reset — unchanged

## Data Flow

```
Physical param slider changed
  -> builder(current_params) produces new LinearSystem
  -> update plant matrices + text buffers
  -> extractTFFromSS updates TF display
  -> needs_recompute = true

Preset selected
  -> copy preset.params to state.current_params
  -> set plant = preset.system

TF slider / Matrix slider changed
  -> updates SS directly (physical params unchanged, last action wins)

Reset All
  -> restore preset[0].params to current_params
```

No back-propagation from TF/matrix sliders to physical params. One-way flow: physical params are the source of truth when used.

## Files Changed

| File | Change |
|---|---|
| `src/analysis/model_library.h` | Add PhysicalParam, DerivationStep structs. Extend ModelEntry. Add `#include <functional>`. |
| `src/analysis/model_library.cpp` | Add params, derivation text, and builder function to all 6 presets. |
| `src/app_state.h` | Add `std::vector<PhysicalParam> current_params` to AppState. |
| `src/panels/model_panel.cpp` | Add Physical Parameters and Derivation UI sections. Update preset loading and Reset to handle current_params. |
