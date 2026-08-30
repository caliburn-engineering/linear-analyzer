// src/panels/model_panel.h
#pragma once

#include "../app_state.h"
#include <vector>

namespace caliburn {

void drawModelPanel(AppState& state, const std::vector<ModelEntry>& presets);
void extractTFFromSS(AppState& state);
void applyTFToSS(AppState& state);

}  // namespace caliburn
