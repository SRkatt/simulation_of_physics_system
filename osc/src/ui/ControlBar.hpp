#pragma once

#include "app/AppState.hpp"

class ControlBar {
public:
    void render(AppState state, bool& requestPlay, bool& requestPause,
                bool& requestReset, bool& requestStep, float& timeScale);
};
