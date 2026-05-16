#pragma once

#include "physics/Parameters.hpp"
#include "app/AppState.hpp"
#include <vector>

class ConfigPanel {
public:
    void render(Parameters& params, AppState state,
                const std::vector<const char*>& integratorNames,
                int& selectedIntegrator,
                bool& requestPlay,
                bool& requestApplyReset);
};
