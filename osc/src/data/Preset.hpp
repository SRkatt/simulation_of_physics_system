#pragma once

#include "physics/Parameters.hpp"
#include <vector>

struct PresetEntry {
    const char* name;
    Parameters params;
};

std::vector<PresetEntry> getPresets();
