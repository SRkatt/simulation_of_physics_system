#include "data/Preset.hpp"

std::vector<PresetEntry> getPresets() {
    std::vector<PresetEntry> out;
    {
        Parameters p; p.resize(2);
        p.initialX = {0.5, 1.0};
        p.setInitialFromCurrent(p.initialX);
        out.push_back({"Two Equal Masses", p});
    }
    {
        Parameters p; p.resize(2);
        p.masses[1] = 5.0;
        p.initialX = {0.5, 1.0};
        p.setInitialFromCurrent(p.initialX);
        out.push_back({"Heavy End", p});
    }
    {
        Parameters p; p.resize(3);
        p.springs[1] = 1.0;
        p.initialX = {0.3, 0.6, 0.9};
        p.setInitialFromCurrent(p.initialX);
        out.push_back({"Weak Middle Spring", p});
    }
    return out;
}
