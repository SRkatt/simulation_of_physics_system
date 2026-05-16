#pragma once

#include "util/CircularBuffer.hpp"
#include <vector>

struct GraphHistory {
    CircularBuffer<double> time{2000}, ke{2000}, pe{2000}, total{2000};
    std::vector<CircularBuffer<double>> perMass;
    void record(double t, double k, double p, double e, const std::vector<double>& massEnergies) {
        time.push(t); ke.push(k); pe.push(p); total.push(e);
        for (std::size_t i = 0; i < massEnergies.size(); ++i) {
            if (i >= perMass.size()) perMass.emplace_back(2000);
            perMass[i].push(massEnergies[i]);
        }
    }
    void clear() {
        time.clear(); ke.clear(); pe.clear(); total.clear(); perMass.clear();
    }
};

class GraphPanel {
public:
    void render(const GraphHistory& hist);
};
