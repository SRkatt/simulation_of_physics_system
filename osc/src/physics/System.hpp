#pragma once
#include "physics/State.hpp"
#include "physics/Parameters.hpp"
#include <tuple>
#include <vector>

struct System {
    static std::tuple<double,double,double> computeEnergies(const State& s, const Parameters& p);
    static std::vector<double> perMassEnergies(const State& s, const Parameters& p);
};
