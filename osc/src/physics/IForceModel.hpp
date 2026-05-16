#pragma once
#include "physics/State.hpp"
#include "physics/Parameters.hpp"
#include <vector>

class IForceModel {
public:
    virtual ~IForceModel() = default;
    virtual std::vector<double> accelerations(const State& state, const Parameters& params) const = 0;
};
