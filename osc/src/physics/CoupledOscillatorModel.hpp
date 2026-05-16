#pragma once
#include "physics/IForceModel.hpp"

class CoupledOscillatorModel : public IForceModel {
public:
    std::vector<double> accelerations(const State& state, const Parameters& params) const override;
};
