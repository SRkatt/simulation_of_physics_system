#pragma once
#include "physics/IIntegrator.hpp"

class RK4Integrator : public IIntegrator {
public:
    void step(State& state, const Parameters& params,
              const IForceModel& model, double dt) override;
    const char* name() const override { return "RK4"; }
};
