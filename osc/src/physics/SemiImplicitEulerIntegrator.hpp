#pragma once
#include "physics/IIntegrator.hpp"

class SemiImplicitEulerIntegrator : public IIntegrator {
public:
    void step(State& state, const Parameters& params,
              const IForceModel& model, double dt) override;
    const char* name() const override { return "Semi-Implicit Euler"; }
};
