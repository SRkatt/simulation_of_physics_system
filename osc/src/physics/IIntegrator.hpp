#pragma once
#include "physics/State.hpp"
#include "physics/Parameters.hpp"
#include "physics/IForceModel.hpp"

class IIntegrator {
public:
    virtual ~IIntegrator() = default;
    virtual void step(State& state, const Parameters& params,
                      const IForceModel& model, double dt) = 0;
    virtual const char* name() const = 0;
};
