#pragma once
#include "physics/State.hpp"
#include "physics/Parameters.hpp"
#include "physics/IIntegrator.hpp"
#include "physics/CoupledOscillatorModel.hpp"
#include <memory>

class Simulator {
public:
    Simulator(Parameters params, std::unique_ptr<IIntegrator> integrator);
    void step(double dt);
    void reset();
    void rebuild(const Parameters& params);
    const State& getState() const { return state_; }
    const Parameters& getParameters() const { return params_; }
    Parameters& getParameters() { return params_; }
    void setPosition(std::size_t i, double x);
    void setVelocity(std::size_t i, double v);
    void setIntegrator(std::unique_ptr<IIntegrator> integrator);
    double getTime() const { return time_; }
private:
    double MassMinX(std::size_t i) const;
    Parameters params_;
    State state_;
    State initialState_;
    std::unique_ptr<IIntegrator> integrator_;
    CoupledOscillatorModel model_;
    double time_ = 0.0;
};
