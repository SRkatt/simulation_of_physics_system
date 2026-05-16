#include "physics/Simulator.hpp"

Simulator::Simulator(Parameters params, std::unique_ptr<IIntegrator> integrator)
    : params_(std::move(params)), integrator_(std::move(integrator))
{
    rebuild(params_);
}

void Simulator::step(double dt) {
    integrator_->step(state_, params_, model_, dt);
    // Wall collision: leftmost mass cannot pass the wall at x=0
    for (std::size_t i = 0; i < state_.size(); ++i) {
        double minX = MassMinX(i);
        if (state_.x[i] < minX) {
            state_.x[i] = minX;
            if (state_.v[i] < 0.0) state_.v[i] = -state_.v[i] * 0.3; // damped bounce
        }
    }
    time_ += dt;
}

double Simulator::MassMinX(std::size_t i) const {
    // Mass i sits to the right of mass i-1 (or the wall for i=0)
    double prevEdge = (i == 0) ? 0.0 : state_.x[i - 1];
    return prevEdge + 0.001; // tiny gap to prevent overlap
}

void Simulator::reset() {
    state_ = initialState_;
    time_ = 0.0;
}

void Simulator::rebuild(const Parameters& params) {
    params_ = params;
    params_.validate();
    state_.resize(params_.count());
    for (std::size_t i=0;i<state_.size();++i) {
        state_.x[i] = params_.initialX[i];
        state_.v[i] = 0.0;
    }
    initialState_ = state_;
    time_ = 0.0;
}

void Simulator::setPosition(std::size_t i, double x) {
    if (i < state_.size()) state_.x[i] = x;
}

void Simulator::setVelocity(std::size_t i, double v) {
    if (i < state_.size()) state_.v[i] = v;
}

void Simulator::setIntegrator(std::unique_ptr<IIntegrator> integrator) {
    integrator_ = std::move(integrator);
}
