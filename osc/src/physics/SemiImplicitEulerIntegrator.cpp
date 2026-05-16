#include "physics/SemiImplicitEulerIntegrator.hpp"

void SemiImplicitEulerIntegrator::step(State& state, const Parameters& params,
                                       const IForceModel& model, double dt) {
    auto a = model.accelerations(state, params);
    for (std::size_t i=0;i<state.size();++i) {
        state.v[i] += a[i]*dt;
        state.x[i] += state.v[i]*dt;
    }
}
