#include "physics/RK4Integrator.hpp"
#include <vector>

static State derivative(const State& s, const std::vector<double>& a) {
    State d; d.resize(s.size());
    for (std::size_t i=0;i<s.size();++i){ d.x[i]=s.v[i]; d.v[i]=a[i]; }
    return d;
}

void RK4Integrator::step(State& state, const Parameters& params,
                         const IForceModel& model, double dt) {
    std::size_t n = state.size();
    auto a1 = model.accelerations(state, params);
    State k1 = derivative(state, a1);

    State s2 = state + k1*(dt*0.5);
    auto a2 = model.accelerations(s2, params);
    State k2 = derivative(s2, a2);

    State s3 = state + k2*(dt*0.5);
    auto a3 = model.accelerations(s3, params);
    State k3 = derivative(s3, a3);

    State s4 = state + k3*dt;
    auto a4 = model.accelerations(s4, params);
    State k4 = derivative(s4, a4);

    for (std::size_t i=0;i<n;++i) {
        state.x[i] += (k1.x[i] + 2.0*k2.x[i] + 2.0*k3.x[i] + k4.x[i])*dt/6.0;
        state.v[i] += (k1.v[i] + 2.0*k2.v[i] + 2.0*k3.v[i] + k4.v[i])*dt/6.0;
    }
}
