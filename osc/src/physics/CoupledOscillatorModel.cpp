#include "physics/CoupledOscillatorModel.hpp"
#include <stdexcept>

std::vector<double> CoupledOscillatorModel::accelerations(const State& state, const Parameters& params) const {
    std::size_t n = params.count();
    if (state.size()!=n) throw std::runtime_error("Size mismatch");
    std::vector<double> a(n,0.0);
    for (std::size_t i=0;i<n;++i) {
        double leftPos = (i==0)?0.0:state.x[i-1];
        double dLeft = (state.x[i]-leftPos)-params.restLengths[i];
        double fLeft = -params.springs[i]*dLeft;
        double fRight = 0.0;
        if (i+1<n) {
            double dRight = (state.x[i+1]-state.x[i])-params.restLengths[i+1];
            fRight = params.springs[i+1]*dRight;
        }
        a[i] = (fLeft + fRight - params.damping[i]*state.v[i]) / params.masses[i];
    }
    return a;
}
