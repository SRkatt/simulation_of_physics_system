#include "physics/State.hpp"
#include "physics/Parameters.hpp"
#include "physics/System.hpp"
#include "physics/Simulator.hpp"
#include "physics/RK4Integrator.hpp"
#include "physics/SemiImplicitEulerIntegrator.hpp"
#include <iostream>
#include <cmath>

bool approx(double a, double b, double tol=1e-4) { return std::abs(a-b)<tol; }

bool test_rk4_energy_conservation() {
    Parameters p; p.resize(1,1.0,10.0,0.0,1.0);
    p.initialX[0]=0.5; p.setInitialFromCurrent(p.initialX);
    Simulator sim(p, std::make_unique<RK4Integrator>());
    auto [ke0,pe0,e0] = System::computeEnergies(sim.getState(), sim.getParameters());
    for (int i=0;i<10000;++i) sim.step(0.001);
    auto [ke1,pe1,e1] = System::computeEnergies(sim.getState(), sim.getParameters());
    return std::abs(e1-e0)/e0 < 0.001;
}

bool test_semieuler_stable() {
    Parameters p; p.resize(1,1.0,10.0,0.0,1.0);
    p.initialX[0]=0.5; p.setInitialFromCurrent(p.initialX);
    Simulator sim(p, std::make_unique<SemiImplicitEulerIntegrator>());
    for (int i=0;i<5000;++i) sim.step(0.001);
    auto [ke,pe,e] = System::computeEnergies(sim.getState(), sim.getParameters());
    return e < 2.0 && std::abs(sim.getState().x[0]) < 2.0;
}

bool test_damping_decay() {
    Parameters p; p.resize(1,1.0,10.0,2.0,1.0);
    p.initialX[0]=1.0; p.setInitialFromCurrent(p.initialX);
    Simulator sim(p, std::make_unique<RK4Integrator>());
    double prev = std::get<2>(System::computeEnergies(sim.getState(), sim.getParameters()));
    for (int i=0;i<500;++i) {
        sim.step(0.001);
        double e = std::get<2>(System::computeEnergies(sim.getState(), sim.getParameters()));
        if (e > prev + 1e-8) return false;
        prev = e;
    }
    return true;
}

bool test_params_validation() {
    Parameters p; p.resize(1);
    try { p.masses[0]=0.0; p.validate(); return false; }
    catch (...) { return true; }
}

int main() {
    int f=0;
    if (!test_rk4_energy_conservation()) { ++f; std::cerr<<"FAIL: RK4 energy conservation\n"; }
    if (!test_semieuler_stable()) { ++f; std::cerr<<"FAIL: SemiImplicitEuler stable\n"; }
    if (!test_damping_decay()) { ++f; std::cerr<<"FAIL: Damping energy decay\n"; }
    if (!test_params_validation()) { ++f; std::cerr<<"FAIL: Parameter validation\n"; }
    if (f==0) std::cout<<"All tests passed.\n";
    else std::cout<<f<<" tests failed.\n";
    return f;
}
