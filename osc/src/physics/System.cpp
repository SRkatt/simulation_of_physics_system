#include "physics/System.hpp"

std::tuple<double,double,double> System::computeEnergies(const State& s, const Parameters& p) {
    double ke=0.0, pe=0.0;
    std::size_t n = p.count();
    for (std::size_t i=0;i<n;++i) {
        ke += 0.5*p.masses[i]*s.v[i]*s.v[i];
        double leftPos = (i==0)?0.0:s.x[i-1];
        double d = (s.x[i]-leftPos)-p.restLengths[i];
        pe += 0.5*p.springs[i]*d*d;
    }
    return {ke, pe, ke+pe};
}

std::vector<double> System::perMassEnergies(const State& s, const Parameters& p) {
    std::size_t n = p.count();
    std::vector<double> e(n,0.0);
    for (std::size_t i=0;i<n;++i) {
        e[i] += 0.5*p.masses[i]*s.v[i]*s.v[i];
        double leftPos = (i==0)?0.0:s.x[i-1];
        double dL = (s.x[i]-leftPos)-p.restLengths[i];
        e[i] += 0.5*p.springs[i]*dL*dL;
        if (i+1<n) {
            double dR = (s.x[i+1]-s.x[i])-p.restLengths[i+1];
            e[i] += 0.5*p.springs[i+1]*dR*dR;
        }
    }
    return e;
}
