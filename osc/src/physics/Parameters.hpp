#pragma once
#include <vector>
#include <stdexcept>

struct Parameters {
    std::vector<double> masses;      // kg, size K
    std::vector<double> springs;     // N/m, size K (spring[i] connects mass i-1 to mass i; spring[0] is wall->mass0)
    std::vector<double> damping;     // kg/s, size K
    std::vector<double> restLengths; // m, size K
    std::vector<double> initialX;    // m, size K

    std::size_t count() const { return masses.size(); }
    void validate() const;
    void resize(std::size_t n, double defaultMass = 1.0, double defaultSpring = 10.0,
                double defaultDamping = 0.1, double defaultRestLength = 0.5);
    void setInitialFromCurrent(const std::vector<double>& pos);
};
