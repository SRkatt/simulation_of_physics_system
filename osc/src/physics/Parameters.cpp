#include "physics/Parameters.hpp"

void Parameters::validate() const {
    std::size_t n = count();
    if (springs.size() != n || damping.size() != n ||
        restLengths.size() != n || initialX.size() != n) {
        throw std::runtime_error("Parameter vector size mismatch");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (masses[i] < 0.01) throw std::runtime_error("Mass too small (minimum 0.01 kg)");
        if (springs[i] <= 0.0) throw std::runtime_error("Spring constant must be positive");
        if (damping[i] < 0.0) throw std::runtime_error("Damping cannot be negative");
        if (restLengths[i] <= 0.0) throw std::runtime_error("Rest length must be positive");
    }
}

void Parameters::resize(std::size_t n, double defaultMass, double defaultSpring,
                        double defaultDamping, double defaultRestLength) {
    masses.assign(n, defaultMass);
    springs.assign(n, defaultSpring);
    damping.assign(n, defaultDamping);
    restLengths.assign(n, defaultRestLength);
    initialX.assign(n, 0.0);
}

void Parameters::setInitialFromCurrent(const std::vector<double>& pos) {
    if (pos.size() != count()) throw std::runtime_error("Initial position size mismatch");
    initialX = pos;
}
