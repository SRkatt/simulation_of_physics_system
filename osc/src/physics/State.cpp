#include "physics/State.hpp"
#include <stdexcept>

void State::resize(std::size_t n) {
    x.resize(n, 0.0);
    v.resize(n, 0.0);
}

State& State::operator+=(const State& rhs) {
    if (x.size() != rhs.x.size()) throw std::runtime_error("State size mismatch");
    for (std::size_t i = 0; i < x.size(); ++i) {
        x[i] += rhs.x[i];
        v[i] += rhs.v[i];
    }
    return *this;
}

State State::operator+(const State& rhs) const {
    State result = *this;
    result += rhs;
    return result;
}

State State::operator*(double scalar) const {
    State result;
    result.resize(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        result.x[i] = x[i] * scalar;
        result.v[i] = v[i] * scalar;
    }
    return result;
}
