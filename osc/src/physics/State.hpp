#pragma once
#include <vector>
#include <cstddef>

struct State {
    std::vector<double> x; // positions
    std::vector<double> v; // velocities

    std::size_t size() const { return x.size(); }
    void resize(std::size_t n);
    State& operator+=(const State& rhs);
    State operator+(const State& rhs) const;
    State operator*(double scalar) const;
};
