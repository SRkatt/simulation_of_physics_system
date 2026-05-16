#pragma once

#include <SFML/System/Vector2.hpp>
#include <optional>
#include <cstddef>

class Simulator;
enum class AppState;

class InputHandler {
public:
    void startDrag(std::size_t massIndex);
    void updateDrag(const sf::Vector2f& worldPos, Simulator& sim, AppState state);
    void endDrag(Simulator& sim);
    bool isDragging() const { return dragged_.has_value(); }
    std::optional<std::size_t> dragged() const { return dragged_; }
    void setThrowScale(double scale) { throwScale_ = scale; }
private:
    std::optional<std::size_t> dragged_;
    bool hasPrevPos_ = false;
    sf::Vector2f prevWorldPos_;
    double releaseVelocityX_ = 0.0;
    double throwScale_ = 15.0;
};
