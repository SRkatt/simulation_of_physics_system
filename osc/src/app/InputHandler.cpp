#include "app/InputHandler.hpp"
#include "physics/Simulator.hpp"
#include "app/AppState.hpp"

void InputHandler::startDrag(std::size_t massIndex) {
    dragged_ = massIndex;
    hasPrevPos_ = false;
    releaseVelocityX_ = 0.0;
}

void InputHandler::updateDrag(const sf::Vector2f& worldPos, Simulator& sim, AppState state) {
    if (!dragged_) return;

    if (hasPrevPos_) {
        releaseVelocityX_ = (worldPos.x - prevWorldPos_.x) * throwScale_;
    } else {
        hasPrevPos_ = true;
    }
    prevWorldPos_ = worldPos;

    if (state == AppState::Editing || state == AppState::Paused) {
        sim.setPosition(*dragged_, worldPos.x);
        sim.setVelocity(*dragged_, 0.0);
    } else if (state == AppState::Running) {
        sim.setPosition(*dragged_, worldPos.x);
    }
}

void InputHandler::endDrag(Simulator& sim) {
    if (dragged_ && *dragged_ < sim.getParameters().count()) {
        if (std::abs(releaseVelocityX_) > 0.01) {
            sim.setVelocity(*dragged_, releaseVelocityX_);
        }
        sim.getParameters().initialX[*dragged_] = sim.getState().x[*dragged_];
    }
    dragged_ = std::nullopt;
    hasPrevPos_ = false;
}
