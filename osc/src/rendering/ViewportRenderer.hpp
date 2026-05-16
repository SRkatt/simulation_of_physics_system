#pragma once

#include <SFML/Graphics.hpp>
#include <optional>
#include <cstddef>

class Simulator;

class ViewportRenderer {
public:
    void updateView(sf::View& view, const Simulator& sim, const sf::Vector2u& winSize,
                    const sf::Vector2f& panOffset = {0.0f, 0.0f}, float zoom = 1.0f);
    void draw(sf::RenderTarget& target, const Simulator& sim);
    std::optional<std::size_t> hitTest(const sf::Vector2f& worldPos, const Simulator& sim) const;
    static constexpr float MASS_RADIUS = 0.08f;
};
