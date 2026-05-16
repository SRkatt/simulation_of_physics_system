#pragma once

#include <SFML/Graphics.hpp>

class Simulator;

class ArrowRenderer {
public:
    void draw(sf::RenderTarget& target, const Simulator& sim, bool drawVel, bool drawForce);
};
