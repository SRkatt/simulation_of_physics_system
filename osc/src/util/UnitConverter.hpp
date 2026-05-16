#pragma once
#include <SFML/System/Vector2.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Graphics/View.hpp>

class UnitConverter {
public:
    static sf::Vector2f screenToWorld(const sf::RenderWindow& window, const sf::View& view, const sf::Vector2i& screenPos);
};
