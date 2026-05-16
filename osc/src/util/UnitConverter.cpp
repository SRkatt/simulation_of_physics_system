#include "util/UnitConverter.hpp"

sf::Vector2f UnitConverter::screenToWorld(const sf::RenderWindow& window, const sf::View& view, const sf::Vector2i& screenPos) {
    return window.mapPixelToCoords(screenPos, view);
}
