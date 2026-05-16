#include "rendering/ViewportRenderer.hpp"
#include "physics/Simulator.hpp"
#include <algorithm>
#include <cmath>

void ViewportRenderer::updateView(sf::View& view, const Simulator& sim, const sf::Vector2u& winSize,
                                  const sf::Vector2f& panOffset, float zoom) {
    const auto& p = sim.getParameters();
    std::size_t n = p.count();
    if (n == 0) return;

    float totalLength = 0.0f;
    for (std::size_t i = 0; i < n; ++i) totalLength += static_cast<float>(p.restLengths[i]);
    float maxInitial = (n == 0) ? 0.0f : static_cast<float>(p.initialX[n - 1]);
    float maxExtent = std::max(totalLength, maxInitial);

    float pad = std::max(1.0f, maxExtent * 0.2f);
    float w = (maxExtent + 2.0f * pad) / zoom;
    float h = w * (static_cast<float>(winSize.y) / winSize.x);
    view.setSize(w, h);
    view.setCenter(maxExtent * 0.5f + panOffset.x, panOffset.y);
}

void ViewportRenderer::draw(sf::RenderTarget& target, const Simulator& sim) {
    const auto& s = sim.getState();
    const auto& p = sim.getParameters();
    std::size_t n = s.size();
    if (n == 0) return;

    // Wall
    sf::RectangleShape wall({0.05f, 0.5f});
    wall.setOrigin(0.025f, 0.25f);
    wall.setPosition(0.0f, 0.0f);
    wall.setFillColor(sf::Color(120, 120, 120));
    target.draw(wall);

    // Springs
    for (std::size_t i = 0; i < n; ++i) {
        float x0 = (i == 0) ? 0.0f : static_cast<float>(s.x[i - 1]);
        float x1 = static_cast<float>(s.x[i]);
        float ext = (x1 - x0) - static_cast<float>(p.restLengths[i]);
        float strain = ext / static_cast<float>(p.restLengths[i]);
        sf::Color c(200, 200, 200);
        if (strain > 0.1f) c = sf::Color(220, 100, 100);
        else if (strain < -0.1f) c = sf::Color(100, 100, 220);

        sf::Vertex line[] = {
            sf::Vertex(sf::Vector2f(x0, 0.0f), c),
            sf::Vertex(sf::Vector2f(x1, 0.0f), c)
        };
        target.draw(line, 2, sf::Lines);

        int segs = 8; float amp = 0.03f;
        std::vector<sf::Vertex> zig;
        zig.reserve(segs + 1);
        for (int j = 0; j <= segs; ++j) {
            float t = j / float(segs);
            float xx = x0 + t * (x1 - x0);
            float yy = ((j % 2 == 0) ? 0.0f : ((j % 4 == 1) ? amp : -amp));
            zig.emplace_back(sf::Vector2f(xx, yy), c);
        }
        target.draw(zig.data(), zig.size(), sf::LineStrip);
    }

    // Masses
    for (std::size_t i = 0; i < n; ++i) {
        float r = MASS_RADIUS;
        sf::CircleShape circle(r);
        circle.setOrigin(r, r);
        circle.setPosition(static_cast<float>(s.x[i]), 0.0f);
        circle.setFillColor(sf::Color(230, 180, 80));
        circle.setOutlineColor(sf::Color(100, 80, 40));
        circle.setOutlineThickness(0.01f);
        target.draw(circle);
    }
}

std::optional<std::size_t> ViewportRenderer::hitTest(const sf::Vector2f& worldPos, const Simulator& sim) const {
    const auto& s = sim.getState();
    const auto& p = sim.getParameters();
    for (std::size_t i = s.size(); i-- > 0; ) {
        float r = MASS_RADIUS;
        float dx = worldPos.x - static_cast<float>(s.x[i]);
        float dy = worldPos.y - 0.0f;
        if (dx * dx + dy * dy <= r * r) return i;
    }
    return std::nullopt;
}
