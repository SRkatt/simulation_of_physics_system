#include "rendering/ArrowRenderer.hpp"
#include "physics/Simulator.hpp"
#include "physics/CoupledOscillatorModel.hpp"

void ArrowRenderer::draw(sf::RenderTarget& target, const Simulator& sim, bool drawVel, bool drawForce) {
    if (!drawVel && !drawForce) return;
    const auto& s = sim.getState();
    const auto& p = sim.getParameters();
    if (s.size() == 0) return;
    CoupledOscillatorModel model;
    std::vector<double> a = model.accelerations(s, p);
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (drawVel) {
            float scale = 0.1f;
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(static_cast<float>(s.x[i]), 0.0f), sf::Color::Green),
                sf::Vertex(sf::Vector2f(static_cast<float>(s.x[i] + s.v[i] * scale), 0.0f), sf::Color::Green)
            };
            target.draw(line, 2, sf::Lines);
        }
        if (drawForce) {
            float scale = 0.05f;
            float f = static_cast<float>(a[i] * p.masses[i] * scale);
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(static_cast<float>(s.x[i]), 0.0f), sf::Color::Red),
                sf::Vertex(sf::Vector2f(static_cast<float>(s.x[i] + f), 0.0f), sf::Color::Red)
            };
            target.draw(line, 2, sf::Lines);
        }
    }
}
