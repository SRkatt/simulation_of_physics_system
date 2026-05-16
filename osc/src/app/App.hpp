#pragma once

#include <SFML/Graphics.hpp>
#include <memory>
#include <vector>
#include "app/AppState.hpp"
#include "app/InputHandler.hpp"
#include "rendering/ViewportRenderer.hpp"
#include "ui/ConfigPanel.hpp"
#include "ui/ControlBar.hpp"
#include "ui/GraphPanel.hpp"
#include "physics/Parameters.hpp"
#include "physics/Simulator.hpp"

class App {
public:
    App();
    ~App();
    void run();
private:
    sf::RenderWindow window_;
    sf::Clock deltaClock_;
    sf::View view_;
    AppState state_ = AppState::Editing;
    InputHandler inputHandler_;
    ViewportRenderer viewportRenderer_;
    ConfigPanel configPanel_;
    ControlBar controlBar_;
    GraphHistory history_;
    GraphPanel graphPanel_;
    Parameters params_;
    std::unique_ptr<Simulator> simulator_;
    int selectedIntegrator_ = 0;
    std::vector<const char*> integratorNames_ = {"Semi-Implicit Euler", "RK4"};
    double accumulator_ = 0.0;
    double timeScale_ = 1.0;
    static constexpr double FIXED_DT = 0.001;
    static constexpr int MAX_SUBSTEPS = 200;
    bool isPanning_ = false;
    sf::Vector2f panStartMouse_;
    sf::Vector2f panOffset_ = {0.0f, 0.0f};
    float zoom_ = 1.0f;
    void buildDefaultParams();
    void applyIntegratorSelection();
    void handleEvents();
    void update(float dt);
    void draw();
};
