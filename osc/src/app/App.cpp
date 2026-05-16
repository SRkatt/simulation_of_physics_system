#include "app/App.hpp"
#include "physics/System.hpp"
#include "physics/RK4Integrator.hpp"
#include "physics/SemiImplicitEulerIntegrator.hpp"
#include "util/UnitConverter.hpp"
#include "data/Preset.hpp"
#include "data/Serializer.hpp"
#include <algorithm>
#include <imgui.h>
#include <imgui-SFML.h>

App::App()
    : window_(sf::VideoMode(1280, 720), "Coupled Oscillator Sandbox")
{
    window_.setFramerateLimit(60);
    buildDefaultParams();
    simulator_ = std::make_unique<Simulator>(params_, std::make_unique<SemiImplicitEulerIntegrator>());
    if (!ImGui::SFML::Init(window_)) { /* ImGui init failed */ }
}

App::~App() { ImGui::SFML::Shutdown(); }

void App::buildDefaultParams() {
    params_.resize(2);
    params_.initialX[0] = 0.5;
    params_.initialX[1] = 1.0;
    params_.setInitialFromCurrent(params_.initialX);
}

void App::applyIntegratorSelection() {
    if (selectedIntegrator_ == 0)
        simulator_->setIntegrator(std::make_unique<SemiImplicitEulerIntegrator>());
    else
        simulator_->setIntegrator(std::make_unique<RK4Integrator>());
}

void App::run() {
    sf::Clock dtClock;
    while (window_.isOpen()) {
        float dt = dtClock.restart().asSeconds();
        handleEvents();
        update(dt);
        draw();
    }
}

void App::handleEvents() {
    // Ensure view_ reflects current pan/zoom before converting coordinates
    viewportRenderer_.updateView(view_, *simulator_, window_.getSize(), panOffset_, zoom_);

    sf::Event e;
    while (window_.pollEvent(e)) {
        ImGui::SFML::ProcessEvent(window_, e);
        if (e.type == sf::Event::Closed) window_.close();

        bool ioWantMouse = ImGui::GetIO().WantCaptureMouse;

        // Zoom
        if (!ioWantMouse && e.type == sf::Event::MouseWheelScrolled) {
            float factor = (e.mouseWheelScroll.delta > 0) ? 1.2f : (1.0f / 1.2f);
            zoom_ *= factor;
            zoom_ = std::clamp(zoom_, 0.1f, 10.0f);
        }

        // Pan start
        if (!ioWantMouse && e.type == sf::Event::MouseButtonPressed
            && e.mouseButton.button == sf::Mouse::Right) {
            isPanning_ = true;
            panStartMouse_ = {static_cast<float>(e.mouseButton.x), static_cast<float>(e.mouseButton.y)};
        }

        // Pan move
        if (e.type == sf::Event::MouseMoved && isPanning_) {
            sf::Vector2f cur = {static_cast<float>(e.mouseMove.x), static_cast<float>(e.mouseMove.y)};
            sf::Vector2f deltaScreen = cur - panStartMouse_;
            sf::Vector2f size = view_.getSize();
            sf::Vector2u winsz = window_.getSize();
            float worldPerPixelX = size.x / winsz.x;
            float worldPerPixelY = size.y / winsz.y;
            panOffset_.x -= deltaScreen.x * worldPerPixelX;
            panOffset_.y -= deltaScreen.y * worldPerPixelY;
            panStartMouse_ = cur;
        }

        // Pan end
        if (e.type == sf::Event::MouseButtonReleased && e.mouseButton.button == sf::Mouse::Right) {
            isPanning_ = false;
        }

        // Drag mass (left click) — always active regardless of ImGui state
        if (e.type == sf::Event::MouseButtonPressed && e.mouseButton.button == sf::Mouse::Left) {
            sf::Vector2f world = UnitConverter::screenToWorld(window_, view_,
                {e.mouseButton.x, e.mouseButton.y});
            auto hit = viewportRenderer_.hitTest(world, *simulator_);
            if (hit) inputHandler_.startDrag(*hit);
        }
        if (e.type == sf::Event::MouseMoved && inputHandler_.isDragging()) {
            sf::Vector2f world = UnitConverter::screenToWorld(window_, view_,
                {e.mouseMove.x, e.mouseMove.y});
            inputHandler_.updateDrag(world, *simulator_, state_);
        }
        if (e.type == sf::Event::MouseButtonReleased && e.mouseButton.button == sf::Mouse::Left) {
            if (inputHandler_.isDragging()) inputHandler_.endDrag(*simulator_);
        }

        if (e.type == sf::Event::KeyPressed) {
            if (e.key.code == sf::Keyboard::Space) {
                if (state_ == AppState::Running) state_ = AppState::Paused;
                else if (state_ == AppState::Paused) state_ = AppState::Running;
                else if (state_ == AppState::Editing) {
                    applyIntegratorSelection();
                    simulator_->rebuild(params_);
                    history_.clear();
                    state_ = AppState::Running;
                }
            }
            if (e.key.code == sf::Keyboard::R) {
                state_ = AppState::Editing;
                simulator_->rebuild(params_);
                history_.clear();
            }
        }
    }
}

void App::update(float realDt) {
    ImGui::SFML::Update(window_, sf::seconds(realDt));

    bool requestPlay = false, requestApplyReset = false, requestPause = false,
         requestReset = false, requestStep = false;

    float ts = static_cast<float>(timeScale_);
    configPanel_.render(params_, state_, integratorNames_, selectedIntegrator_,
                        requestPlay, requestApplyReset);
    controlBar_.render(state_, requestPlay, requestPause, requestReset, requestStep, ts);
    timeScale_ = static_cast<double>(ts);

    if (requestApplyReset) {
        applyIntegratorSelection();
        simulator_->rebuild(params_);
        history_.clear();
        state_ = AppState::Editing;
    }
    if (requestPlay) {
        if (state_ == AppState::Editing) {
            applyIntegratorSelection();
            simulator_->rebuild(params_);
            history_.clear();
        }
        state_ = AppState::Running;
    }
    if (requestPause) state_ = AppState::Paused;
    if (requestReset) {
        state_ = AppState::Editing;
        simulator_->rebuild(params_);
        history_.clear();
    }

    // Physics stepping
    if (state_ == AppState::Running) {
        accumulator_ += realDt * timeScale_;
        int steps = 0;
        while (accumulator_ >= FIXED_DT && steps < MAX_SUBSTEPS) {
            simulator_->step(FIXED_DT);
            auto [ke, pe, tot] = System::computeEnergies(simulator_->getState(), simulator_->getParameters());
            auto perMass = System::perMassEnergies(simulator_->getState(), simulator_->getParameters());
            history_.record(simulator_->getTime(), ke, pe, tot, perMass);
            accumulator_ -= FIXED_DT;
            ++steps;
        }
    }
    if (requestStep && state_ == AppState::Paused) {
        simulator_->step(FIXED_DT);
        auto [ke, pe, tot] = System::computeEnergies(simulator_->getState(), simulator_->getParameters());
        auto perMass = System::perMassEnergies(simulator_->getState(), simulator_->getParameters());
        history_.record(simulator_->getTime(), ke, pe, tot, perMass);
    }

    // Live parameter sync during Running/Paused
    if (state_ != AppState::Editing) {
        auto& sp = simulator_->getParameters();
        sp.masses = params_.masses;
        sp.springs = params_.springs;
        sp.damping = params_.damping;
    }
}

void App::draw() {
    window_.clear(sf::Color(30, 30, 30));
    viewportRenderer_.updateView(view_, *simulator_, window_.getSize(), panOffset_, zoom_);
    window_.setView(view_);
    viewportRenderer_.draw(window_, *simulator_);
    window_.setView(window_.getDefaultView());
    graphPanel_.render(history_);
    ImGui::SFML::Render(window_);
    window_.display();
}
