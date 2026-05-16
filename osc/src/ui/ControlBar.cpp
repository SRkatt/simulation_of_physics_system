#include "ui/ControlBar.hpp"
#include <imgui.h>

void ControlBar::render(AppState state, bool& requestPlay, bool& requestPause,
                        bool& requestReset, bool& requestStep, float& timeScale) {
    requestPlay = requestPause = requestReset = requestStep = false;
    if (ImGui::BeginMainMenuBar()) {
        if ((state == AppState::Editing || state == AppState::Paused) && ImGui::Button("Play"))
            requestPlay = true;
        if ((state == AppState::Running || state == AppState::Paused) && ImGui::Button("Reset"))
            requestReset = true;
        if (state == AppState::Running && ImGui::Button("Pause"))
            requestPause = true;
        if (state == AppState::Paused && ImGui::Button("Step"))
            requestStep = true;
        ImGui::Text("| State: %s", state == AppState::Editing ? "Editing" :
                     (state == AppState::Running ? "Running" : "Paused"));
        ImGui::SliderFloat("Speed", &timeScale, 0.05f, 5.0f, "%.2fx");
        ImGui::EndMainMenuBar();
    }
}
