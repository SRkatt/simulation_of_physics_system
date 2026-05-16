#include "ui/ConfigPanel.hpp"
#include <imgui.h>

void ConfigPanel::render(Parameters& params, AppState state,
                         const std::vector<const char*>& integratorNames,
                         int& selectedIntegrator,
                         bool& requestPlay,
                         bool& requestApplyReset) {
    requestPlay = false;
    requestApplyReset = false;
    ImGui::Begin("Configuration");
    bool editing = (state == AppState::Editing);

    if (editing) {
        if (ImGui::Button("Add Mass")) {
            std::size_t n = params.count();
            params.resize(n + 1);
            if (n > 0) {
                params.initialX[n] = params.initialX[n - 1] + params.restLengths[n - 1];
                params.setInitialFromCurrent(params.initialX);
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Remove Last") && params.count() > 1) {
            std::size_t n = params.count() - 1;
            params.resize(n);
            params.setInitialFromCurrent(params.initialX);
        }
    }

    ImGui::Separator();
    for (std::size_t i = 0; i < params.count(); ++i) {
        ImGui::PushID(static_cast<int>(i));
        ImGui::Text("Mass %zu", i + 1);
        bool canEdit = editing || (state == AppState::Paused);
        float m = static_cast<float>(params.masses[i]);
        if (ImGui::SliderFloat("Mass (kg)", &m, 0.01f, 50.0f)) params.masses[i] = m;
        float c = static_cast<float>(params.damping[i]);
        if (ImGui::SliderFloat("Damping (Ns/m)", &c, 0.0f, 20.0f)) params.damping[i] = c;
        if (canEdit) {
            float x0 = static_cast<float>(params.initialX[i]);
            if (ImGui::SliderFloat("Initial X (m)", &x0, -5.0f, 5.0f)) params.initialX[i] = x0;
        }
        ImGui::PopID();
    }

    ImGui::Separator();
    ImGui::Text("Springs");
    for (std::size_t i = 0; i < params.count(); ++i) {
        ImGui::PushID(static_cast<int>(i) + 1000);
        float k = static_cast<float>(params.springs[i]);
        if (ImGui::SliderFloat("k (N/m)", &k, 0.1f, 500.0f)) params.springs[i] = k;
        float rl = static_cast<float>(params.restLengths[i]);
        if (ImGui::SliderFloat("Rest L (m)", &rl, 0.01f, 2.0f)) params.restLengths[i] = rl;
        ImGui::PopID();
    }

    ImGui::Separator();
    if (!integratorNames.empty()) {
        ImGui::Combo("Integrator", &selectedIntegrator, integratorNames.data(),
                     static_cast<int>(integratorNames.size()));
    }

    if (editing) {
        if (ImGui::Button("Apply & Reset")) requestApplyReset = true;
        ImGui::SameLine();
        if (ImGui::Button("Play")) requestPlay = true;
    }
    ImGui::End();
}
