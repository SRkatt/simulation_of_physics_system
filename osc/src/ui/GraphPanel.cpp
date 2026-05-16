#include "ui/GraphPanel.hpp"
#include <imgui.h>
#include <algorithm>
#include <cmath>

static ImU32 HSLToRGB(float h, float s, float l) {
    float c = (1.0f - std::abs(2.0f * l - 1.0f)) * s;
    float hh = h / 60.0f;
    float x = c * (1.0f - std::abs(std::fmod(hh, 2.0f) - 1.0f));
    float r, g, b;
    if      (hh < 1) { r = c; g = x; b = 0; }
    else if (hh < 2) { r = x; g = c; b = 0; }
    else if (hh < 3) { r = 0; g = c; b = x; }
    else if (hh < 4) { r = 0; g = x; b = c; }
    else if (hh < 5) { r = x; g = 0; b = c; }
    else             { r = c; g = 0; b = x; }
    float m = l - c * 0.5f;
    return IM_COL32(
        static_cast<int>((r + m) * 255),
        static_cast<int>((g + m) * 255),
        static_cast<int>((b + m) * 255), 255);
}

void GraphPanel::render(const GraphHistory& hist) {
    ImGui::Begin("Energy History");
    ImVec2 pos = ImGui::GetWindowPos();
    ImVec2 sz = ImGui::GetWindowSize();
    float padL = 50, padB = 25, padT = 20, padR = 10;
    ImVec2 tl = {pos.x + padL, pos.y + padT};
    ImVec2 br = {pos.x + sz.x - padR, pos.y + sz.y - padB};
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddRectFilled(tl, br, IM_COL32(20, 20, 20, 255));
    if (hist.time.size() < 2) { ImGui::End(); return; }

    double maxE = 0.0;
    for (std::size_t i = 0; i < hist.total.size(); ++i) maxE = std::max(maxE, hist.total[i]);
    for (const auto& buf : hist.perMass)
        for (std::size_t i = 0; i < buf.size(); ++i) maxE = std::max(maxE, buf[i]);
    if (maxE < 1e-6) maxE = 1.0;

    double t0 = hist.time[0];
    double t1 = hist.time[hist.time.size() - 1];
    double dt = t1 - t0; if (dt < 0.001) dt = 1.0;

    auto map = [&](double t, double e)->ImVec2 {
        float x = tl.x + static_cast<float>((t - t0) / dt) * (br.x - tl.x);
        float y = br.y - static_cast<float>(e / maxE) * (br.y - tl.y);
        return {x, y};
    };

    for (int i = 0; i <= 4; ++i) {
        float y = tl.y + i * (br.y - tl.y) / 4.0f;
        dl->AddLine({tl.x, y}, {br.x, y}, IM_COL32(50, 50, 50, 255));
    }

    auto drawSeries = [&](const CircularBuffer<double>& buf, ImU32 col) {
        for (std::size_t i = 1; i < buf.size(); ++i)
            dl->AddLine(map(hist.time[i - 1], buf[i - 1]), map(hist.time[i], buf[i]), col, 2.0f);
    };

    // Total / KE / PE
    drawSeries(hist.total, IM_COL32(255, 255, 255, 255));
    drawSeries(hist.ke, IM_COL32(255, 100, 100, 255));
    drawSeries(hist.pe, IM_COL32(100, 100, 255, 255));

    // Per-mass energies
    for (std::size_t i = 0; i < hist.perMass.size(); ++i) {
        ImU32 col = HSLToRGB(static_cast<float>(i) * (360.0f / std::max(1u, static_cast<unsigned>(hist.perMass.size()))), 0.6f, 0.5f);
        drawSeries(hist.perMass[i], col);
    }

    dl->AddText({tl.x, br.y + 2}, IM_COL32(200, 200, 200, 255), "Time (s)");
    dl->AddText({tl.x - 35, tl.y}, IM_COL32(200, 200, 200, 255), "Energy");

    // Legend
    ImGui::Text("White=Total  Red=KE  Blue=PE");
    for (std::size_t i = 0; i < hist.perMass.size(); ++i) {
        ImU32 col = HSLToRGB(static_cast<float>(i) * (360.0f / std::max(1u, static_cast<unsigned>(hist.perMass.size()))), 0.6f, 0.5f);
        ImGui::SameLine();
        ImGui::TextColored(ImColor(col), "Mass %zu", i + 1);
    }

    ImGui::End();
}
