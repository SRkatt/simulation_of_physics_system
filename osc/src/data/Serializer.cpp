#include "data/Serializer.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

bool saveToFile(const std::string& path, const Parameters& p, const std::string& integratorName) {
    json j;
    j["version"] = 1;
    j["params"]["masses"] = p.masses;
    j["params"]["springs"] = p.springs;
    j["params"]["damping"] = p.damping;
    j["params"]["restLengths"] = p.restLengths;
    j["params"]["initialX"] = p.initialX;
    j["settings"]["integrator"] = integratorName;
    std::ofstream f(path);
    if (!f) return false;
    f << j.dump(2);
    return true;
}

bool loadFromFile(const std::string& path, Parameters& outP, std::string& outName) {
    std::ifstream f(path);
    if (!f) return false;
    json j;
    try { f >> j; } catch (...) { return false; }
    if (!j.contains("version") || j["version"] != 1) return false;
    try {
        outP.masses = j["params"]["masses"].get<std::vector<double>>();
        outP.springs = j["params"]["springs"].get<std::vector<double>>();
        outP.damping = j["params"]["damping"].get<std::vector<double>>();
        outP.restLengths = j["params"]["restLengths"].get<std::vector<double>>();
        outP.initialX = j["params"]["initialX"].get<std::vector<double>>();
        outP.setInitialFromCurrent(outP.initialX);
        outName = j.value("/settings/integrator"_json_pointer, "Semi-Implicit Euler");
        outP.validate();
        return true;
    } catch (...) { return false; }
}
