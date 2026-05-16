#pragma once

#include "physics/Parameters.hpp"
#include <string>

bool saveToFile(const std::string& path, const Parameters& p, const std::string& integratorName);
bool loadFromFile(const std::string& path, Parameters& outParams, std::string& outIntegratorName);
