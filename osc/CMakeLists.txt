cmake_minimum_required(VERSION 3.16)
project(OscillatorSandbox VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

FetchContent_Declare(SFML
    GIT_REPOSITORY https://github.com/SFML/SFML.git
    GIT_TAG 2.6.1
    GIT_SHALLOW TRUE
)
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG v1.89.9
    GIT_SHALLOW TRUE
)
FetchContent_Declare(imgui-sfml
    GIT_REPOSITORY https://github.com/SFML/imgui-sfml.git
    GIT_TAG v2.6
    GIT_SHALLOW TRUE
)
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.3
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(imgui json SFML)
set(IMGUI_DIR ${imgui_SOURCE_DIR})
set(IMGUI_SFML_FIND_SFML OFF)
FetchContent_MakeAvailable(imgui-sfml)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)
add_executable(oscillator_sandbox ${SOURCES})
target_include_directories(oscillator_sandbox PRIVATE src)
target_link_libraries(oscillator_sandbox PRIVATE
    sfml-graphics sfml-window sfml-system
    ImGui-SFML::ImGui-SFML
    nlohmann_json::nlohmann_json
)

file(GLOB_RECURSE TEST_SOURCES CONFIGURE_DEPENDS tests/*.cpp src/physics/*.cpp)
add_executable(test_physics ${TEST_SOURCES})
target_include_directories(test_physics PRIVATE src)
