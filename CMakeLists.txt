# projects/linear-analyzer/CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(linear-analyzer LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(FetchContent)

# --- Eigen ---
FetchContent_Declare(eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.0
    GIT_SHALLOW    TRUE
)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen)

# --- GLFW ---
FetchContent_Declare(glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
    GIT_SHALLOW    TRUE
)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(glfw)

# --- GLAD (vendored) ---
add_library(glad STATIC vendor/glad/src/glad.c)
target_include_directories(glad PUBLIC vendor/glad/include)

# --- ImGui ---
FetchContent_Declare(imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        docking
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(imgui)
add_library(imgui STATIC
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_demo.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)
target_include_directories(imgui PUBLIC
    ${imgui_SOURCE_DIR}
    ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PUBLIC glfw glad)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLAD)

# --- ImPlot ---
FetchContent_Declare(implot
    GIT_REPOSITORY https://github.com/epezent/implot.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(implot)
add_library(implot STATIC
    ${implot_SOURCE_DIR}/implot.cpp
    ${implot_SOURCE_DIR}/implot_items.cpp
)
target_include_directories(implot PUBLIC ${implot_SOURCE_DIR})
target_link_libraries(implot PUBLIC imgui)

# --- Analysis library (sources added incrementally) ---
add_library(analysis_lib STATIC
    src/rk4.cpp
    src/analysis/model_library.cpp
    src/analysis/system_properties.cpp
    src/analysis/frequency_response.cpp
    src/analysis/pole_zero.cpp
)
target_include_directories(analysis_lib PUBLIC src)
target_link_libraries(analysis_lib PUBLIC Eigen3::Eigen)

# --- Tests ---
enable_testing()

add_executable(test_model_library tests/test_model_library.cpp)
target_link_libraries(test_model_library PRIVATE analysis_lib)
target_include_directories(test_model_library PRIVATE tests)
add_test(NAME test_model_library COMMAND test_model_library)

add_executable(test_system_properties tests/test_system_properties.cpp)
target_link_libraries(test_system_properties PRIVATE analysis_lib)
target_include_directories(test_system_properties PRIVATE tests)
add_test(NAME test_system_properties COMMAND test_system_properties)

add_executable(test_frequency_response tests/test_frequency_response.cpp)
target_link_libraries(test_frequency_response PRIVATE analysis_lib)
target_include_directories(test_frequency_response PRIVATE tests)
add_test(NAME test_frequency_response COMMAND test_frequency_response)

add_executable(test_pole_zero tests/test_pole_zero.cpp)
target_link_libraries(test_pole_zero PRIVATE analysis_lib)
target_include_directories(test_pole_zero PRIVATE tests)
add_test(NAME test_pole_zero COMMAND test_pole_zero)
