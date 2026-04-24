cmake_minimum_required(VERSION 3.18)
project(qsp_sim_core CXX)

# Standalone build entry point. In practice consumers never build this
# directory on its own — they call find_package(qsp_sim_core CONFIG)
# which pulls in this tree via FetchContent or adds the package dir to
# CMAKE_PREFIX_PATH from the Python helper. We still support a standalone
# build to smoke-test the static library.
#
# Dependency sourcing (SUNDIALS, Boost::serialization, yaml-cpp) happens
# via include(qsp_sim_core/QspSimCoreDeps). Consumers include the same
# helper from their own CMakeLists so the whole dep graph lives in one
# place.

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(QspSimCoreDeps)

add_library(qsp_sim_core STATIC
    src/CVODEBase.cpp
    src/ParamBase.cpp
    src/default_hooks.cpp
)

target_include_directories(qsp_sim_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(qsp_sim_core PUBLIC
    ${SUNDIALS_LIBRARIES}
    Boost::serialization
    yaml-cpp::yaml-cpp
)

target_compile_features(qsp_sim_core PUBLIC cxx_std_17)

if(USE_KLU)
    target_compile_definitions(qsp_sim_core PUBLIC USE_KLU)
endif()

# Expose the driver main as a variable so consumers can compile it
# alongside their generated ODE + model-init sources.
set(QSP_SIM_CORE_DRIVER_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/qsp_sim_main.cpp"
    CACHE INTERNAL "Path to the generic qsp_sim driver main")

add_library(qsp_sim_core::qsp_sim_core ALIAS qsp_sim_core)
