# Shared dependency sourcing for qsp_sim_core and its consumers.
#
# Pulls SUNDIALS, Boost::serialization, and yaml-cpp onto the project
# with a single include(). Consumers that need anything beyond the
# qsp_sim_core static library (e.g. a separate profile_jacobian target
# that also links SUNDIALS directly) can reuse the same helper so the
# version pins stay aligned.
#
# Variables set after this module runs:
#   SUNDIALS_LIBRARIES        — linkable target(s) for CVODE (+ KLU if enabled)
#   yaml-cpp (imported target `yaml-cpp::yaml-cpp`)
#   Boost::serialization
#
# Options surfaced at project scope:
#   USE_KLU                   — opt-in sparse linear solver via SuiteSparse

include_guard(GLOBAL)
include(FetchContent)

# --- KLU opt-in ------------------------------------------------------
option(USE_KLU "Enable KLU sparse linear solver via SuiteSparse (WIP)" OFF)
if(USE_KLU)
    find_package(SuiteSparse QUIET COMPONENTS KLU)
    if(NOT SuiteSparse_FOUND)
        find_path(KLU_INCLUDE_DIR klu.h PATH_SUFFIXES suitesparse)
        find_library(KLU_LIBRARY klu)
        find_library(AMD_LIBRARY amd)
        find_library(COLAMD_LIBRARY colamd)
        find_library(BTF_LIBRARY btf)
        find_library(SUITESPARSECONFIG_LIBRARY suitesparseconfig)
        if(KLU_INCLUDE_DIR AND KLU_LIBRARY)
            set(SuiteSparse_FOUND TRUE)
            set(KLU_LIBRARIES ${KLU_LIBRARY} ${AMD_LIBRARY} ${COLAMD_LIBRARY}
                              ${BTF_LIBRARY} ${SUITESPARSECONFIG_LIBRARY})
        else()
            message(WARNING "USE_KLU=ON but SuiteSparse/KLU not found. "
                            "Falling back to dense linear solver.")
            set(USE_KLU OFF)
        endif()
    endif()
endif()

# --- SUNDIALS --------------------------------------------------------
set(SUNDIALS_DIR "" CACHE PATH "SUNDIALS installation")
if(NOT SUNDIALS_DIR AND DEFINED ENV{SUNDIALS_DIR})
    set(SUNDIALS_DIR "$ENV{SUNDIALS_DIR}")
endif()

set(_SUNDIALS_READY FALSE)
if(SUNDIALS_DIR)
    list(APPEND CMAKE_PREFIX_PATH "${SUNDIALS_DIR}")
    find_package(SUNDIALS 7.0 QUIET COMPONENTS cvode)
    if(SUNDIALS_FOUND)
        set(_SUNDIALS_READY TRUE)
        set(SUNDIALS_LIBRARIES SUNDIALS::cvode)
        if(USE_KLU)
            list(APPEND SUNDIALS_LIBRARIES SUNDIALS::sunlinsolklu)
        endif()
    endif()
endif()

if(NOT _SUNDIALS_READY)
    message(STATUS "Fetching SUNDIALS 7.6.0...")
    FetchContent_Declare(sundials
        GIT_REPOSITORY https://github.com/LLNL/sundials.git
        GIT_TAG v7.6.0
        GIT_SHALLOW ON
    )
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(SUNDIALS_BUILD_WITH_MONITORING OFF CACHE BOOL "" FORCE)
    set(EXAMPLES_ENABLE_C OFF CACHE BOOL "" FORCE)
    set(EXAMPLES_ENABLE_CXX OFF CACHE BOOL "" FORCE)
    if(USE_KLU)
        set(ENABLE_KLU ON CACHE BOOL "" FORCE)
    endif()
    FetchContent_MakeAvailable(sundials)
    set(SUNDIALS_LIBRARIES sundials_cvode)
    if(USE_KLU)
        list(APPEND SUNDIALS_LIBRARIES sundials_sunlinsolklu)
    endif()
endif()

# --- Boost.serialization --------------------------------------------
find_package(Boost 1.70 REQUIRED COMPONENTS serialization)

# --- yaml-cpp --------------------------------------------------------
find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND)
    message(STATUS "Fetching yaml-cpp 0.8.0...")
    FetchContent_Declare(yaml-cpp
        GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
        GIT_TAG 0.8.0
        GIT_SHALLOW ON
    )
    set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
    set(YAML_CPP_BUILD_CONTRIB OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(yaml-cpp)
endif()
