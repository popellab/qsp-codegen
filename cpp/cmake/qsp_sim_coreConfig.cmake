# qsp_sim_coreConfig.cmake — consumed by `find_package(qsp_sim_core CONFIG)`.
#
# Does not build anything on its own — `add_subdirectory` compiles the
# static library in-tree inside the consumer's build. We go through
# add_subdirectory instead of shipping a prebuilt install tree because
# the library must be compiled against the consumer's SUNDIALS /
# yaml-cpp configuration (ABI match) and because the qsp-codegen wheel
# ships sources only.
#
# After find_package completes, the consumer has:
#   target qsp_sim_core::qsp_sim_core      — static lib + transitive deps
#   variable QSP_SIM_CORE_DRIVER_SOURCE     — absolute path to qsp_sim_main.cpp
# and can include(QspSimCoreDeps) itself to reuse the shared dep graph.

get_filename_component(_qsp_sim_core_this_dir "${CMAKE_CURRENT_LIST_FILE}" DIRECTORY)
set(_qsp_sim_core_root "${_qsp_sim_core_this_dir}/..")
get_filename_component(_qsp_sim_core_root "${_qsp_sim_core_root}" ABSOLUTE)

list(APPEND CMAKE_MODULE_PATH "${_qsp_sim_core_root}/cmake")

# Canonicalize embedded build paths so the consumer's binary is byte-
# deterministic across checkouts/worktrees of the same source. Without this,
# absolute paths in __FILE__ macros and DWARF debug info embed the consumer's
# project / build directory; multiplied across compilation units this
# cascades through PIE relocations to ~65% of the binary, which forces any
# build-byte-keyed cache (e.g. the SBI pool-hash) to invalidate per worktree.
#
# Set via add_compile_options BEFORE qsp_sim_core's add_subdirectory below so
# the flags propagate into sundials/yaml-cpp FetchContent builds too.
# Opt out via QSP_SIM_CORE_NO_CANONICAL_PATHS=ON if you actually need the
# real paths embedded (debugger source resolution against an unrelocatable
# build tree, etc.).
option(QSP_SIM_CORE_NO_CANONICAL_PATHS
    "Disable -ffile-prefix-map / -no_uuid path canonicalization"
    OFF)
if(NOT QSP_SIM_CORE_NO_CANONICAL_PATHS)
    get_filename_component(_qsp_consumer_src "${CMAKE_SOURCE_DIR}" ABSOLUTE)
    get_filename_component(_qsp_consumer_bin "${CMAKE_BINARY_DIR}" ABSOLUTE)
    # Use CMAKE_{C,CXX}_FLAGS strings rather than add_compile_options so the
    # flags propagate into FetchContent sub-builds (sundials, yaml-cpp). Their
    # CMakeLists call `project()` which resets directory-scoped COMPILE_OPTIONS
    # but inherits the FLAGS strings.
    set(_qsp_canonical_flags
        "-ffile-prefix-map=${_qsp_consumer_src}=/qsp_project_src"
        "-ffile-prefix-map=${_qsp_consumer_bin}=/qsp_project_build"
        "-ffile-prefix-map=${_qsp_sim_core_root}=/qsp_sim_core")
    foreach(_flag IN LISTS _qsp_canonical_flags)
        if(NOT CMAKE_C_FLAGS MATCHES "${_flag}")
            string(APPEND CMAKE_C_FLAGS " ${_flag}")
        endif()
        if(NOT CMAKE_CXX_FLAGS MATCHES "${_flag}")
            string(APPEND CMAKE_CXX_FLAGS " ${_flag}")
        endif()
    endforeach()
    if(APPLE)
        # Strip the random Mach-O LC_UUID slot so two builds of the same
        # source don't differ in that 16-byte field.
        if(NOT CMAKE_EXE_LINKER_FLAGS MATCHES "no_uuid")
            string(APPEND CMAKE_EXE_LINKER_FLAGS " -Wl,-no_uuid")
        endif()
    endif()
endif()

if(NOT TARGET qsp_sim_core::qsp_sim_core)
    add_subdirectory("${_qsp_sim_core_root}"
                     "${CMAKE_BINARY_DIR}/_qsp_sim_core_build")
endif()

set(qsp_sim_core_FOUND TRUE)
