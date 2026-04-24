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

if(NOT TARGET qsp_sim_core::qsp_sim_core)
    add_subdirectory("${_qsp_sim_core_root}"
                     "${CMAKE_BINARY_DIR}/_qsp_sim_core_build")
endif()

set(qsp_sim_core_FOUND TRUE)
