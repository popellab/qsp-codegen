"""Locate the C++ tree shipped with qsp-codegen.

Consumers wire this into their CMake by reading these paths from a
small Python shim in their own CMakeLists invocation, typically via::

    execute_process(COMMAND python -m qsp_codegen.cmake --prefix
                    OUTPUT_VARIABLE QSP_SIM_CORE_CMAKE_DIR
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    list(APPEND CMAKE_PREFIX_PATH "${QSP_SIM_CORE_CMAKE_DIR}")
    find_package(qsp_sim_core CONFIG REQUIRED)

We keep this helper tiny and stdlib-only so it can run before any
project-level Python env is activated.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def cpp_root() -> Path:
    """Return the directory holding CMakeLists.txt for qsp_sim_core."""
    return Path(__file__).resolve().parent / "cpp"


def cmake_dir() -> Path:
    """Return the directory containing qsp_sim_coreConfig.cmake.

    `find_package` looks for `<name>Config.cmake` under a path on
    CMAKE_PREFIX_PATH, so consumers add this directory to the prefix
    path before the find_package call.
    """
    return cpp_root() / "cmake"


def driver_source() -> Path:
    return cpp_root() / "src" / "qsp_sim_main.cpp"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Print paths into qsp-codegen's shipped C++ tree."
    )
    parser.add_argument(
        "--prefix", action="store_true",
        help="Print the CMake package prefix (dir with qsp_sim_coreConfig.cmake)",
    )
    parser.add_argument(
        "--cpp-root", action="store_true",
        help="Print the C++ source root (dir with CMakeLists.txt)",
    )
    parser.add_argument(
        "--driver-source", action="store_true",
        help="Print the absolute path to qsp_sim_main.cpp",
    )
    args = parser.parse_args()
    if args.prefix:
        print(cmake_dir())
    elif args.cpp_root:
        print(cpp_root())
    elif args.driver_source:
        print(driver_source())
    else:
        parser.error("pick exactly one of --prefix / --cpp-root / --driver-source")


if __name__ == "__main__":
    main()
