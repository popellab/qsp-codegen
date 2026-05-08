#!/usr/bin/env python3
"""Compare MATLAB and C++ ODE trajectories.

Usage:
    python compare_trajectories.py <matlab_csv> <cpp_csv> [--rtol 1e-3] [--atol 1e-9]
"""

import argparse
import sys

import numpy as np


def load_csv(path):
    with open(path) as f:
        header = f.readline().strip().split(",")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    times = data[:, 0]
    values = data[:, 1:]
    names = header[1:]
    return times, values, names


def compare(matlab_csv, cpp_csv, rtol=1e-3, atol=1e-9):
    """Compare two trajectory files. Returns (pass, report_string)."""
    t_m, v_m, names_m = load_csv(matlab_csv)
    t_c, v_c, names_c = load_csv(cpp_csv)

    lines = []
    lines.append(f"MATLAB: {len(t_m)} time points, {len(names_m)} species")
    lines.append(f"C++:    {len(t_c)} time points, {len(names_c)} species")

    # MATLAB: "V_T_CD8", C++: "V_T.CD8" — normalize to underscores.
    def normalize(name):
        return name.replace(".", "_")

    common = []
    m_norm = {normalize(n): (n, i) for i, n in enumerate(names_m)}
    c_idx = {n: i for i, n in enumerate(names_c)}

    for name in names_c:
        key = normalize(name)
        if key in m_norm:
            m_name, mi = m_norm[key]
            common.append((name, mi, c_idx[name]))

    lines.append(f"Matched species: {len(common)} / {len(names_c)}")

    if len(common) == 0:
        lines.append("ERROR: No species matched between files!")
        return False, "\n".join(lines)

    # The C++ side runs CV_ONE_STEP with a min-cadence floor (qsp-codegen
    # v3 schema), so its sample times are non-uniform and don't line up
    # with MATLAB's fixed-grid output. Compare on MATLAB's grid using
    # linear interpolation of each C++ species column. Restrict to MATLAB
    # times within the C++ trajectory's span to avoid extrapolation.
    t_lo = max(t_m.min(), t_c.min())
    t_hi = min(t_m.max(), t_c.max())
    mask = (t_m >= t_lo) & (t_m <= t_hi)
    t_common = t_m[mask]
    if len(t_common) == 0:
        lines.append("ERROR: MATLAB and C++ time spans do not overlap.")
        return False, "\n".join(lines)
    lines.append(
        f"Comparing on MATLAB grid (linear-interp C++ → MATLAB times, "
        f"{len(t_common)} points in [{t_lo:.3f}, {t_hi:.3f}])"
    )

    n_fail = 0
    max_rdiff = 0.0
    worst_species = ""
    worst_time = 0.0

    for name, mi, ci in common:
        vm_col = v_m[mask, mi]
        vc_col = np.interp(t_common, t_c, v_c[:, ci])

        denom = np.maximum.reduce(
            [np.abs(vm_col), np.abs(vc_col), np.full_like(vm_col, atol)]
        )
        rdiff = np.abs(vm_col - vc_col) / denom

        # Track worst across all timepoints for this species.
        worst_idx = int(np.argmax(rdiff))
        if rdiff[worst_idx] > max_rdiff:
            max_rdiff = float(rdiff[worst_idx])
            worst_species = name
            worst_time = float(t_common[worst_idx])

        # Flag any timepoint exceeding both rtol and atol budgets.
        fail_mask = (rdiff > rtol) & (np.abs(vm_col - vc_col) > atol)
        n_fail += int(fail_mask.sum())
        for fi in np.where(fail_mask)[0][: max(0, 10 - n_fail + int(fail_mask.sum()))]:
            t = float(t_common[fi])
            lines.append(
                f"  FAIL: {name} at t={t:.2f}: "
                f"MATLAB={vm_col[fi]:.6e}, C++={vc_col[fi]:.6e}, "
                f"rdiff={rdiff[fi]:.2e}"
            )

    total_comparisons = len(common) * len(t_common)
    lines.append(f"Comparisons: {total_comparisons}")
    lines.append(f"Failures: {n_fail} (rtol={rtol}, atol={atol})")
    lines.append(
        f"Worst relative diff: {max_rdiff:.2e} ({worst_species} at t={worst_time:.2f})"
    )

    passed = n_fail == 0
    lines.append(f"Result: {'PASS' if passed else 'FAIL'}")

    return passed, "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Compare MATLAB vs C++ ODE trajectories"
    )
    parser.add_argument("matlab_csv")
    parser.add_argument("cpp_csv")
    parser.add_argument("--rtol", type=float, default=1e-3)
    parser.add_argument("--atol", type=float, default=1e-9)
    args = parser.parse_args()

    passed, report = compare(args.matlab_csv, args.cpp_csv, args.rtol, args.atol)
    print(report)
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
