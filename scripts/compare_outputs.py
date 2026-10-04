#!/usr/bin/env python3
"""Compare two whitespace-separated numeric files element by element.

Usage: compare_outputs.py RESULT REFERENCE [--rtol R] [--atol A]

An element passes when |result - reference| <= atol + rtol * |reference|
(the numpy.isclose criterion). Prints max absolute / relative differences
and the relative 2-norm error, and exits non-zero on any failure or shape
mismatch. Standard library only, so it runs anywhere python3 does.
"""
import argparse
import math
import sys


def read_values(path):
    with open(path) as fh:
        rows = [line.split() for line in fh if line.strip()]
    return [len(r) for r in rows], [float(x) for r in rows for x in r]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("result")
    ap.add_argument("reference")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--atol", type=float, default=1e-9)
    args = ap.parse_args()

    shape_res, res = read_values(args.result)
    shape_ref, ref = read_values(args.reference)
    if shape_res != shape_ref:
        print(f"FAIL: shape mismatch ({len(shape_res)} rows vs {len(shape_ref)} rows, "
              f"{len(res)} vs {len(ref)} values)")
        return 1

    max_abs = max_rel = 0.0
    n_fail = 0
    diff_sq = ref_sq = 0.0
    for a, b in zip(res, ref):
        d = abs(a - b)
        max_abs = max(max_abs, d)
        if b != 0.0:
            max_rel = max(max_rel, d / abs(b))
        if not d <= args.atol + args.rtol * abs(b):  # also catches NaN
            n_fail += 1
        diff_sq += (a - b) ** 2
        ref_sq += b ** 2
    rel_norm = math.sqrt(diff_sq / ref_sq) if ref_sq > 0 else math.sqrt(diff_sq)

    status = "PASS" if n_fail == 0 else "FAIL"
    print(f"{status}: {len(ref)} values, max |diff| = {max_abs:.3e}, "
          f"max rel diff = {max_rel:.3e}, rel 2-norm error = {rel_norm:.3e} "
          f"(rtol={args.rtol:g}, atol={args.atol:g}, {n_fail} out of tolerance)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
