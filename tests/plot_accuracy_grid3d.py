#!/usr/bin/env python3
"""
Visualize the accuracy of the 3D rectilinear grid methods (non OpenCL).

Reads accuracy_grid3d.csv (produced by accuracy_grid3d) and, for each model,
draws grouped bar charts of the mean relative error and the raytracing time.
Bars are grouped by (resolution, precision) so the convergence from medium to
fine resolution -- and the float/double comparison -- can be read directly.

The script adapts to whatever combinations are present in the CSV, so partial
runs (e.g. medium only, or a single method) plot fine.

Usage:
    python plot_accuracy_grid3d.py [accuracy_grid3d.csv]
"""

import csv
import sys

import numpy as np
import matplotlib.pyplot as plt


def load(filename):
    rows = []
    with open(filename, newline="") as f:
        for row in csv.DictReader(f):
            row["error"] = float(row["error"])
            row["time"] = float(row["time"])
            # tolerate an old CSV without a resolution column
            row.setdefault("resolution", "medium")
            rows.append(row)
    return rows


def ordered_unique(values, preferred):
    """unique values, ordered by `preferred` first then appearance"""
    seen = list(dict.fromkeys(values))
    out = [p for p in preferred if p in seen]
    out += [v for v in seen if v not in out]
    return out


def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "accuracy_grid3d.csv"
    rows = load(filename)

    models = ordered_unique([r["model"] for r in rows],
                            ["layers", "gradient", "constant"])
    methods = ordered_unique([r["method"] for r in rows],
                             ["FAST_SWEEPING", "SHORTEST_PATH",
                              "DYNAMIC_SHORTEST_PATH"])
    resolutions = ordered_unique([r["resolution"] for r in rows],
                                 ["medium", "fine"])
    precisions = ordered_unique([r["precision"] for r in rows],
                                ["double", "float"])

    # the bar groups: one per (resolution, precision) combination present
    groups = [(res, prec) for res in resolutions for prec in precisions]

    # lookup[(model, method, resolution, precision)] -> row
    lookup = {(r["model"], r["method"], r["resolution"], r["precision"]): r
              for r in rows}

    # colour by resolution, lighten for float so both stay distinguishable
    base = {"medium": "#4C72B0", "fine": "#DD8452"}
    extra = ["#55A868", "#C44E52", "#8172B3", "#937860"]

    def color_for(res, prec):
        c = base.get(res)
        if c is None:
            c = extra[resolutions.index(res) % len(extra)]
        if prec != precisions[0]:           # lighten the second precision
            rgb = np.array(plt.matplotlib.colors.to_rgb(c))
            c = tuple(rgb + (1.0 - rgb) * 0.45)
        return c

    nmodels = len(models)
    fig, axes = plt.subplots(2, nmodels, figsize=(5.5 * nmodels, 9),
                             squeeze=False)

    x = np.arange(len(methods))
    width = 0.8 / max(len(groups), 1)

    for col, model in enumerate(models):
        ax_err = axes[0][col]
        ax_time = axes[1][col]

        for gi, (res, prec) in enumerate(groups):
            errs, times = [], []
            for m in methods:
                r = lookup.get((model, m, res, prec))
                errs.append(r["error"] if r else np.nan)
                times.append(r["time"] if r else np.nan)
            offset = (gi - (len(groups) - 1) / 2) * width
            label = f"{res} / {prec}"
            color = color_for(res, prec)
            ax_err.bar(x + offset, errs, width, label=label, color=color)
            ax_time.bar(x + offset, times, width, label=label, color=color)

        ax_err.set_title(f"{model} : mean relative error")
        ax_err.set_yscale("log")
        ax_err.set_ylabel("mean relative error")
        ax_err.set_xticks(x)
        ax_err.set_xticklabels(methods, rotation=20, ha="right")
        ax_err.grid(True, axis="y", which="both", ls=":", alpha=0.5)
        ax_err.legend(fontsize=8)

        ax_time.set_title(f"{model} : raytracing time")
        ax_time.set_yscale("log")
        ax_time.set_ylabel("time (s)")
        ax_time.set_xticks(x)
        ax_time.set_xticklabels(methods, rotation=20, ha="right")
        ax_time.grid(True, axis="y", which="both", ls=":", alpha=0.5)
        ax_time.legend(fontsize=8)

    fig.suptitle("Accuracy of 3D rectilinear grid methods (non OpenCL)",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = "accuracy_grid3d.png"
    fig.savefig(out, dpi=150)
    print(f"Figure saved to {out}")
    plt.show()


if __name__ == "__main__":
    main()
