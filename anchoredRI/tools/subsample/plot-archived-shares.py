"""Render the Task E archived component-share figure from the published fit table."""

import csv
from math import isfinite
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ARCHIVE = Path(__file__).resolve().parent / "archived-shares"
METHODS = (("anchoredRI", "Baseline", 3), ("armB", "Arm B", 5))


def main():
    with (ARCHIVE / "fits.csv").open(newline="") as source:
        records = list(csv.DictReader(source))
    panels = [
        (method, name, k)
        for method, name, cap in METHODS
        for k in range(1, cap + 1)
        if any(r["method"] == method and int(r["components"]) == k for r in records)
    ]
    if len(panels) != 6:
        raise ValueError("Expected the six selected method-by-component groups in Task E.")

    fig, axes = plt.subplots(2, 3, figsize=(11, 6.25), sharey=True)
    fig.subplots_adjust(left=0.077, right=0.985, bottom=0.152, top=0.875, hspace=0.49, wspace=0.10)
    ordinary, extreme = "#355e74", "#ae441a"
    for ax, (method, name, k) in zip(axes.flat, panels):
        rows = [r for r in records if r["method"] == method and int(r["components"]) == k]
        for r in rows:
            if not r["d"] or not isfinite(float(r["d"])) or float(r["d"]) <= 0:
                raise ValueError("The published log-scale plot requires a positive finite error per fit.")
        flagged = [r for r in rows if max(float(r[f] or 0) for f in ("abs_z_lower", "abs_z_upper")) > 5]
        flagged_ids = {r["Index"] for r in flagged}
        unflagged = [r for r in rows if r["Index"] not in flagged_ids]
        ax.scatter([float(r["gap_pp"]) for r in unflagged], [float(r["d"]) for r in unflagged],
                   color=ordinary, linewidths=0, s=11, alpha=0.57)
        ax.scatter([float(r["gap_pp"]) for r in flagged], [float(r["d"]) for r in flagged],
                   marker="^", color=extreme, linewidths=0, s=24, alpha=0.88)
        ax.set_title(f"{name} · {k} components", fontsize=11.5, fontweight="bold", loc="left", pad=6)
        ax.text(.98, 1.03, f"n = {len(rows):,}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=9, color="#464b4d")
        ax.set_yscale("log")
        ax.set_ylim(0.00045, 100)
        ax.set_yticks([0.001, 0.01, 0.1, 1, 10, 100], ["0.001", "0.01", "0.1", "1", "10", "100"])
        if k == 2:
            ax.set_xlim(-0.16, 0.025)
            ax.set_xticks([-.15, -.10, -.05, 0], ["−0.15", "−0.10", "−0.05", "0"])
        else:
            ax.set_xlim(-100, 55)
            ax.set_xticks([-80, -40, 0, 40])
        ax.axvline(0, color="#677278", linestyle=":", lw=0.8, zorder=0)
        ax.grid(axis="y", color="#d8dcde", lw=0.6)
        ax.set_axisbelow(True)
        ax.minorticks_off()
        ax.tick_params(labelsize=9, length=3, color="#9da4a9")
        for edge in ("left", "bottom"):
            ax.spines[edge].set_color("#b8bec1")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Component 1 share and reference-limit error", fontsize=15, fontweight="bold",
                 x=.077, ha="left", y=.991)
    fig.supylabel("Mean absolute endpoint z error", x=.012, fontsize=11)
    fig.supxlabel("Fitted share minus designed healthy share (percentage points)", y=.070, fontsize=11)
    fig.legend(handles=[
        Line2D([], [], marker="o", linestyle="none", color=ordinary, markersize=5,
               label="Every requested endpoint |z| ≤ 5"),
        Line2D([], [], marker="^", linestyle="none", color=extreme, markersize=6,
               label="Any requested endpoint |z| > 5"),
    ], loc="lower center", bbox_to_anchor=(.52,.009), frameon=False, ncol=2, fontsize=9)
    fig.savefig(ARCHIVE / "share-gap-vs-error.png", dpi=180, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
