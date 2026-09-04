"""Single source of truth for every figure in the MIMIC-IV validation paper.
Mirrors make_charts.py's role for the FAERS paper: make_mimic4_paper.py
embeds the PNGs written here rather than generating its own matplotlib code.

Palette: the same validated colorblind-safe categorical palette already in
use (and already checked) for the sibling FAERS paper's charts
(make_charts.py) -- reused as-is for visual consistency between the two
companion papers, not re-validated since the hex values are unchanged.

Usage:
    python make_mimic4_charts.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

CH = "mimic4_charts"
os.makedirs(CH, exist_ok=True)
OUT = "mimic4_pipeline/output"

# ---- validated colorblind-safe palette (same constants as make_charts.py) --
BLUE, AQUA, YELLOW, GREEN, VIOLET = "#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7"
RED = "#e34948"
GRAY_MID = "#f0efec"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "text.color": INK,
        "axes.edgecolor": GRID,
        "axes.labelcolor": SECONDARY_INK,
        "xtick.color": SECONDARY_INK,
        "ytick.color": SECONDARY_INK,
        "axes.grid": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)

facts = json.load(open("/private/tmp/claude-501/-Users-praffulg-projects-Drugs-Research-Paper/"
                        "9224774d-5b97-413b-ac4a-557519d945e6/scratchpad/paper_facts.json"))

threshold = pd.read_csv(f"{OUT}/analysis_clinical_threshold_crossing.csv")
adjusted = pd.read_csv(f"{OUT}/analysis_severity_adjusted.csv")
matched_threshold = pd.read_csv(f"{OUT}/analysis_matched_threshold_crossing.csv")
balance = pd.read_csv(f"{OUT}/propensity_balance.csv")

ANALYTES = ["ALT", "AST", "Total Bilirubin", "Sodium"]


def odds_ratio_from_counts(a_events, a_n, b_events, b_n):
    """Simple (unconditional) 2x2 odds ratio with a Wald CI on the log
    scale; Haldane-Anscombe continuity correction (+0.5 to every cell) if
    any cell is zero, so this never divides by zero."""
    a1, a0, b1, b0 = a_events, a_n - a_events, b_events, b_n - b_events
    if min(a1, a0, b1, b0) == 0:
        a1, a0, b1, b0 = a1 + 0.5, a0 + 0.5, b1 + 0.5, b0 + 0.5
    or_ = (a1 / a0) / (b1 / b0)
    se = np.sqrt(1 / a1 + 1 / a0 + 1 / b1 + 1 / b0)
    log_or = np.log(or_)
    return or_, np.exp(log_or - 1.96 * se), np.exp(log_or + 1.96 * se)


def thresh_counts(df, analyte, group):
    row = df[(df["analyte"] == analyte) & (df["group"] == group)].iloc[0]
    return int(row["n_crossing"]), int(row["n_evaluable"])


# ============================================================ FIGURE 1
# Forest plot: OR (Tricyclic vs SSRI+SNRI) for crossing the pre-specified
# clinical threshold, at three progressively more rigorous analytic layers.
def figure_1_forest_plot():
    layers = ["Unadjusted", "Charlson-adjusted", "Propensity-matched"]
    colors = {"Unadjusted": MUTED, "Charlson-adjusted": AQUA, "Propensity-matched": BLUE}

    rows = []
    for analyte in ANALYTES:
        a_ev, a_n = thresh_counts(threshold, analyte, "Tricyclic")
        b_ev, b_n = thresh_counts(threshold, analyte, "SSRI+SNRI")
        or_u, lo_u, hi_u = odds_ratio_from_counts(a_ev, a_n, b_ev, b_n)
        rows.append({"analyte": analyte, "layer": "Unadjusted", "or": or_u, "lo": lo_u, "hi": hi_u})

        adj_row = adjusted[adjusted["analyte"] == analyte].iloc[0]
        rows.append({"analyte": analyte, "layer": "Charlson-adjusted", "or": adj_row["adjusted_or"],
                     "lo": adj_row["adjusted_or_ci_low"], "hi": adj_row["adjusted_or_ci_high"]})

        a_ev, a_n = thresh_counts(matched_threshold, analyte, "Tricyclic")
        b_ev, b_n = thresh_counts(matched_threshold, analyte, "SSRI+SNRI")
        or_m, lo_m, hi_m = odds_ratio_from_counts(a_ev, a_n, b_ev, b_n)
        rows.append({"analyte": analyte, "layer": "Propensity-matched", "or": or_m, "lo": lo_m, "hi": hi_m})

    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    y_positions, y_labels = [], []
    y = 0
    group_gap = 1.0
    within_gap = 0.9
    for analyte in reversed(ANALYTES):
        for layer in reversed(layers):
            r = df[(df["analyte"] == analyte) & (df["layer"] == layer)].iloc[0]
            ax.plot([r["lo"], r["hi"]], [y, y], color=colors[layer], linewidth=2, solid_capstyle="round", zorder=2)
            ax.plot(r["or"], y, "o", color=colors[layer], markersize=7, zorder=3,
                     markeredgecolor="white", markeredgewidth=0.8)
            y_positions.append(y)
            y_labels.append(f"{analyte} -- {layer}")
            y += within_gap
        y += group_gap - within_gap

    ax.axvline(1.0, color=SECONDARY_INK, linewidth=1, linestyle="--", zorder=1)
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio (Tricyclic vs. SSRI+SNRI), log scale")
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.set_xlim(0.5, 2.2)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Effect size across three analytic layers, by outcome", fontsize=12, color=INK, pad=12)

    handles = [plt.Line2D([0], [0], color=colors[l], marker="o", linewidth=2, markersize=7, label=l) for l in layers]
    ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{CH}/forest_plot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================ FIGURE 2
# Love plot: covariate balance (standardized mean difference) before vs.
# after propensity-score matching.
def figure_2_love_plot():
    b = balance.copy()
    b = b.reindex(b["smd_before"].abs().sort_values(ascending=True).index)
    labels = b["covariate"].str.replace("_", " ").str.capitalize().tolist()
    y = np.arange(len(b))

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.axvline(0, color=SECONDARY_INK, linewidth=1, zorder=1)
    ax.axvline(0.1, color=GRID, linewidth=1, linestyle="--", zorder=1)
    ax.axvline(-0.1, color=GRID, linewidth=1, linestyle="--", zorder=1)

    for yi, (_, r) in zip(y, b.iterrows()):
        ax.plot([r["smd_before"], r["smd_after"]], [yi, yi], color=GRID, linewidth=1, zorder=1)
    ax.plot(b["smd_before"], y, "o", color=MUTED, markersize=6, label="Before matching", zorder=2,
             markeredgecolor="white", markeredgewidth=0.6)
    ax.plot(b["smd_after"], y, "o", color=BLUE, markersize=6, label="After matching", zorder=3,
             markeredgecolor="white", markeredgewidth=0.6)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Standardized mean difference (Tricyclic - SSRI+SNRI)")
    ax.set_title("Covariate balance before and after propensity-score matching", fontsize=12, color=INK, pad=12)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{CH}/love_plot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================ FIGURE 3
# CONSORT-style cohort flow diagram.
def figure_3_consort_diagram():
    fig, ax = plt.subplots(figsize=(8.5, 10))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 130)
    ax.axis("off")

    def box(cx, cy, w, h, text, fc=GRAY_MID, ec=SECONDARY_INK, fontsize=8.5, textcolor=INK):
        b = FancyBboxPatch((cx - w / 2, cy - h / 2), w, h, boxstyle="round,pad=0.6,rounding_size=2",
                            linewidth=1.1, edgecolor=ec, facecolor=fc, zorder=2)
        ax.add_patch(b)
        ax.text(cx, cy, text, ha="center", va="center", fontsize=fontsize, color=textcolor, zorder=3, linespacing=1.4)
        return (cx, cy - h / 2), (cx, cy + h / 2), (cx - w / 2, cy), (cx + w / 2, cy)

    def arrow(p1, p2, color=SECONDARY_INK):
        a = FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=12, linewidth=1.1,
                             color=color, zorder=1, shrinkA=0, shrinkB=2)
        ax.add_patch(a)

    f = facts
    main_w = 62

    bottom_a, top_a, _, _ = box(50, 122, main_w, 10, f"All hospital admissions\nn = {f['n_total_admissions']:,}")
    bottom_b, top_b, _, _ = box(50, 106, main_w, 10,
                                 f"Systemic exposure to a curated antidepressant\nn = {f['n_exposed']:,}")
    arrow(bottom_a, top_b)

    _, _, left_c, right_c = box(50, 90, main_w, 12,
                                 f">1 antidepressant class present (raw)\nn = {f['n_raw_multiclass']:,}")
    arrow(bottom_b, (50, 96))

    box(84, 90, 30, 10, f"Excluded:\ngenuine multi-class\nn = {f['n_excluded_multiclass']:,}", fc="white", fontsize=8)
    arrow(right_c, (69, 90))
    box(84, 76, 30, 10, f"Excluded: every class\nlow-dose-ambiguous only\nn = {f['n_excluded_ambiguous_low_dose_only']:,}",
        fc="white", fontsize=8)
    arrow((84, 85), (84, 81), color=SECONDARY_INK)

    bottom_c = (50, 84)
    bottom_d, top_d, _, _ = box(50, 62, main_w, 12,
                                 f"Single-class exposed\n(primary analytic population)\nn = {f['n_single_class']:,}")
    ax.text(15, 74, f"resolved by dose-based\nreclassification, n = {f['n_dose_reclassified']:,}",
            ha="right", va="center", fontsize=7.3, color=SECONDARY_INK, style="italic")
    arrow(bottom_c, top_d)

    classes = [("Tricyclic", f["by_class"]["Tricyclic"]), ("SSRI", f["by_class"]["SSRI"]),
               ("SNRI", f["by_class"]["SNRI"]), ("Atypical", f["by_class"]["Atypical"]),
               ("RIMA", f["by_class"]["RIMA"])]
    xs = np.linspace(12, 88, len(classes))
    leaf_colors = {"Tricyclic": RED, "SSRI": BLUE, "SNRI": VIOLET, "Atypical": AQUA, "RIMA": YELLOW}
    for (cls, n), x in zip(classes, xs):
        _, leaf_top, _, _ = box(x, 34, 15.5, 12, f"{cls}\nn = {n:,}", fc="white", ec=leaf_colors[cls], fontsize=8)
        arrow((50, 56), (x, 40), color=MUTED)

    fig.suptitle("Cohort construction flow", fontsize=12, color=INK, y=0.97)
    fig.tight_layout()
    fig.savefig(f"{CH}/consort_diagram.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    figure_1_forest_plot()
    figure_2_love_plot()
    figure_3_consort_diagram()
    print(f"Wrote 3 figures -> {CH}/")
