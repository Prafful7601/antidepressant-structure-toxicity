"""
Single source of truth for every figure in the manuscript. Both make_paper.py
and make_results_doc.py embed the PNGs written here rather than each
generating their own (previously duplicated) matplotlib code.

Palette: validated colorblind-safe categorical/diverging hues (dataviz skill
reference palette, references/palette.md) instead of the original red/green
scheme, which is not colorblind-safe.

Usage:
    python make_charts.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import spearmanr

CH = "charts"
os.makedirs(CH, exist_ok=True)
ORGANS = ["cardiac", "hepatic", "neurologic", "metabolic"]

# ---- validated colorblind-safe palette (dataviz skill reference palette) ----
BLUE, AQUA, YELLOW, GREEN, VIOLET = "#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7"
RED = "#e34948"
GRAY_MID = "#f0efec"
INK = "#0b0b0b"
SECONDARY_INK = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"

CLASS_COLORS = {"SSRI": BLUE, "TCA": RED, "Atypical": AQUA,
               "SNRI": VIOLET, "RIMA": YELLOW}
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "blue_gray_red", [BLUE, GRAY_MID, RED])

plt.rcParams.update({
    "font.family": "sans-serif", "axes.edgecolor": MUTED,
    "axes.labelcolor": INK, "text.color": INK, "xtick.color": SECONDARY_INK,
    "ytick.color": SECONDARY_INK, "axes.grid": False,
})


def chart_trs_validation():
    """Two panels: in-sample (circular) fit vs honest leave-one-out prediction.
    Showing both, side by side, is the point - a reader should never see the
    in-sample number without the LOOCV number next to it."""
    trs = pd.read_csv("drug_trs.csv")
    rho_in, p_in = spearmanr(trs["TRS_fitted"], trs["observed_tox_index"])
    rho_loo, p_loo = spearmanr(trs["TRS_loocv"], trs["observed_tox_index"])

    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6))

    ax = axes[0]
    x, y = trs["TRS_fitted"], trs["observed_tox_index"]
    ax.scatter(x, y, s=42, color=BLUE, edgecolor="white", linewidth=0.6, zorder=3)
    m, b = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, m * xs + b, color=RED, lw=1.8, zorder=2)
    ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
    ax.set_xlabel("Fitted TRS (structure only)")
    ax.set_ylabel("Observed toxicity index (mean ln ROR)")
    ax.set_title(f"In-sample fit (circular)\nSpearman rho = {rho_in:+.3f}, "
                 f"p = {p_in:.4f} (n={len(trs)})", fontsize=10)

    ax = axes[1]
    x2, y2 = trs["TRS_loocv"], trs["observed_tox_index"]
    ax.scatter(x2, y2, s=42, color=VIOLET, edgecolor="white", linewidth=0.6, zorder=3)
    ax.axhline(0, color=GRID, lw=0.8); ax.axvline(0, color=GRID, lw=0.8)
    ax.set_xlabel("Leave-one-drug-out predicted TRS")
    ax.set_ylabel("Observed toxicity index (mean ln ROR)")
    ax.set_title(f"Honest out-of-sample prediction\nSpearman rho = {rho_loo:+.3f}, "
                 f"p = {p_loo:.4f} (n={len(trs)})", fontsize=10)

    fig.suptitle("Structure-based Toxicity Risk Score: circular fit vs. honest "
                 "cross-validation", fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(f"{CH}/trs_validation.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return rho_in, p_in, rho_loo, p_loo


def chart_ror_heatmap():
    """Diverging blue (under-reported) <-> gray (parity) <-> red (over-reported)
    ROR heatmap. Cells meeting the PRIMARY signal (ROR lower 95% bound > 2 and
    >=3 cases) are boxed; a footnote marks the (much larger) legacy ROR>1 set."""
    sig = pd.read_csv("faers_signals.csv")
    wide = (sig.pivot_table(index="drug_name", columns="organ_system",
                            values="ror", aggfunc="first").reindex(columns=ORGANS))
    wide = wide.loc[wide.mean(axis=1).sort_values(ascending=False).index]
    sig_wide = (sig.pivot_table(index="drug_name", columns="organ_system",
                                values="signal", aggfunc="first")
                .reindex(index=wide.index, columns=ORGANS))
    legacy_wide = (sig.pivot_table(index="drug_name", columns="organ_system",
                                   values="signal_ror1", aggfunc="first")
                  .reindex(index=wide.index, columns=ORGANS))

    fig, ax = plt.subplots(figsize=(6.2, 8.6))
    norm = TwoSlopeNorm(vmin=0.3, vcenter=1.0, vmax=3.0)
    im = ax.imshow(wide.values, cmap=DIVERGING_CMAP, norm=norm, aspect="auto")
    ax.set_xticks(range(len(ORGANS)))
    ax.set_xticklabels([o.capitalize() for o in ORGANS], rotation=30, ha="right")
    ax.set_yticks(range(len(wide))); ax.set_yticklabels(wide.index, fontsize=7.5)
    for i in range(len(wide)):
        for j in range(len(ORGANS)):
            v = wide.values[i, j]
            if np.isnan(v):
                continue
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.2,
                    color=INK if 0.55 < v < 1.9 else "white")
            if sig_wide.values[i, j] == 1:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor=INK, linewidth=1.8))
            elif legacy_wide.values[i, j] == 1:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                           edgecolor=MUTED, linewidth=0.8,
                                           linestyle=(0, (2, 1))))
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("ROR (blue = under-reported, red = over-reported)", fontsize=8)
    ax.set_title("Reporting Odds Ratio by drug and organ system\n"
                 "Solid box = primary signal (ROR lower 95% CI > 2, >=3 cases); "
                 "dashed = legacy ROR>1 only", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{CH}/ror_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def chart_trs_weights():
    """TRS feature weights. Blue/red diverging instead of the original
    red/green (red-green is not colorblind-safe)."""
    wts = pd.read_csv("trs_weights.csv").sort_values("weight")
    colors = [BLUE if w < 0 else RED for w in wts["weight"]]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.barh(wts["feature"], wts["weight"], color=colors)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("TRS weight  (red = toxicophore, blue = safeophore)")
    ax.set_title("Fitted structural-feature weights (Toxicity Risk Score)\n"
                 "None survive Benjamini-Hochberg FDR correction at q<0.05 - "
                 "exploratory", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{CH}/trs_weights.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def chart_class_toxicity():
    """Per-drug toxicity index, fixed categorical hue order (never cycled)."""
    rank = pd.read_csv("class_safety_ranking.csv")
    rk = rank.sort_values(["class", "toxicity_index"])
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    ypos = range(len(rk))
    ax.barh(list(ypos), rk["toxicity_index"],
            color=[CLASS_COLORS.get(c, MUTED) for c in rk["class"]])
    ax.set_yticks(list(ypos)); ax.set_yticklabels(rk["drug_name"], fontsize=7)
    ax.invert_yaxis(); ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("Toxicity index (mean ln ROR; < 0 = safer than class avg)")
    ax.set_title("Per-drug toxicity index, coloured by pharmacological class")
    handles = [plt.Rectangle((0, 0), 1, 1, color=CLASS_COLORS[c])
              for c in CLASS_COLORS]
    ax.legend(handles, CLASS_COLORS.keys(), fontsize=7.5, loc="lower right",
             frameon=False)
    fig.tight_layout()
    fig.savefig(f"{CH}/class_toxicity.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def chart_feature_stability():
    """Rank-biserial effect size per feature x organ (small-multiple heatmap).
    Directly shows features whose role FLIPS between organs (e.g. carboxylic
    acid: hepatic toxicophore, neurologic safeophore) - a feature cannot be
    labelled 'the' toxicophore/safeophore without specifying the organ."""
    corr = pd.read_csv("structure_toxicity_correlation_continuous.csv")
    wide = corr.pivot_table(index="feature", columns="organ_system",
                            values="rank_biserial").reindex(columns=ORGANS)
    wide = wide.loc[wide.abs().max(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(5.6, 6.2))
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    im = ax.imshow(wide.values, cmap=DIVERGING_CMAP, norm=norm, aspect="auto")
    ax.set_xticks(range(len(ORGANS)))
    ax.set_xticklabels([o.capitalize() for o in ORGANS], rotation=30, ha="right")
    ax.set_yticks(range(len(wide))); ax.set_yticklabels(wide.index, fontsize=7.5)
    for i in range(len(wide)):
        for j in range(len(ORGANS)):
            v = wide.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=6,
                        color=INK if abs(v) < 0.6 else "white")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("Rank-biserial effect (red = toxicophore, blue = safeophore)",
                fontsize=7.5)
    ax.set_title("Structural-feature effect direction by organ system\n"
                 "(sign flips across organs = feature is NOT globally safe/unsafe)",
                fontsize=8.5)
    fig.tight_layout()
    fig.savefig(f"{CH}/feature_stability.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def generate_all():
    rho_in, p_in, rho_loo, p_loo = chart_trs_validation()
    chart_ror_heatmap()
    chart_trs_weights()
    chart_class_toxicity()
    chart_feature_stability()
    print(f"Wrote 5 charts to {os.path.abspath(CH)}")
    return {"rho_in": rho_in, "p_in": p_in, "rho_loo": rho_loo, "p_loo": p_loo}


if __name__ == "__main__":
    generate_all()
