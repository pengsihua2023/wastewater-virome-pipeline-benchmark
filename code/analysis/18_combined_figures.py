#!/usr/bin/env python3
"""VIRAL figures for the combined cross-modality manuscript.

Everything is scored against the 151 non-phage viral ground truth; non-viral
detections and non-viral false positives are excluded throughout.

Seven figures written to figures_combined/ as 600-dpi PNG and vector PDF:

  Figure 1  study design
  Figure 2  viral recall / precision / F1 on both modalities
  Figure 3  viral F1 on both modalities, per pipeline
  Figure 4  viral TP / FN / FP on both modalities
  Figure 5  viral recall by designed abundance, both modalities
  Figure 6  viral recall-precision trade-off across thresholds, both modalities
  Figure 7  viral detection overlap between modalities

Changes made in response to the figure-only review (Peer-review-figures-V3.md):
  F1  Figure 3 was titled "Viral detection does not transfer between
      modalities" with rho as its evidence. The manuscript withdrew that
      inference, so the title is now descriptive and rho carries the same
      qualification as the legend.
  F2  Figure 1's annotation was struck through by an arrow, and Figure 6's
      shared x-axis label overwrote the tick labels. Both are rerouted.
  F3  Figure 5 used fixed 0-1 axes that left ~85% of the Illumina panel empty
      and collided the threshold annotations. Panels are now scaled to their
      data and only the primary operating point is annotated.
  m1  Figures used internal identifiers; they now use the manuscript's
      pipeline names.
  m2/m7  Zero counts are labelled explicitly instead of being drawn as stubs.
  m3  Series that are identically zero are offset slightly so they are all
      visible; stated in the manuscript legend.
  m4  The right-hand panel repeats its y tick labels.
  m5  The FP series is hatched so it does not rely on the green/red contrast.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"
FIG = AN / "figures_combined"
FIG.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight",
})

COL = {"TaxProfiler": "#0072B2", "TaxTriage": "#56B4E9", "mag-nf": "#009E73",
       "KrakenMetaReads-nf": "#E69F00", "GOTTCHA2": "#D55E00",
       "CLARK": "#CC79A7", "sourmash": "#7F7F7F"}
STYLE = {"TaxTriage": dict(ls="--", marker="s", mfc="none")}
# the names used in the manuscript text and tables, so figures and tables agree
DISPLAY = {"TaxProfiler": "nf-core/taxprofiler", "TaxTriage": "nf-core/taxtriage",
           "mag-nf": "nf-core/mag", "KrakenMetaReads-nf": "KrakenMetaReads-nf",
           "GOTTCHA2": "GOTTCHA2", "CLARK": "CLARK", "sourmash": "sourmash"}
# fixed order (by Illumina viral F1) so rank changes are visible across panels
ORDER = ["TaxProfiler", "TaxTriage", "mag-nf", "GOTTCHA2",
         "KrakenMetaReads-nf", "sourmash", "CLARK"]
LABELS = [DISPLAY[p] for p in ORDER]
BINS = ["<1e-4%", "1e-4-1e-3%", "1e-3-1e-2%", ">=1e-2%"]
BIN_N = {"<1e-4%": 15, "1e-4-1e-3%": 75, "1e-3-1e-2%": 43, ">=1e-2%": 18}
PLATS = ["Illumina", "PacBio"]
N_TRUTH = 151


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=600)
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print("  wrote", name)


def main():
    vm = pd.read_csv(INT / "viral_sample_metrics.tsv", sep="\t")
    vm = vm[vm.ground_truth == "non-phage virus (n=151)"]
    prim = vm[vm.min_reads == 1]
    agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    ab = pd.read_csv(INT / "viral_abundance_recall.tsv", sep="\t")
    ov = pd.read_csv(INT / "viral_detection_overlap.tsv", sep="\t")
    xs = np.arange(len(ORDER))

    # ------------------------------------------------------------- Figure 1
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.axis("off")
    boxes = [
        (0.02, 0.66, 0.20, 0.28, "Wastewater sample\nHTP-2024-01-15\nGrimm et al. 2026\nKraken2 clade reads"),
        (0.26, 0.66, 0.20, 0.28, "151-member\nviral ground truth\n(+ non-viral\nbackground)"),
        (0.50, 0.66, 0.20, 0.28, "One set of\nabundance files\ndistribution_0/1/2"),
        (0.75, 0.81, 0.23, 0.15, "CAMISIM 2.0\n→ 3 Illumina samples\n~16.5 M read pairs each"),
        (0.75, 0.62, 0.23, 0.15, "PBSIM3 + ccs\n→ 3 PacBio HiFi samples\n~479 k reads each"),
        (0.06, 0.26, 0.40, 0.20, "7 pipelines × 6 datasets\n= 42 scored runs"),
        (0.54, 0.26, 0.40, 0.20, "Viral-only scoring by NCBI taxid\nat ≥1 / ≥10 / ≥100 reads"),
        (0.20, 0.02, 0.60, 0.16, "Viral recall, precision and F1 for each pipeline,\nreported separately for each modality"),
    ]
    for x, y, w, h, t in boxes:
        fc = "#FDE9D9" if "reported separately" in t else "#F2F2F2"
        ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc,
                                   edgecolor="#333333", lw=0.9))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=6.6)
    for x0, y0, x1, y1 in [(0.22, 0.80, 0.26, 0.80), (0.46, 0.80, 0.50, 0.80),
                           (0.70, 0.80, 0.75, 0.885), (0.70, 0.80, 0.75, 0.695),
                           (0.46, 0.36, 0.54, 0.36), (0.50, 0.26, 0.50, 0.18)]:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    # F2: route the simulator -> pipelines arrow as an elbow BELOW the caption,
    # instead of a diagonal that struck the caption through
    ax.plot([0.865, 0.865, 0.30], [0.62, 0.555, 0.555], lw=0.9, color="#333333",
            solid_capstyle="butt")
    ax.annotate("", xy=(0.30, 0.463), xytext=(0.30, 0.545),
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    ax.text(0.30, 0.605, "the same abundance files drive both simulators,\n"
            "calibrated to the same 12% viral + phage read share (Section 2.2)",
            ha="center", va="center", fontsize=5.6, style="italic",
            color="#8B2500")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    save(fig, "Figure1_study_design")

    # ------------------------------------------------------------- Figure 2
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
    for axi, plat in zip(axes, PLATS):
        sub = prim[prim.platform == plat]
        for k, (m, c, h) in enumerate([("recall", "#BBBBBB", ""),
                                       ("precision", "#777777", "///"),
                                       ("f1", "#333333", "...")]):
            means = [np.nan_to_num(sub[sub.pipeline == p][m].mean()) for p in ORDER]
            axi.bar(xs + (k - 1) * 0.26, means, 0.26, color=c,
                    edgecolor="black", lw=0.4, hatch=h,
                    label=["Recall", "Precision", "F1"][k])
            for xi, p, mv in zip(xs + (k - 1) * 0.26, ORDER, means):
                v = sub[sub.pipeline == p][m].values
                v = v[np.isfinite(v)]
                axi.plot([xi] * len(v), v, "o", ms=2, color="#D55E00", zorder=3)
                if mv == 0:                      # m7: mark true zeros
                    axi.text(xi, 0.015, "0", ha="center", va="bottom",
                             fontsize=5, color="#333333")
        axi.set_xticks(xs)
        axi.set_xticklabels(LABELS, rotation=35, ha="right")
        axi.set_title(plat)
        axi.set_ylim(0, 1.05)
        axi.tick_params(labelleft=True)          # m4
    axes[0].set_ylabel(f"Score (viral species, n = {N_TRUTH})")
    axes[0].legend(frameon=False, ncol=3, loc="upper center")
    save(fig, "Figure2_viral_performance_both")

    # ------------------------------------------------------------- Figure 3
    p = agg.pivot_table(index="pipeline", columns="platform", values="f1_mean")
    rho = stats.spearmanr(p["Illumina"], p["PacBio"])
    fig, ax = plt.subplots(figsize=(5.0, 4.0))

    def spread(vals, gap):
        idx = sorted(range(len(vals)), key=lambda i: vals[i])
        out = list(vals)
        for k in range(1, len(idx)):
            a, b = idx[k - 1], idx[k]
            if out[b] - out[a] < gap:
                out[b] = out[a] + gap
        return out

    left = [p.loc[q, "Illumina"] for q in ORDER]
    right = [p.loc[q, "PacBio"] for q in ORDER]
    gap = (max(right + left) - min(right + left)) * 0.055
    l_lab, r_lab = spread(left, gap), spread(right, gap)
    for pipe, y0, y1, l0, l1 in zip(ORDER, left, right, l_lab, r_lab):
        ls = STYLE.get(pipe, {}).get("ls", "-")
        ax.plot([0, 1], [y0, y1], ls=ls, marker="o", ms=5, lw=1.6, color=COL[pipe])
        ax.annotate(DISPLAY[pipe], xy=(0, y0), xytext=(-0.06, l0),
                    textcoords="data", ha="right", va="center", fontsize=6.5,
                    color=COL[pipe],
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=COL[pipe],
                                    shrinkA=0, shrinkB=2))
        ax.annotate(f"{y1:.3f}", xy=(1, y1), xytext=(1.06, l1), textcoords="data",
                    ha="left", va="center", fontsize=6.5, color=COL[pipe],
                    arrowprops=dict(arrowstyle="-", lw=0.5, color=COL[pipe],
                                    shrinkA=0, shrinkB=2))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Illumina", "PacBio HiFi"])
    ax.set_xlim(-0.95, 1.35)
    ax.set_ylabel("Viral F1 (species level, ≥1 clade read)")
    # F1: descriptive title. The rank correlation is uninformative at n = 7 and
    # changes sign when the byte-identical Illumina duplicate is removed, so it
    # is reported as context rather than as the finding.
    ax.set_title("Viral F1 for each pipeline on both modalities\n"
                 f"Spearman ρ = {rho.statistic:+.2f} "
                 f"(p = {rho.pvalue:.2f}, n = 7; uninformative — see Table 5)",
                 fontsize=8.5)
    save(fig, "Figure3_viral_rank_transfer")

    # ------------------------------------------------------------- Figure 4
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
    for axi, plat in zip(axes, PLATS):
        c = prim[prim.platform == plat].groupby("pipeline")[["TP", "FN", "FP"]].mean()
        for k, (m, col, h) in enumerate([("TP", "#009E73", ""),
                                         ("FN", "#E69F00", ""),
                                         ("FP", "#D55E00", "///")]):   # m5
            vals = [c.loc[p_, m] for p_ in ORDER]
            pos = xs + (k - 1) * 0.27
            # m2: draw a bar only where the count is non-zero; label true zeros
            axi.bar([x for x, v in zip(pos, vals) if v > 0],
                    [v for v in vals if v > 0], 0.27, color=col,
                    edgecolor="black", lw=0.4, hatch=h, label=m)
            for x, v in zip(pos, vals):
                if v == 0:
                    axi.text(x, 0.45, "0", ha="center", va="bottom", fontsize=5,
                             color=col, fontweight="bold")
        axi.set_yscale("log")
        axi.set_ylim(0.4, 3e4)
        axi.axhline(N_TRUTH, ls="--", lw=0.7, color="#333333")
        axi.set_xticks(xs)
        axi.set_xticklabels(LABELS, rotation=35, ha="right")
        axi.set_title(plat)
        axi.tick_params(labelleft=True)          # m4
    axes[0].set_ylabel("Viral taxa (mean of 3 samples, log scale)")
    axes[0].legend(frameon=False, ncol=3, loc="upper left")
    save(fig, "Figure4_viral_TP_FN_FP_both")

    # ------------------------------------------------------------- Figure 5
    # F3: scale each panel to its own data instead of forcing 0-1, and annotate
    # only the primary operating point so the threshold labels stop colliding.
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3))
    for axi, plat in zip(axes, PLATS):
        sweep = (vm[vm.platform == plat]
                 .groupby(["pipeline", "min_reads"])[["recall", "precision"]]
                 .mean().reset_index())
        rmax = pmax = 0.0
        for oi, pipe in enumerate(ORDER):
            g = sweep[sweep.pipeline == pipe].sort_values("min_reads")
            g = g[np.isfinite(g.precision) & np.isfinite(g.recall)]
            if g.empty:
                continue
            rmax = max(rmax, g.recall.max()); pmax = max(pmax, g.precision.max())
            st = STYLE.get(pipe, dict(ls="-", marker="o"))
            axi.plot(g.recall, g.precision, ms=4, lw=1.2, color=COL[pipe],
                     label=DISPLAY[pipe], **st)
            one = g[g.min_reads == 1]
            if len(one):
                r0 = one.iloc[0]
                # stagger the offset so labels on near-coincident points
                # (several pipelines sit at recall 0 on PacBio) do not stack
                axi.annotate("≥1", (r0.recall, r0.precision), fontsize=5.5,
                             xytext=(4, 3 + 5 * (oi % 2)),
                             textcoords="offset points", color=COL[pipe])
        axi.set_xlabel("Viral recall")
        axi.set_title(plat)
        axi.set_xlim(-0.03 * max(rmax, 0.1), rmax * 1.14 + 0.03)
        axi.set_ylim(-0.03 * max(pmax, 0.1), pmax * 1.18 + 0.03)
        axi.tick_params(labelleft=True)
    axes[0].set_ylabel("Viral precision")
    axes[1].legend(frameon=False, fontsize=6, loc="upper left",
                   bbox_to_anchor=(1.02, 1))
    # NOTE (V5): numbered Figure 6, not 5 -- the manuscript now places figures
    # inline at first reference, and the abundance figure below (Section 3.5)
    # precedes this threshold figure (Section 3.6) in body order, so the
    # numbering swaps relative to the original Figure-Legends-only ordering.
    save(fig, "Figure6_viral_threshold_both")

    # ------------------------------------------------------------- Figure 5
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
    for axi, plat in zip(axes, PLATS):
        s = ab[ab.platform == plat]
        zero_seen = 0
        for pipe in ORDER:
            y = s[s.pipeline == pipe].set_index("bin").reindex(BINS).recall_mean
            yv = y.values.astype(float)
            # m3: series that are identically zero would hide one another
            if np.nanmax(yv) == 0:
                yv = yv + zero_seen * 0.012
                zero_seen += 1
            st = STYLE.get(pipe, dict(ls="-", marker="o"))
            axi.plot(range(len(BINS)), yv, ms=3.5, lw=1.2,
                     color=COL[pipe], label=DISPLAY[pipe], **st)
        axi.set_xticks(range(len(BINS)))
        axi.set_xticklabels([f"{b}\n(n={BIN_N[b]})" for b in BINS], fontsize=6)
        axi.set_title(plat)
        axi.set_ylim(-0.04, 1.05)
        axi.tick_params(labelleft=True)          # m4
    axes[0].set_ylabel("Viral recall")
    axes[1].legend(frameon=False, fontsize=6, loc="upper left",
                   bbox_to_anchor=(1.02, 1))
    # F2: make room below the tick labels instead of drawing the shared x-axis
    # label on top of them
    fig.subplots_adjust(bottom=0.30)
    fig.supxlabel("Designed relative abundance of the virus", fontsize=8, y=0.02)
    save(fig, "Figure5_viral_abundance_recall_both")

    # ------------------------------------------------------------- Figure 7
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    cats = ["Both\nmodalities", "PacBio\nonly", "Illumina\nonly", "Neither"]
    order = ["both modalities", "PacBio only", "Illumina only", "neither"]
    vals = [int(ov.loc[ov.category == c, "n_viruses"].iloc[0]) for c in order]
    cols = ["#009E73", "#D55E00", "#0072B2", "#333333"]
    bars = ax.bar(cats, vals, color=cols, edgecolor="black", lw=0.5)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 2.5, f"{v}\n({v/N_TRUTH:.1%})",
                ha="center", fontsize=7)
    ax.set_ylabel(f"Viruses (of {N_TRUTH})")
    ax.set_ylim(0, max(vals) * 1.35)
    ax.set_title("Viruses detected by at least one pipeline")
    save(fig, "Figure7_viral_detection_overlap")

    # drop the previous whole-community figure set
    for stale in ["Figure2_performance_both_modalities", "Figure3_rank_transfer",
                  "Figure4_domain_recall_both", "Figure5_TP_FN_FP_both",
                  "Figure6_threshold_tradeoff_both", "Figure7_abundance_recall_both",
                  "Figure8_detection_overlap"]:
        for ext in (".png", ".pdf"):
            f = FIG / f"{stale}{ext}"
            if f.exists():
                f.unlink()

    print(f"\nviral rank-transfer Spearman rho = {rho.statistic:+.3f} (p = {rho.pvalue:.3f})")
    print("overlap:", dict(zip(order, vals)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
