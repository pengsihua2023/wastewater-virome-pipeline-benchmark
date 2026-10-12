#!/usr/bin/env python3
"""Publication figures. Every figure is written as 600-dpi PNG and vector PDF."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"
FIG = AN / "_superseded" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight",
})

ORDER = ["TaxProfiler", "TaxTriage", "mag-nf", "KrakenMetaReads-nf",
         "GOTTCHA2", "CLARK", "sourmash"]
# Okabe-Ito colour-blind-safe palette
COL = {"TaxProfiler": "#0072B2", "TaxTriage": "#56B4E9", "mag-nf": "#009E73",
       "KrakenMetaReads-nf": "#E69F00", "GOTTCHA2": "#D55E00",
       "CLARK": "#CC79A7", "sourmash": "#7F7F7F"}
PLATCOL = {"Illumina": "#0072B2", "PacBio": "#D55E00"}


def save(fig, name):
    fig.savefig(FIG / f"{name}.png", dpi=600)
    fig.savefig(FIG / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name)


def main():
    prim = pd.read_csv(INT / "sample_level_metrics.tsv", sep="\t")
    met = prim[prim.min_reads == 1].copy()
    agg = pd.read_csv(INT / "aggregate_metrics.tsv", sep="\t")
    plat = pd.read_csv(INT / "platform_contrast.tsv", sep="\t")
    dom = pd.read_csv(INT / "domain_recall.tsv", sep="\t")
    ab = pd.read_csv(INT / "abundance_recall.tsv", sep="\t")
    cor = pd.read_csv(INT / "abundance_correlation.tsv", sep="\t")
    det = pd.read_csv(INT / "detection_detail.tsv", sep="\t", low_memory=False)
    comp = pd.read_csv(INT / "realised_composition.tsv", sep="\t")

    # ---------------------------------------------------------------- Figure 1
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    ax.axis("off")
    boxes = [
        (0.02, 0.55, 0.20, 0.32, "Reference community\nHTP-2024-01-15\n(Grimm et al. 2026)\nKraken2 clade reads"),
        (0.26, 0.55, 0.20, 0.32, "708-member design\n490 bacteria/archaea\n67 phage\n151 virus"),
        (0.50, 0.55, 0.20, 0.32, "CAMISIM 2.0\nsimulation"),
        (0.74, 0.75, 0.24, 0.20, "Illumina\n3 samples\n16.6M reads each"),
        (0.74, 0.52, 0.24, 0.20, "PacBio HiFi\n3 samples\n0.45-0.52M reads"),
        (0.02, 0.08, 0.44, 0.32, "7 classification pipelines\nTaxProfiler · TaxTriage · mag-nf\nKrakenMetaReads-nf · GOTTCHA2\nCLARK · sourmash"),
        (0.52, 0.08, 0.46, 0.32, "Species-level scoring by NCBI taxid\nTP / FP / FN → recall, precision, F1\nstratified by domain and abundance"),
    ]
    for x, y, w, h, t in boxes:
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=True, facecolor="#F2F2F2",
                                   edgecolor="#333333", lw=0.8))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=6.6)
    for x0, y0, x1, y1 in [(0.22, 0.71, 0.26, 0.71), (0.46, 0.71, 0.50, 0.71),
                           (0.70, 0.71, 0.74, 0.85), (0.70, 0.71, 0.74, 0.62),
                           (0.50, 0.52, 0.30, 0.40), (0.50, 0.24, 0.52, 0.24)]:
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    save(fig, "Figure1_study_design")

    # ---------------------------------------------------------------- Figure 2
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for axi, platform in zip(axes, ["Illumina", "PacBio"]):
        sub = met[met.platform == platform]
        xs = np.arange(len(ORDER))
        w = 0.26
        for k, (m, hatch) in enumerate([("recall", ""), ("precision", "///"), ("f1", "...")]):
            means = [sub[sub.pipeline == p][m].mean() for p in ORDER]
            pts = [sub[sub.pipeline == p][m].values for p in ORDER]
            axi.bar(xs + (k - 1) * w, means, w, color=["#BBBBBB", "#777777", "#333333"][k],
                    edgecolor="black", lw=0.4, hatch=hatch,
                    label=["Recall", "Precision", "F1"][k])
            for xi, v in zip(xs + (k - 1) * w, pts):
                axi.plot([xi] * len(v), v, "o", ms=1.8, color="#D55E00", zorder=3)
        axi.set_xticks(xs)
        axi.set_xticklabels(ORDER, rotation=35, ha="right")
        axi.set_title(platform)
        axi.set_ylim(0, 1.05)
    axes[0].set_ylabel("Score (species level)")
    axes[0].legend(frameon=False, loc="upper right")
    save(fig, "Figure2_recall_precision_f1")

    # ---------------------------------------------------------------- Figure 3
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for axi, platform in zip(axes, ["Illumina", "PacBio"]):
        sub = met[met.platform == platform].groupby("pipeline")[["TP", "FP", "FN"]].mean()
        xs = np.arange(len(ORDER))
        for k, (m, c) in enumerate([("TP", "#009E73"), ("FN", "#E69F00"), ("FP", "#D55E00")]):
            axi.bar(xs + (k - 1) * 0.27, [max(sub.loc[p, m], 0.5) for p in ORDER], 0.27,
                    color=c, edgecolor="black", lw=0.4, label=m)
        axi.set_yscale("log")
        axi.set_xticks(xs)
        axi.set_xticklabels(ORDER, rotation=35, ha="right")
        axi.set_title(platform)
        axi.axhline(708, ls="--", lw=0.7, color="#333333")
    axes[0].set_ylabel("Taxa (mean of 3 samples, log scale)")
    axes[0].legend(frameon=False, ncol=3, loc="upper left")
    save(fig, "Figure3_TP_FP_FN")

    # ---------------------------------------------------------------- Figure 4
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0))
    for axi, metric in zip(axes, ["recall", "precision", "f1"]):
        mat = (met.assign(col=met.platform + "\n" + met["sample"].str[-1])
               .pivot_table(index="pipeline", columns="col", values=metric)
               .reindex(ORDER))
        im = axi.imshow(mat.values, cmap="viridis", vmin=0, vmax=1, aspect="auto")
        axi.set_xticks(range(mat.shape[1]))
        axi.set_xticklabels(mat.columns, rotation=90, fontsize=6)
        axi.set_yticks(range(len(ORDER)))
        axi.set_yticklabels(ORDER if metric == "recall" else [], fontsize=6.5)
        axi.set_title(metric.capitalize())
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                axi.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5,
                         color="white" if v < 0.6 else "black")
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="Score")
    save(fig, "Figure4_heatmap_pipeline_sample")

    # ---------------------------------------------------------------- Figure 5
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.8))
    for axi, metric in zip(axes, ["recall", "precision", "f1"]):
        for p in ORDER:
            il = met[(met.pipeline == p) & (met.platform == "Illumina")][metric].mean()
            pb = met[(met.pipeline == p) & (met.platform == "PacBio")][metric].mean()
            axi.plot([0, 1], [il, pb], "-o", ms=3.5, lw=1.2, color=COL[p], label=p)
        axi.set_xticks([0, 1])
        axi.set_xticklabels(["Illumina", "PacBio"])
        axi.set_xlim(-0.25, 1.25)
        axi.set_ylim(0, 1.05)
        axi.set_title(metric.capitalize())
    axes[0].set_ylabel("Score (mean of 3 samples)")
    axes[-1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1))
    save(fig, "Figure5_platform_contrast")

    # ---------------------------------------------------------------- Figure 6
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    order_bins = ["<1e-4%", "1e-4-1e-3%", "1e-3-1e-2%", "1e-2-0.1%", ">=0.1%"]
    for axi, platform in zip(axes, ["Illumina", "PacBio"]):
        sub = ab[ab.platform == platform]
        for p in ORDER:
            s = sub[sub.pipeline == p].set_index("abund_bin").reindex(order_bins)
            axi.plot(range(len(order_bins)), s.recall_mean.values, "-o", ms=3,
                     lw=1.1, color=COL[p], label=p)
        axi.set_xticks(range(len(order_bins)))
        axi.set_xticklabels(order_bins, rotation=35, ha="right")
        axi.set_title(platform)
        axi.set_ylim(-0.03, 1.05)
        axi.set_xlabel("Designed relative abundance")
    axes[0].set_ylabel("Recall")
    axes[1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1))
    save(fig, "Figure6_abundance_stratified_recall")

    # ---------------------------------------------------------------- Figure 7
    tp = det[(det.status == "TP") & (det.observed_reads > 0) &
             (det.rel_abundance_pct > 0) & (det["sample"] == "sample_0")]
    fig, axes = plt.subplots(2, 7, figsize=(7.2, 3.2), sharex=True, sharey=True)
    for i, platform in enumerate(["Illumina", "PacBio"]):
        for j, p in enumerate(ORDER):
            axi = axes[i, j]
            g = tp[(tp.platform == platform) & (tp.pipeline == p)]
            axi.scatter(g.rel_abundance_pct, g.observed_reads, s=1.5, alpha=0.35,
                        color=COL[p], edgecolors="none")
            axi.set_xscale("log")
            axi.set_yscale("log")
            r = cor[(cor.platform == platform) & (cor.pipeline == p) &
                    (cor["sample"] == "sample_0")]
            lbl = f"ρ={r.spearman_rho.iloc[0]:.2f}" if len(r) else "n/a"
            axi.text(0.05, 0.92, lbl, transform=axi.transAxes, fontsize=5.5, va="top")
            if i == 0:
                axi.set_title(p, fontsize=5.8)
            if j == 0:
                axi.set_ylabel(f"{platform}\nobserved reads", fontsize=6)
            axi.tick_params(labelsize=5)
    fig.supxlabel("Designed relative abundance (%)", fontsize=8)
    save(fig, "Figure7_abundance_correlation")

    # ---------------------------------------------------------------- Figure 8
    doms = ["Bacteria/Archaea", "Bacteriophage", "Virus (non-phage)"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    for axi, platform in zip(axes, ["Illumina", "PacBio"]):
        mat = (dom[dom.platform == platform]
               .pivot_table(index="pipeline", columns="domain", values="recall_mean")
               .reindex(index=ORDER, columns=doms))
        im = axi.imshow(mat.values, cmap="magma", vmin=0, vmax=1, aspect="auto")
        axi.set_xticks(range(3))
        axi.set_xticklabels([f"{d}\n(n={int(dom[dom.domain == d].n_truth.iloc[0])})"
                             for d in doms], fontsize=6)
        axi.set_yticks(range(len(ORDER)))
        axi.set_yticklabels(ORDER if platform == "Illumina" else [], fontsize=6.5)
        axi.set_title(platform)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                axi.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.5,
                         color="white" if v < 0.55 else "black")
    fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02, label="Recall")
    save(fig, "Figure8_domain_recall")

    # ---------------------------------------------------------------- Figure 9
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    for axi, platform in zip(axes, ["Illumina", "PacBio"]):
        sub = prim[prim.platform == platform].groupby(["pipeline", "min_reads"])[
            ["recall", "precision"]].mean().reset_index()
        for p in ORDER:
            s = sub[sub.pipeline == p].sort_values("min_reads")
            axi.plot(s.recall, s.precision, "-o", ms=3, lw=1.1, color=COL[p], label=p)
            for _, r in s.iterrows():
                axi.annotate(f"≥{int(r.min_reads)}", (r.recall, r.precision),
                             fontsize=4.5, xytext=(2, 2), textcoords="offset points")
        axi.set_xlabel("Recall")
        axi.set_title(platform)
        axi.set_xlim(-0.03, 1.03)
        axi.set_ylim(-0.03, 1.03)
    axes[0].set_ylabel("Precision")
    axes[1].legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1))
    save(fig, "Figure9_read_threshold_tradeoff")

    # --------------------------------------------------------------- Figure 10
    diag = pd.read_csv(INT / "abundance_scale_diagnosis.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(5.2, 3.0))
    xs = np.arange(6)
    lab = comp.platform + "\n" + comp["sample"].str[-1]
    ax.bar(xs, comp.pct_virus, color="#D55E00", edgecolor="black", lw=0.4,
           label="Observed viral read share")
    ax.axhline(12.0, ls="--", lw=1.0, color="#0072B2",
               label="Designed viral genome-copy share (12%)")
    # length-weighted prediction of the viral READ share, per platform
    for platform, x0, x1, col in [("Illumina", -0.45, 2.45, "#009E73"),
                                  ("PacBio", 2.55, 5.45, "#7F2704")]:
        d = diag[diag.platform == platform].groupby("domain").ratio_vs_bacteria.mean()
        ab = diag[diag.platform == platform].groupby("domain").designed_abundance_pct.mean()
        num = sum(ab[k] * d[k] for k in ["Bacteriophage", "Virus"])
        den = num + ab["Bacteria/Archaea"] * d["Bacteria/Archaea"]
        ax.plot([x0, x1], [100 * num / den] * 2, ls=":", lw=1.4, color=col,
                label=f"{platform}: length-weighted prediction")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(lab, fontsize=6)
    ax.set_ylabel("Reads assigned to Viruses (% of total)")
    ax.set_ylim(5e-3, 3e2)
    ax.legend(frameon=False, fontsize=5.5, loc="lower center", ncol=2,
              bbox_to_anchor=(0.5, -0.42))
    save(fig, "Figure10_realised_viral_fraction")
    return 0


if __name__ == "__main__":
    sys.exit(main())
