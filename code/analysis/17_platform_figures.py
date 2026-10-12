#!/usr/bin/env python3
"""Platform-specific VIRAL figures for the two per-platform papers.

Everything is scored against the 151 non-phage viral ground truth; non-viral
detections and non-viral false positives are excluded. Sources:
interim/viral_sample_metrics.tsv, viral_aggregate_metrics.tsv (script 14),
viral_abundance_recall.tsv, viral_abundance_correlation.tsv, viral_tp_reads.tsv
(script 16).

Five figures per platform, written to figures_illumina/ and figures_pacbio/ as
600-dpi PNG and vector PDF:

  Figure 1  viral recall / precision / F1 by pipeline (per-sample points)
  Figure 2  viral TP / FN / FP counts (log scale)
  Figure 3  viral recall by designed abundance decade
  Figure 4  viral recall-precision trade-off across detection thresholds
  Figure 5  observed reads vs designed abundance for viral true positives
  Figure 6  detection vs the reads each virus actually contributed (Illumina only)
  Figure 7  CLARK: unfiltered output vs its own gamma>=0.2 filter
  Figure 8  viral recall grouped by the input each pipeline was given (PacBio only)
  Figure 8  assembly loss + fastp effect (Illumina only)
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "data"      # all data live under data/
AN = ROOT / "analysis_708"
INT = AN / "interim"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.bbox": "tight",
})

COL = {"TaxProfiler": "#0072B2", "TaxTriage": "#56B4E9", "mag-nf": "#009E73",
       "KrakenMetaReads-nf": "#E69F00", "GOTTCHA2": "#D55E00",
       "CLARK": "#CC79A7", "sourmash": "#7F7F7F"}
# TaxTriage tracks TaxProfiler closely; dashed + open marker keeps both visible
STYLE = {"TaxTriage": dict(ls="--", marker="s", mfc="none")}
BINS = ["<1e-4%", "1e-4-1e-3%", "1e-3-1e-2%", ">=1e-2%"]
BIN_N = {"<1e-4%": 15, "1e-4-1e-3%": 75, "1e-3-1e-2%": 43, ">=1e-2%": 18}
N_TRUTH = 151


def save(fig, outdir, name):
    fig.savefig(outdir / f"{name}.png", dpi=600)
    fig.savefig(outdir / f"{name}.pdf")
    plt.close(fig)
    print("  wrote", name)


def main():
    vm = pd.read_csv(INT / "viral_sample_metrics.tsv", sep="\t")
    vm = vm[vm.ground_truth == "non-phage virus (n=151)"]
    prim = vm[vm.min_reads == 1]
    agg = pd.read_csv(INT / "viral_aggregate_metrics.tsv", sep="\t")
    ab = pd.read_csv(INT / "viral_abundance_recall.tsv", sep="\t")
    cor = pd.read_csv(INT / "viral_abundance_correlation.tsv", sep="\t")
    tpr = pd.read_csv(INT / "viral_tp_reads.tsv", sep="\t")
    dvd = pd.read_csv(INT / "detection_vs_depth.tsv", sep="\t")
    clk = pd.read_csv(INT / "clark_gamma_comparison.tsv", sep="\t")
    bity = pd.read_csv(INT / "viral_recall_by_input_type.tsv", sep="\t")
    asmv = pd.read_csv(INT / "assembly_vs_reads.tsv", sep="\t")
    fst = pd.read_csv(INT / "fastp_effect.tsv", sep="\t")

    for platform, tag in [("Illumina", "illumina"), ("PacBio", "pacbio")]:
        outdir = AN / f"figures_{tag}"
        outdir.mkdir(exist_ok=True)
        print(f"--- {platform} -> {outdir.name}")
        ORDER = list(agg[agg.platform == platform]
                     .sort_values("f1_mean", ascending=False).pipeline)
        sub = prim[prim.platform == platform]
        xs = np.arange(len(ORDER))

        # ---------------------------------------------------------- Figure 1
        fig, ax = plt.subplots(figsize=(6.8, 3.2))
        w = 0.26
        for k, (m, c, hatch) in enumerate([("recall", "#BBBBBB", ""),
                                           ("precision", "#777777", "///"),
                                           ("f1", "#333333", "...")]):
            means = [sub[sub.pipeline == p][m].mean() for p in ORDER]
            ax.bar(xs + (k - 1) * w, means, w, color=c, edgecolor="black", lw=0.4,
                   hatch=hatch, label=["Recall", "Precision", "F1"][k])
            for xi, p in zip(xs + (k - 1) * w, ORDER):
                v = sub[sub.pipeline == p][m].values
                v = v[np.isfinite(v)]
                ax.plot([xi] * len(v), v, "o", ms=2.2, color="#D55E00", zorder=3)
        ax.set_xticks(xs)
        ax.set_xticklabels(ORDER, rotation=30, ha="right")
        ax.set_ylabel(f"Score (viral species, n = {N_TRUTH})")
        ax.set_ylim(0, 1.05)
        ax.set_title(f"{platform}: viral detection at ≥1 clade read")
        ax.legend(frameon=False, ncol=3, loc="upper right")
        save(fig, outdir, "Figure1_viral_recall_precision_f1")

        # ---------------------------------------------------------- Figure 2
        fig, ax = plt.subplots(figsize=(6.8, 3.2))
        counts = sub.groupby("pipeline")[["TP", "FN", "FP"]].mean()
        for k, (m, c) in enumerate([("TP", "#009E73"), ("FN", "#E69F00"),
                                    ("FP", "#D55E00")]):
            vals = [max(counts.loc[p, m], 0.5) for p in ORDER]
            bars = ax.bar(xs + (k - 1) * 0.27, vals, 0.27, color=c,
                          edgecolor="black", lw=0.4, label=m)
            for b, p in zip(bars, ORDER):
                ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.18,
                        f"{counts.loc[p, m]:.0f}", ha="center", fontsize=5, rotation=90)
        ax.set_yscale("log")
        ax.set_ylim(0.4, 3e4)
        ax.axhline(N_TRUTH, ls="--", lw=0.7, color="#333333")
        ax.text(len(ORDER) - 0.45, N_TRUTH * 1.25, f"{N_TRUTH} viruses",
                fontsize=6, color="#333333")
        ax.set_xticks(xs)
        ax.set_xticklabels(ORDER, rotation=30, ha="right")
        ax.set_ylabel("Viral taxa (mean of 3 samples, log scale)")
        ax.set_title(f"{platform}: viral true positives, false negatives, false positives")
        ax.legend(frameon=False, ncol=3, loc="upper left")
        save(fig, outdir, "Figure2_viral_TP_FN_FP")

        # ---------------------------------------------------------- Figure 3
        fig, ax = plt.subplots(figsize=(5.8, 3.2))
        s = ab[ab.platform == platform]
        for p in ORDER:
            y = s[s.pipeline == p].set_index("bin").reindex(BINS).recall_mean
            st = STYLE.get(p, dict(ls="-", marker="o"))
            ax.plot(range(len(BINS)), y.values, ms=3.8, lw=1.2,
                    color=COL[p], label=p, **st)
        ax.set_xticks(range(len(BINS)))
        ax.set_xticklabels([f"{b}\n(n={BIN_N[b]})" for b in BINS], fontsize=6.5)
        ax.set_xlabel("Designed relative abundance of the virus")
        ax.set_ylabel("Viral recall")
        ax.set_ylim(-0.03, 1.05)
        ax.set_title(f"{platform}: viral recall by designed abundance")
        ax.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1))
        save(fig, outdir, "Figure3_viral_abundance_recall")

        # ---------------------------------------------------------- Figure 4
        fig, ax = plt.subplots(figsize=(5.8, 3.4))
        sweep = (vm[vm.platform == platform]
                 .groupby(["pipeline", "min_reads"])[["recall", "precision"]]
                 .mean().reset_index())
        for p in ORDER:
            g = sweep[sweep.pipeline == p].sort_values("min_reads")
            g = g[np.isfinite(g.precision) & np.isfinite(g.recall)]
            if g.empty:
                continue
            st = STYLE.get(p, dict(ls="-", marker="o"))
            ax.plot(g.recall, g.precision, ms=4.5, lw=1.3, color=COL[p], label=p, **st)
            for _, r in g.iterrows():
                ax.annotate(f"≥{int(r.min_reads)}", (r.recall, r.precision),
                            fontsize=5, xytext=(3, 3), textcoords="offset points",
                            color=COL[p])
        ax.set_xlabel("Viral recall")
        ax.set_ylabel("Viral precision")
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_title(f"{platform}: viral recall–precision across thresholds")
        ax.legend(frameon=False, fontsize=6, loc="upper left", bbox_to_anchor=(1.02, 1))
        save(fig, outdir, "Figure4_viral_threshold_tradeoff")

        # ---------------------------------------------------------- Figure 5
        g0 = tpr[(tpr.platform == platform) & (tpr["sample"] == "sample_0")]
        n = len(ORDER)
        fig, axes = plt.subplots(1, n, figsize=(7.2, 2.3), sharex=True, sharey=True)
        for axi, p in zip(np.atleast_1d(axes), ORDER):
            g = g0[g0.pipeline == p]
            if len(g):
                axi.scatter(g.abund, g.reads, s=4, alpha=0.45,
                            color=COL[p], edgecolors="none")
            r = cor[(cor.platform == platform) & (cor.pipeline == p) &
                    (cor["sample"] == "sample_0")]
            if len(r) and np.isfinite(r.spearman_rho.iloc[0]):
                lbl = f"ρ={r.spearman_rho.iloc[0]:.2f}"
            else:
                lbl = "no detections" if not len(g) else "n<10"
            axi.set_xscale("log")
            axi.set_yscale("log")
            axi.set_title(f"{p}\n{lbl}", fontsize=5.8)
            axi.tick_params(labelsize=5)
        np.atleast_1d(axes)[0].set_ylabel("Observed reads", fontsize=7)
        fig.supxlabel("Designed relative abundance of the virus (%)", fontsize=8)
        save(fig, outdir, "Figure5_viral_abundance_correlation")

        # ---------------------------------------------------------- Figure 6
        # how much of the recall shortfall is simply absent signal
        if platform == "Illumina":
            dd = dvd[dvd.platform == platform]
            order = ["<0.5", "0.5-1", "1-2", "2-5", "5-20", ">20"]
            fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
            g = (dd.groupby(["pipeline", "sample", "depth_bin"])
                 .agg(n=("t", "size"), det=("detected", "sum")).reset_index())
            g["r"] = g.det / g.n
            m = g.groupby(["pipeline", "depth_bin"]).r.mean().unstack().reindex(
                columns=order)
            nb = g.groupby("depth_bin").n.first().reindex(order)
            for p_ in ORDER:
                st = STYLE.get(p_, dict(ls="-", marker="o"))
                axes[0].plot(range(len(order)), m.loc[p_].values, ms=3.5, lw=1.2,
                             color=COL[p_], label=p_, **st)
            axes[0].set_xticks(range(len(order)))
            # NOTE (2026-09, regenerated data): the low-expected-reads bins
            # (<0.5, 0.5-1, 1-2, 2-5) are now empty -- every ground-truth virus
            # falls in 5-20 or >20 once the CAMISIM abundance-unit bug is fixed.
            # This figure's whole premise (recall degrading at near-zero expected
            # depth) no longer has data to show; it needs a redesigned bin scheme
            # for V4 rather than just this crash fix. Guarding here only so the
            # rest of the script can still run.
            axes[0].set_xticklabels(
                [f"{b}\n(n={0 if pd.isna(nb[b]) else int(nb[b])})" for b in order],
                fontsize=6)
            axes[0].set_xlabel("Reads the virus was expected to contribute")
            axes[0].set_ylabel("Viral recall")
            axes[0].set_ylim(-0.03, 1.05)
            axes[0].set_title("Recall vs available signal")

            thrs = [0, 1, 2, 5]
            xs2 = np.arange(len(thrs))
            for p_ in ORDER:
                ys, ns = [], []
                for t_ in thrs:
                    sub = dd[(dd.pipeline == p_) & (dd.expected_reads >= t_)]
                    gg = sub.groupby("sample").agg(n=("t", "size"), det=("detected", "sum"))
                    ys.append((gg.det / gg.n).mean())
                    ns.append(int(gg.n.iloc[0]))
                st = STYLE.get(p_, dict(ls="-", marker="o"))
                axes[1].plot(xs2, ys, ms=4, lw=1.3, color=COL[p_], label=p_, **st)
            axes[1].set_xticks(xs2)
            axes[1].set_xticklabels([f"all\n(n=151)", f"≥1\n(n={ns[1]})",
                                     f"≥2\n(n={ns[2]})", f"≥5\n(n={ns[3]})"], fontsize=6.5)
            axes[1].set_xlabel("Restricted to viruses with ≥ N expected reads")
            axes[1].set_ylim(-0.03, 1.05)
            axes[1].set_title("Recall on viruses that had signal")
            axes[1].legend(frameon=False, fontsize=6, loc="upper left",
                           bbox_to_anchor=(1.02, 1))
            save(fig, outdir, "Figure6_recall_vs_available_signal")

        # ---------------------------------------------------------- Figure 7
        # CLARK carries its confidence in gamma, not in read count
        cl = clk[(clk.platform == platform) & (clk.min_reads == 1)].set_index("stream")
        fig, ax = plt.subplots(figsize=(4.6, 3.0))
        xs3 = np.arange(3)
        w3 = 0.36
        for k, (stream, col, hatch) in enumerate(
                [("unfiltered", "#CC79A7", ""), ("gamma>=0.2", "#0072B2", "///")]):
            vals = [cl.loc[stream, "recall"], cl.loc[stream, "precision"],
                    cl.loc[stream, "f1"]]
            ax.bar(xs3 + (k - 0.5) * w3, vals, w3, color=col, edgecolor="black",
                   lw=0.4, hatch=hatch,
                   label=("unfiltered" if k == 0 else "γ ≥ 0.2 (CLARK's own filter)"))
            for x_, v_ in zip(xs3 + (k - 0.5) * w3, vals):
                ax.text(x_, v_ + 0.012, f"{v_:.3f}", ha="center", fontsize=5.5)
        ax.set_xticks(xs3)
        ax.set_xticklabels(["Recall", "Precision", "F1"])
        ax.set_ylim(0, max(0.6, cl["recall"].max() * 1.25))
        ax.set_ylabel("Viral score (n = 151)")
        ax.set_title(f"{platform}: CLARK, unfiltered vs γ ≥ 0.2\n"
                     f"false positives {cl.loc['unfiltered','FP']:.0f} → "
                     f"{cl.loc['gamma>=0.2','FP']:.0f}")
        ax.legend(frameon=False, fontsize=6)
        save(fig, outdir, ("Figure7_clark_gamma_filter" if platform == "Illumina"
                           else "Figure6_clark_gamma_filter"))

        # ---------------------------------------------------------- Figure 8
        # on PacBio the pipelines split by the input they were given, not by
        # classifier family
        if platform == "PacBio":
            b = bity[bity.platform == platform].copy()
            b = b.sort_values(["input_type", "recall_mean"], ascending=[True, False])
            fig, ax = plt.subplots(figsize=(5.6, 3.2))
            xs4 = np.arange(len(b))
            cols = ["#D55E00" if t == "contigs" else "#0072B2" for t in b.input_type]
            bars = ax.bar(xs4, b.recall_mean, 0.68, color=cols,
                          edgecolor="black", lw=0.5)
            for x_, v_ in zip(xs4, b.recall_mean):
                ax.text(x_, v_ + 0.02, f"{v_:.3f}", ha="center", fontsize=6.5)
            ax.set_xticks(xs4)
            ax.set_xticklabels(b.pipeline, rotation=35, ha="right", fontsize=7)
            ax.set_ylabel("Viral recall (n = 151)")
            ax.set_ylim(0, 1.05)
            nc = (b.input_type == "contigs").sum()
            ax.axvline(nc - 0.5, ls="--", lw=0.9, color="#333333")
            ax.text((nc - 1) / 2, 0.97, "given contigs", ha="center", fontsize=7,
                    color="#D55E00")
            ax.text(nc + (len(b) - nc - 1) / 2, 0.97, "given reads", ha="center",
                    fontsize=7, color="#0072B2")
            ax.set_title("PacBio: viral recall tracks the input, not the classifier")
            save(fig, outdir, "Figure7_recall_by_input_type")

        # ---------------------------------------------------------- Figure 8
        # two controlled comparisons: what assembly costs, and what one QC step costs
        if platform == "Illumina":
            fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
            a = asmv[asmv.platform == "Illumina"].set_index("evidence")
            names = ["reads (nf-core/taxprofiler)", "contigs (metaSPAdes)",
                     "contigs (MEGAHIT)"]
            short = ["reads", "metaSPAdes\ncontigs", "MEGAHIT\ncontigs"]
            vals = [a.loc[n, "gt_viruses_hit"] for n in names]
            cols = ["#0072B2", "#D55E00", "#D55E00"]
            bars = axes[0].bar(range(3), vals, 0.6, color=cols, edgecolor="black", lw=0.5)
            for x_, v_ in zip(range(3), vals):
                axes[0].text(x_, v_ + 0.8, f"{v_:.1f}", ha="center", fontsize=7)
            axes[0].set_xticks(range(3))
            axes[0].set_xticklabels(short, fontsize=7)
            axes[0].set_ylabel("Ground-truth viruses recovered (of 151)")
            axes[0].set_ylim(0, max(vals) * 1.25)
            axes[0].set_title("Assembly discards viral signal\n(same Kraken2, same database)")

            m = fst.mean(numeric_only=True)
            labels = ["Viral TP", "Viral FP"]
            tv = [m["taxprofiler_TP"], m["taxprofiler_FP"]]
            mv = [m["mag-nf_TP"], m["mag-nf_FP"]]
            xs5 = np.arange(2)
            axes[1].bar(xs5 - 0.19, tv, 0.38, color="#0072B2", edgecolor="black",
                        lw=0.5, label="nf-core/taxprofiler")
            axes[1].bar(xs5 + 0.19, mv, 0.38, color="#009E73", edgecolor="black",
                        lw=0.5, hatch="///", label="nf-core/mag (+ fastp)")
            for x_, v_ in zip(xs5 - 0.19, tv):
                axes[1].text(x_, v_ + 8, f"{v_:.0f}", ha="center", fontsize=7)
            for x_, v_ in zip(xs5 + 0.19, mv):
                axes[1].text(x_, v_ + 8, f"{v_:.0f}", ha="center", fontsize=7)
            axes[1].set_xticks(xs5)
            axes[1].set_xticklabels(labels)
            axes[1].set_ylabel("Viral taxa")
            axes[1].set_ylim(0, max(mv) * 1.3)
            axes[1].set_title("One QC step doubles false positives\n(same reads, same database)")
            axes[1].legend(frameon=False, fontsize=6.5, loc="upper left")
            save(fig, outdir, "Figure8_assembly_and_preprocessing")

        # remove figures from the previous whole-community figure set
        for stale in ["Figure8_recall_by_input_type", "Figure1_recall_precision_f1", "Figure2_TP_FN_FP",
                      "Figure3_domain_recall", "Figure4_abundance_stratified_recall",
                      "Figure5_threshold_tradeoff", "Figure6_abundance_correlation"]:
            for ext in (".png", ".pdf"):
                f = outdir / f"{stale}{ext}"
                if f.exists():
                    f.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
