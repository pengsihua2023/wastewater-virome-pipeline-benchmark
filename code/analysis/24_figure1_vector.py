#!/usr/bin/env python3
"""Vector redraw of the author's PowerPoint Figure 1 (figure-1.png).

Layout, colours and wording follow figure-1.png (1706 x 922 px canvas, y down).
Writes Figure1_study_design_vector.{svg,pdf} (text converted to outlines in the
SVG, embedded TrueType in the PDF) and a 600 dpi PNG for Word/PowerPoint.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = Path(__file__).resolve().parents[2] / "data" / "analysis_708" / "figures_combined"
plt.rcParams.update({"svg.fonttype": "path", "pdf.fonttype": 42, "font.family": "DejaVu Sans"})

W_PX, H_PX = 1706, 922
FIG_W = 7.2
PT = FIG_W * 72 / W_PX                     # points per canvas pixel

BLUE = dict(fc="#E0EEFD", ec="#0B3D7A", tc="#111111", hc="#111111")
TEAL = dict(fc="#E0F2F1", ec="#237D7D", tc="#1A1A1A", hc="#0F6E6E")
ORNG = dict(fc="#FEEAD2", ec="#E48607", tc="#1A1A1A", hc="#111111")
ARROW = "#0F4E83"

fig = plt.figure(figsize=(FIG_W, FIG_W * H_PX / W_PX))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W_PX)
ax.set_ylim(H_PX, 0)
ax.axis("off")
fig.patch.set_facecolor("white")

HEAD, BODY = 26 * PT, 24 * PT


def box(x0, y0, x1, y1, style, head, head_y, body, body_y, line=35):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=16",
                                fc=style["fc"], ec=style["ec"], lw=4 * PT))
    cx = (x0 + x1) / 2
    for i, t in enumerate(head):
        ax.text(cx, head_y + i * line, t, ha="center", va="center", fontsize=HEAD,
                fontweight="bold", color=style["hc"])
    for i, t in enumerate(body):
        ax.text(cx, body_y + i * line, t, ha="center", va="center", fontsize=BODY, color=style["tc"])


def arrow(pts):
    xs, ys = zip(*pts)
    if len(pts) > 2:
        ax.plot(xs[:-1], ys[:-1], color=ARROW, lw=4 * PT, solid_capstyle="butt", solid_joinstyle="miter")
    ax.annotate("", xy=pts[-1], xytext=pts[-2],
                arrowprops=dict(arrowstyle="-|>,head_length=0.55,head_width=0.28", color=ARROW,
                                lw=4 * PT, shrinkA=0, shrinkB=0))


box(36, 92, 387, 306, BLUE, ["Wastewater sample"], 155,
    ["HTP-2024-01-15", "Grimm et al. 2026", "Kraken2 clade reads"], 191)
box(437, 92, 794, 306, BLUE, ["151-member viral", "ground truth"], 155, ["(+ non-viral background)"], 232)
box(843, 92, 1195, 306, BLUE, ["One set of abundance", "files"], 159, ["distribution_0/1/2"], 238)
box(1291, 57, 1670, 206, TEAL, ["CAMISIM 2.0"], 95, ["→ 3 Illumina samples", "~16.5 M read pairs each"], 133)
box(1291, 232, 1670, 380, TEAL, ["PBSIM3 + CCS"], 271, ["→ 3 PacBio HiFi samples", "~479 k reads each"], 311)
box(172, 475, 813, 631, BLUE, ["7 pipelines × 6 datasets"], 534, ["42 scored runs"], 573)
box(893, 475, 1536, 631, BLUE, ["Viral-only scoring by NCBI taxid"], 534, ["at ≥1 / ≥10 / ≥100 reads"], 573)
box(252, 728, 1453, 870, ORNG, ["Viral recall, precision, and F1 for each pipeline"], 779,
    ["Reported separately for each modality"], 820)

arrow([(389, 196), (436, 196)])
arrow([(796, 196), (842, 196)])
ax.plot([1196, 1232], [188, 188], color=ARROW, lw=4 * PT, solid_capstyle="butt")
ax.plot([1232, 1232], [115, 293], color=ARROW, lw=4 * PT, solid_capstyle="butt")
arrow([(1232, 115), (1288, 115)])
arrow([(1232, 293), (1288, 293)])
arrow([(1480, 381), (1480, 416), (520, 416), (520, 466)])
arrow([(815, 553), (888, 553)])
arrow([(1213, 633), (1213, 720)])

for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
    fig.savefig(OUT / f"Figure1_study_design_vector.{ext}", facecolor="white", **kw)
print("wrote", [p.name for p in sorted(OUT.glob("Figure1_study_design_vector.*"))])
