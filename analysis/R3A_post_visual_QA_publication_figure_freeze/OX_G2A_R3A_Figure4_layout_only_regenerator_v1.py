from pathlib import Path
import sys
import csv
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from PIL import Image

(
    design_s,
    global_s,
    components_s,
    subjects_s,
    perm_s,
    outbase_s,
) = sys.argv[1:]

design_path = Path(design_s)
global_path = Path(global_s)
components_path = Path(components_s)
subjects_path = Path(subjects_s)
perm_path = Path(perm_s)
outbase = Path(outbase_s)


def read_tsv(path):
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as f:
        return list(
            csv.DictReader(
                f,
                delimiter="\t",
            )
        )


def fmt_p(p):
    p = float(p)
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}".rstrip("0").rstrip(".")


design = {
    r["property"]: r["value"]
    for r in read_tsv(design_path)
}

required = {
    "canvas_width_mm": "180",
    "font_family": "DejaVu Sans",
    "final_raster_dpi": "600",
    "bbox_policy": "FIXED_CANVAS_NO_TIGHT_DIMENSION_DRIFT",
}

for k, v in required.items():
    if design.get(k) != v:
        raise RuntimeError(
            f"design-system mismatch {k}: "
            f"{design.get(k)!r}"
        )

ANALYTES = [
    "24S_OHC",
    "25_OHC",
    "27_OHC",
    "7B_OHC",
    "7_KC",
]

LABELS = [
    "24S-OHC",
    "25-OHC",
    "27-OHC",
    "7β-OHC",
    "7-KC",
]

BLUE = "#0072B2"
POSITIVE = "#D55E00"
NEUTRAL = "#6F6F6F"
LIGHT = "#E5E5E5"
DARK = "#202020"

WIDTH_MM = 180.0
HEIGHT_MM = 154.0
DPI = 600
MM_TO_IN = 1.0 / 25.4

BASE_FONT = float(design["base_font_pt"])
AXIS_FONT = float(design["axis_label_pt"])
TICK_FONT = float(design["tick_label_pt"])
PANEL_FONT = float(design["panel_label_pt"])
ANNOT_FONT = float(design["annotation_pt"])

matplotlib.rcParams.update({
    "font.family": design["font_family"],
    "font.size": BASE_FONT,
    "axes.labelsize": AXIS_FONT,
    "xtick.labelsize": TICK_FONT,
    "ytick.labelsize": TICK_FONT,
    "axes.linewidth": float(design["axis_line_width_pt"]),
    "lines.linewidth": float(design["data_line_width_pt"]),
    "pdf.fonttype": int(design["pdf_fonttype"]),
    "savefig.facecolor": "white",
    "figure.facecolor": "white",
})


def panel_label(ax, label, x=-0.10):
    ax.text(
        x,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=PANEL_FONT,
        fontweight="bold",
        clip_on=False,
    )


def panel_heading(ax, text):
    ax.text(
        0.0,
        1.03,
        text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=BASE_FONT,
        fontweight="bold",
        clip_on=False,
    )


def clean_axis(ax, grid=False):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if grid:
        ax.grid(
            axis="x",
            color=LIGHT,
            linewidth=0.6,
            zorder=0,
        )


components = read_tsv(
    components_path
)

component_by = {
    r["analyte"]: r
    for r in components
}

subjects = read_tsv(
    subjects_path
)

global_row = read_tsv(
    global_path
)[0]

perm = read_tsv(
    perm_path
)

fig = plt.figure(
    figsize=(
        WIDTH_MM * MM_TO_IN,
        HEIGHT_MM * MM_TO_IN,
    ),
    facecolor="white",
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.08,
    right=0.98,
    bottom=0.09,
    top=0.95,
    wspace=0.30,
    hspace=0.38,
    height_ratios=[
        1.0,
        0.85,
    ],
)

axA = fig.add_subplot(
    gs[0, 0]
)

axB = fig.add_subplot(
    gs[0, 1]
)

axC = fig.add_subplot(
    gs[1, :]
)

# ----------------------------------------------------------------------
# Panel A
# ----------------------------------------------------------------------

vals = np.asarray([
    float(
        component_by[a][
            "hyper_minus_control"
        ]
    )
    for a in ANALYTES
])

y = np.arange(5)

for yi, value in zip(
    y,
    vals,
):
    axA.plot(
        [0, value],
        [yi, yi],
        color=LIGHT,
        linewidth=1.3,
    )

    axA.scatter(
        value,
        yi,
        color=(
            BLUE
            if value < 0
            else POSITIVE
        ),
        s=30,
        zorder=3,
    )

axA.axvline(
    0,
    color=DARK,
    linewidth=0.8,
)

axA.set_yticks(
    y,
    LABELS,
)

axA.invert_yaxis()

axA.set_xlabel(
    "Hyper − control mean Δrank"
)

panel_label(
    axA,
    "A",
)

panel_heading(
    axA,
    "Five-component S2 contrast",
)

clean_axis(
    axA,
    grid=True,
)

# ----------------------------------------------------------------------
# Panel B
# ----------------------------------------------------------------------

matrix = np.asarray([
    [
        float(
            row[
                f"{a}_delta_rank"
            ]
        )
        for a in ANALYTES
    ]
    for row in subjects
])

v = max(
    1e-12,
    float(
        np.max(
            np.abs(matrix)
        )
    ),
)

im = axB.imshow(
    matrix,
    cmap="RdBu_r",
    vmin=-v,
    vmax=v,
    aspect="auto",
    interpolation="nearest",
)

axB.set_xticks(
    range(5),
    LABELS,
    rotation=35,
    ha="right",
)

subject_labels = []

for row in subjects:

    prefix = (
        "H"
        if row["clinical_group"]
        == "HYPERCHOLESTEROLAEMIC"
        else "C"
    )

    subject_labels.append(
        prefix
        + str(
            row["subject_id"]
        )
    )

axB.set_yticks(
    range(
        len(subjects)
    ),
    subject_labels,
)

axB.axhline(
    9.5,
    color=DARK,
    linewidth=1.0,
)

panel_label(
    axB,
    "B",
)

panel_heading(
    axB,
    "Subject longitudinal rank displacement",
)

cax = inset_axes(
    axB,
    width="4%",
    height="60%",
    loc="center right",
    borderpad=1.1,
)

cb = fig.colorbar(
    im,
    cax=cax,
)

cb.set_label(
    "Δrank (B − A)",
    fontsize=ANNOT_FONT,
)

# ----------------------------------------------------------------------
# Panel C
#
# ONLY SCIENTIFIC-FIGURE CHANGE:
# panel label moved inward from x=-0.10 to x=-0.055.
# ----------------------------------------------------------------------

perm_values = np.asarray([
    float(
        row["global_statistic"]
    )
    for row in perm
])

observed = float(
    global_row[
        "observed_statistic"
    ]
)

exact_p = float(
    global_row[
        "exact_p"
    ]
)

axC.hist(
    perm_values,
    bins=55,
    color="#D9D9D9",
    edgecolor="white",
    linewidth=0.35,
)

axC.axvline(
    observed,
    color=DARK,
    linewidth=1.3,
)

axC.text(
    0.98,
    0.96,
    "Exact p="
    + fmt_p(
        exact_p
    ),
    transform=axC.transAxes,
    ha="right",
    va="top",
    fontsize=ANNOT_FONT,
)

axC.set_xlabel(
    "S2 global rank-displacement statistic"
)

axC.set_ylabel(
    "Assignments"
)

clean_axis(
    axC,
    grid=False,
)

panel_label(
    axC,
    "C",
    x=-0.055,
)

panel_heading(
    axC,
    "Complete 10-of-20 label-permutation reference",
)

# ----------------------------------------------------------------------
# Save requested final formats
# ----------------------------------------------------------------------

png = outbase.with_suffix(
    ".png"
)

tiff = outbase.with_suffix(
    ".tiff"
)

pdf = outbase.with_suffix(
    ".pdf"
)

fig.canvas.draw()

fig.savefig(
    png,
    dpi=DPI,
    format="png",
    bbox_inches=None,
)

fig.savefig(
    tiff,
    dpi=DPI,
    format="tiff",
    bbox_inches=None,
    pil_kwargs={
        "compression":
            "tiff_lzw",
    },
)

fig.savefig(
    pdf,
    format="pdf",
    bbox_inches=None,
)

plt.close(fig)

expected_w = round(
    WIDTH_MM
    * MM_TO_IN
    * DPI
)

expected_h = round(
    HEIGHT_MM
    * MM_TO_IN
    * DPI
)

for path in [
    png,
    tiff,
]:
    with Image.open(
        path
    ) as im:
        if (
            abs(
                im.size[0]
                - expected_w
            ) > 2
            or
            abs(
                im.size[1]
                - expected_h
            ) > 2
        ):
            raise RuntimeError(
                f"raster-size mismatch "
                f"{path}: {im.size}"
            )

for path in [
    png,
    tiff,
    pdf,
]:
    if (
        not path.exists()
        or path.stat().st_size <= 0
    ):
        raise RuntimeError(
            f"empty Figure 4 output: "
            f"{path}"
        )

print(
    "FIGURE4_LAYOUT_ONLY_REGENERATION=PASS"
)
print(
    f"FIGURE4_PNG={png}"
)
print(
    f"FIGURE4_TIFF={tiff}"
)
print(
    f"FIGURE4_PDF={pdf}"
)
print(
    "FIGURE4_NEW_P_VALUES_COMPUTED=NO"
)
print(
    "FIGURE4_SCIENTIFIC_RESULT_CHANGED=NO"
)
