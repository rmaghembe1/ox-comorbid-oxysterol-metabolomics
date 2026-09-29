from pathlib import Path
import sys
import csv
import math
import hashlib
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import (
    FancyBboxPatch,
    FancyArrowPatch,
    Circle,
)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from PIL import Image


# ==========================================================================
# Arguments
# ==========================================================================

(
    design_summary_s,
    p1_subject_s,
    p3_subject_s,
    source_manifest_s,
    design_system_s,
    panel_spec_s,
    no_recalc_s,
    expected_outputs_s,
    r2y_receipt_s,
    r2y_freeze_s,
    claims_s,
    limitations_s,
    figdir_s,
    out_s,
    gov_s,
    expected_design_summary,
    expected_p1_subject,
    expected_p3_subject,
    expected_source_manifest,
    expected_design_system,
    expected_panel_spec,
    expected_no_recalc,
    expected_outputs,
    expected_r2y_receipt,
    expected_r2y_freeze,
    expected_claims,
    expected_limitations,
) = sys.argv[1:]

design_summary = Path(design_summary_s)
p1_subject = Path(p1_subject_s)
p3_subject = Path(p3_subject_s)
source_manifest = Path(source_manifest_s)
design_system = Path(design_system_s)
panel_spec = Path(panel_spec_s)
no_recalc = Path(no_recalc_s)
expected_outputs_path = Path(expected_outputs_s)
r2y_receipt = Path(r2y_receipt_s)
r2y_freeze = Path(r2y_freeze_s)
claims = Path(claims_s)
limitations = Path(limitations_s)
figdir = Path(figdir_s)
out = Path(out_s)
gov = Path(gov_s)

expected_binding = {
    "DESIGN_SUMMARY": expected_design_summary,
    "P1_SUBJECT": expected_p1_subject,
    "P3_SUBJECT": expected_p3_subject,
    "SOURCE_MANIFEST": expected_source_manifest,
    "DESIGN_SYSTEM": expected_design_system,
    "PANEL_SPEC": expected_panel_spec,
    "NO_RECALC": expected_no_recalc,
    "EXPECTED_OUTPUTS": expected_outputs,
    "R2Y_RECEIPT": expected_r2y_receipt,
    "R2Y_FREEZE": expected_r2y_freeze,
    "CLAIMS": expected_claims,
    "LIMITATIONS": expected_limitations,
}

binding_paths = {
    "DESIGN_SUMMARY": design_summary,
    "P1_SUBJECT": p1_subject,
    "P3_SUBJECT": p3_subject,
    "SOURCE_MANIFEST": source_manifest,
    "DESIGN_SYSTEM": design_system,
    "PANEL_SPEC": panel_spec,
    "NO_RECALC": no_recalc,
    "EXPECTED_OUTPUTS": expected_outputs_path,
    "R2Y_RECEIPT": r2y_receipt,
    "R2Y_FREEZE": r2y_freeze,
    "CLAIMS": claims,
    "LIMITATIONS": limitations,
}


# ==========================================================================
# Utilities
# ==========================================================================

def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b"",
        ):
            h.update(block)

    return h.hexdigest()


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


def parse_receipt(path):

    d = {}

    for raw in path.read_text(
        encoding="utf-8",
        errors="strict",
    ).splitlines():

        if "=" not in raw:
            continue

        k, v = raw.split(
            "=",
            1,
        )

        d[k.strip()] = v.strip()

    return d


def as_float(x):

    return float(
        str(x).strip()
    )


def can_float(x):

    try:
        float(
            str(x).strip()
        )
        return True

    except Exception:
        return False


def first_existing(row, names):

    for name in names:

        if name in row:
            return row[name]

    raise KeyError(
        "none of fields found: "
        + "|".join(names)
    )


def fmt_p(p):

    p = float(p)

    if p < 0.001:
        return f"{p:.2e}"

    return f"{p:.4f}".rstrip("0").rstrip(".")


def pretty_analyte(x):

    mapping = {
        "24S_OHC": "24S-OHC",
        "25_OHC": "25-OHC",
        "27_OHC": "27-OHC",
        "7B_OHC": "7β-OHC",
        "7beta_OHC": "7β-OHC",
        "7_KC": "7-KC",
    }

    return mapping.get(
        x,
        x.replace("_", "-"),
    )


ANALYTES = [
    "24S_OHC",
    "25_OHC",
    "27_OHC",
    "7B_OHC",
    "7_KC",
]

ANALYTE_LABELS = [
    pretty_analyte(x)
    for x in ANALYTES
]


# ==========================================================================
# 1. Frozen R2Y binding
# ==========================================================================

print("=" * 112)
print("OX-COMORBID G2A R2Z PUBLICATION FIGURE GENERATION")
print("VISUALIZATION ONLY / SIX MAIN FIGURES")
print("NO NEW P-VALUES / NO NEW HYPOTHESIS TESTS / NO NEW INFERENTIAL FAMILIES")
print("=" * 112)

print()
print("----- R2Y BINDING -----")

for label, path in binding_paths.items():

    if not path.exists():

        print(
            f"REFUSE_{label}_MISSING={path}"
        )

        raise SystemExit(2)

    actual = sha256(
        path
    )

    print(
        f"{label}={path}"
    )
    print(
        f"{label}_SHA256={actual}"
    )

    if (
        actual
        !=
        expected_binding[label]
    ):

        print(
            f"REFUSE_{label}_SHA256_MISMATCH"
        )

        raise SystemExit(3)

print("R2Y_BINDING=PASS")


# ==========================================================================
# 2. R2Y authorization state
# ==========================================================================

receipt_state = parse_receipt(
    r2y_receipt
)

required_states = {
    "R2X_STATE": "PASS",
    "GRANULAR_FIGURE_SOURCE_BINDING":
        "PASS",
    "FIGURE_ARCHITECTURE_IDENTITY":
        "PASS",
    "FIGURE_DEPENDENCIES":
        "PASS",
    "P1_SUBJECT_DIRECTION_DECOMPOSITION_IDENTITY":
        "PASS",
    "P3_SUBJECT_DISPLACEMENT_DECOMPOSITION_IDENTITY":
        "PASS",
    "FIGURE1_STUDY_DESIGN_SUMMARY":
        "PASS",
    "PUBLICATION_FIGURE_COUNT":
        "6",
    "PUBLICATION_OUTPUT_FORMATS":
        "TIFF_600DPI|SVG|PDF",
    "QA_OUTPUT_FORMAT":
        "PNG_200DPI",
    "NEW_INFERENTIAL_FAMILIES":
        "HOLD",
    "NEW_P_VALUES_COMPUTED":
        "NO",
    "NEW_HYPOTHESIS_TESTS_EXECUTED":
        "NO",
    "ABSOLUTE_CONCENTRATION_INTERPRETATION":
        "NOT_AUTHORIZED",
    "CONCENTRATION_FOLD_CHANGE_INTERPRETATION":
        "NOT_AUTHORIZED",
    "HOMEOSTASIS_RESTORATION_CLAIM":
        "NOT_AUTHORIZED",
    "CAUSAL_TREATMENT_EFFECT_INTERPRETATION":
        "NOT_AUTHORIZED",
    "PUBLICATION_FIGURES_GENERATED":
        "NO",
    "FIGURE_IMPLEMENTATION":
        "AUTHORIZED_FOR_NEXT_GATE",
}

print()
print("----- R2Y AUTHORIZATION STATE -----")

for key, required in required_states.items():

    actual = receipt_state.get(
        key
    )

    print(
        f"{key}={actual}"
    )

    if actual != required:

        print(
            f"REFUSE_R2Y_STATE="
            f"{key}:{actual}"
        )

        raise SystemExit(4)

print("R2Y_AUTHORIZATION_STATE=PASS")


# ==========================================================================
# 3. Source manifest revalidation
# ==========================================================================

manifest_rows = read_tsv(
    source_manifest
)

sources = {}

print()
print("----- FIGURE SOURCE REVALIDATION -----")

for row in manifest_rows:

    role = row[
        "source_role"
    ]

    path = Path(
        row["path"]
    )

    expected_hash = row[
        "sha256"
    ]

    if not path.exists():

        print(
            f"REFUSE_SOURCE_MISSING="
            f"{role}:{path}"
        )

        raise SystemExit(5)

    actual_hash = sha256(
        path
    )

    if actual_hash != expected_hash:

        print(
            f"REFUSE_SOURCE_HASH="
            f"{role}"
        )

        raise SystemExit(6)

    if role in sources:

        print(
            f"REFUSE_DUPLICATE_SOURCE_ROLE="
            f"{role}"
        )

        raise SystemExit(7)

    sources[
        role
    ] = path

    print(
        f"SOURCE={role} "
        f"SHA256={actual_hash} "
        f"STATUS=PASS"
    )

print(
    "FIGURE_SOURCE_REVALIDATION=PASS"
)


# ==========================================================================
# 4. Expected-output identity
# ==========================================================================

expected_rows = read_tsv(
    expected_outputs_path
)

if len(
    expected_rows
) != 6:

    print(
        f"REFUSE_EXPECTED_OUTPUT_COUNT="
        f"{len(expected_rows)}"
    )

    raise SystemExit(8)

expected_stems = [
    r["figure"]
    for r in expected_rows
]

required_stems = [
    "Figure_1_provenance_and_analysis_architecture",
    "Figure_2_primary_longitudinal_topology",
    "Figure_3_rank_network_rewiring",
    "Figure_4_between_block_rank_displacement",
    "Figure_5_control_anchored_multivariate_departure",
    "Figure_6_response_scale_sensitivity",
]

if expected_stems != required_stems:

    print(
        "REFUSE_EXPECTED_FIGURE_STEMS"
    )

    raise SystemExit(9)

print()
print("EXPECTED_FIGURE_COUNT=6")
print("EXPECTED_OUTPUT_IDENTITY=PASS")


# ==========================================================================
# 5. Frozen visual design
# ==========================================================================

design_rows = read_tsv(
    design_system
)

design = {
    r["property"]: r["value"]
    for r in design_rows
}

required_design = {
    "canvas_width_mm": "180",
    "maximum_canvas_height_mm": "205",
    "font_family": "DejaVu Sans",
    "final_raster_dpi": "600",
    "final_formats": "TIFF|SVG|PDF",
    "qa_format": "PNG",
    "qa_dpi": "200",
    "bbox_policy":
        "FIXED_CANVAS_NO_TIGHT_DIMENSION_DRIFT",
}

for key, value in required_design.items():

    if design.get(key) != value:

        print(
            f"REFUSE_DESIGN_SYSTEM="
            f"{key}:{design.get(key)}"
        )

        raise SystemExit(10)

WIDTH_MM = float(
    design["canvas_width_mm"]
)

MAX_HEIGHT_MM = float(
    design["maximum_canvas_height_mm"]
)

MM_TO_IN = (
    1.0 / 25.4
)

WIDTH_IN = (
    WIDTH_MM * MM_TO_IN
)

BASE_FONT = float(
    design["base_font_pt"]
)

AXIS_FONT = float(
    design["axis_label_pt"]
)

TICK_FONT = float(
    design["tick_label_pt"]
)

PANEL_FONT = float(
    design["panel_label_pt"]
)

ANNOT_FONT = float(
    design["annotation_pt"]
)

matplotlib.rcParams.update({
    "font.family":
        design["font_family"],
    "font.size":
        BASE_FONT,
    "axes.labelsize":
        AXIS_FONT,
    "xtick.labelsize":
        TICK_FONT,
    "ytick.labelsize":
        TICK_FONT,
    "axes.linewidth":
        float(
            design[
                "axis_line_width_pt"
            ]
        ),
    "lines.linewidth":
        float(
            design[
                "data_line_width_pt"
            ]
        ),
    "svg.fonttype":
        design[
            "svg_fonttype"
        ],
    "pdf.fonttype":
        int(
            design[
                "pdf_fonttype"
            ]
        ),
    "savefig.facecolor":
        "white",
    "figure.facecolor":
        "white",
})


# ==========================================================================
# 6. Consistent colorblind-safe role palette
# ==========================================================================

HYPER = "#0072B2"
CONTROL = "#D55E00"
NEGATIVE = "#0072B2"
POSITIVE = "#D55E00"
SUPPORTED = "#009E73"
NOT_SUPPORTED = "#B8B8B8"
NEUTRAL = "#6F6F6F"
LIGHT = "#E5E5E5"
DARK = "#202020"


# ==========================================================================
# Plot helpers
# ==========================================================================

def panel_label(
    ax,
    label,
):

    ax.text(
        -0.10,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=PANEL_FONT,
        fontweight="bold",
        clip_on=False,
    )


def panel_heading(
    ax,
    text,
):

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


def clean_axis(
    ax,
    grid=False,
):

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)

    if grid:

        ax.grid(
            axis="x",
            color=LIGHT,
            linewidth=0.6,
            zorder=0,
        )


def new_figure(
    height_mm,
):

    if height_mm > MAX_HEIGHT_MM:

        raise RuntimeError(
            "figure height exceeds "
            "frozen maximum"
        )

    fig = plt.figure(
        figsize=(
            WIDTH_IN,
            height_mm * MM_TO_IN,
        ),
        facecolor="white",
    )

    return fig


def save_figure(
    fig,
    stem,
    height_mm,
):

    publication_dpi = int(
        design[
            "final_raster_dpi"
        ]
    )

    qa_dpi = int(
        design[
            "qa_dpi"
        ]
    )

    tiff = figdir / (
        stem + ".tiff"
    )

    svg = figdir / (
        stem + ".svg"
    )

    pdf = figdir / (
        stem + ".pdf"
    )

    png = figdir / (
        stem + "_QA.png"
    )

    fig.canvas.draw()

    fig.savefig(
        tiff,
        dpi=publication_dpi,
        format="tiff",
        bbox_inches=None,
        pil_kwargs={
            "compression":
                "tiff_lzw",
        },
    )

    fig.savefig(
        svg,
        format="svg",
        bbox_inches=None,
    )

    fig.savefig(
        pdf,
        format="pdf",
        bbox_inches=None,
    )

    fig.savefig(
        png,
        dpi=qa_dpi,
        format="png",
        bbox_inches=None,
    )

    plt.close(
        fig
    )

    expected_pub_w = round(
        WIDTH_IN
        * publication_dpi
    )

    expected_pub_h = round(
        height_mm
        * MM_TO_IN
        * publication_dpi
    )

    expected_qa_w = round(
        WIDTH_IN
        * qa_dpi
    )

    expected_qa_h = round(
        height_mm
        * MM_TO_IN
        * qa_dpi
    )

    with Image.open(
        tiff
    ) as im:

        pub_size = im.size
        compression = (
            getattr(
                im,
                "info",
                {},
            ).get(
                "compression",
                ""
            )
        )

    with Image.open(
        png
    ) as im:

        qa_size = im.size

    if abs(
        pub_size[0]
        - expected_pub_w
    ) > 2 or abs(
        pub_size[1]
        - expected_pub_h
    ) > 2:

        raise RuntimeError(
            f"publication raster size "
            f"mismatch for {stem}: "
            f"{pub_size} versus "
            f"{expected_pub_w}x"
            f"{expected_pub_h}"
        )

    if abs(
        qa_size[0]
        - expected_qa_w
    ) > 2 or abs(
        qa_size[1]
        - expected_qa_h
    ) > 2:

        raise RuntimeError(
            f"QA raster size mismatch "
            f"for {stem}: {qa_size}"
        )

    outputs = [
        (
            "TIFF",
            tiff,
            pub_size[0],
            pub_size[1],
            publication_dpi,
            compression,
        ),
        (
            "SVG",
            svg,
            "",
            "",
            "",
            "",
        ),
        (
            "PDF",
            pdf,
            "",
            "",
            "",
            "",
        ),
        (
            "PNG_QA",
            png,
            qa_size[0],
            qa_size[1],
            qa_dpi,
            "",
        ),
    ]

    for _, path, *_ in outputs:

        if (
            not path.exists()
            or
            path.stat().st_size <= 0
        ):

            raise RuntimeError(
                f"empty output: {path}"
            )

    print(
        f"FIGURE={stem} "
        f"STATUS=GENERATED "
        f"PUB_PIXELS="
        f"{pub_size[0]}x{pub_size[1]} "
        f"QA_PIXELS="
        f"{qa_size[0]}x{qa_size[1]}"
    )

    return outputs


def matrix_from_tsv(
    path,
):

    rows = read_tsv(
        path
    )

    if not rows:

        raise RuntimeError(
            f"empty matrix table: {path}"
        )

    keys = list(
        rows[0].keys()
    )

    row_key = None

    for key in keys:

        values = [
            r[key]
            for r in rows
        ]

        if not all(
            can_float(v)
            for v in values
        ):

            row_key = key
            break

    if row_key is None:

        row_key = keys[0]

    row_labels = [
        r[row_key]
        for r in rows
    ]

    numeric_cols = [
        key
        for key in keys
        if (
            key != row_key
            and
            all(
                can_float(
                    r[key]
                )
                for r in rows
            )
        )
    ]

    if all(
        label in numeric_cols
        for label in row_labels
    ):

        cols = row_labels

    else:

        if len(
            numeric_cols
        ) < len(
            rows
        ):

            raise RuntimeError(
                f"cannot identify square "
                f"numeric matrix: {path}"
            )

        cols = numeric_cols[
            :len(rows)
        ]

    arr = np.asarray([
        [
            float(
                row[col]
            )
            for col in cols
        ]
        for row in rows
    ])

    if (
        arr.shape[0]
        !=
        arr.shape[1]
    ):

        raise RuntimeError(
            f"matrix not square: "
            f"{path}:{arr.shape}"
        )

    return (
        row_labels,
        cols,
        arr,
    )


def find_family_stat_col(
    rows,
    family,
):

    keys = list(
        rows[0].keys()
    )

    family_l = (
        family.lower()
    )

    preferred_tokens = [
        "statistic",
        "t_dir",
        "t_net",
        "t_centroid",
        "t_",
        "stat",
    ]

    candidates = [
        key
        for key in keys
        if (
            family_l
            in key.lower()
            and
            all(
                can_float(
                    r[key]
                )
                for r in rows
            )
        )
    ]

    for token in preferred_tokens:

        matches = [
            key
            for key in candidates
            if token in key.lower()
        ]

        if len(
            matches
        ) == 1:

            return matches[0]

    if len(
        candidates
    ) == 1:

        return candidates[0]

    raise RuntimeError(
        f"cannot identify "
        f"{family} permutation "
        f"statistic column; "
        f"candidates={candidates}"
    )


def histogram_panel(
    ax,
    values,
    observed,
    p_text,
    xlabel,
    bins=40,
):

    ax.hist(
        values,
        bins=bins,
        color="#D9D9D9",
        edgecolor="white",
        linewidth=0.35,
    )

    ax.axvline(
        observed,
        color=DARK,
        linewidth=1.3,
    )

    ax.text(
        0.98,
        0.96,
        p_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=ANNOT_FONT,
    )

    ax.set_xlabel(
        xlabel
    )

    ax.set_ylabel(
        "Assignments"
    )

    clean_axis(
        ax,
        grid=False,
    )


def draw_network(
    ax,
    matrix,
    labels,
    heading,
):

    n = len(
        labels
    )

    angles = (
        np.linspace(
            np.pi / 2,
            np.pi / 2
            - 2 * np.pi,
            n,
            endpoint=False,
        )
    )

    positions = np.column_stack([
        np.cos(
            angles
        ),
        np.sin(
            angles
        ),
    ])

    for i in range(n):

        for j in range(
            i + 1,
            n,
        ):

            rho = float(
                matrix[i, j]
            )

            color = (
                POSITIVE
                if rho >= 0
                else NEGATIVE
            )

            width = (
                0.55
                + 2.6
                * abs(rho)
            )

            alpha = (
                0.25
                + 0.65
                * abs(rho)
            )

            ax.plot(
                [
                    positions[i, 0],
                    positions[j, 0],
                ],
                [
                    positions[i, 1],
                    positions[j, 1],
                ],
                color=color,
                linewidth=width,
                alpha=min(
                    0.95,
                    alpha,
                ),
                zorder=1,
            )

    ax.scatter(
        positions[:, 0],
        positions[:, 1],
        s=430,
        facecolor="white",
        edgecolor=DARK,
        linewidth=1.0,
        zorder=3,
    )

    for i, label in enumerate(
        labels
    ):

        ax.text(
            positions[i, 0],
            positions[i, 1],
            pretty_analyte(
                label
            ),
            ha="center",
            va="center",
            fontsize=ANNOT_FONT,
            fontweight="bold",
            zorder=4,
        )

    ax.set_xlim(
        -1.40,
        1.40,
    )

    ax.set_ylim(
        -1.30,
        1.30,
    )

    ax.set_aspect(
        "equal"
    )

    ax.axis(
        "off"
    )

    panel_heading(
        ax,
        heading,
    )


# ==========================================================================
# Load commonly used frozen surfaces
# ==========================================================================

claim_rows = read_tsv(
    claims
)

design_summary_rows = read_tsv(
    design_summary
)

design_summary_dict = {
    r["property"]: r["value"]
    for r in design_summary_rows
}

primary_rows = read_tsv(
    sources[
        "primary_global_results"
    ]
)

primary_by_id = {
    r["family_id"]: r
    for r in primary_rows
}

p1_component_rows = read_tsv(
    sources[
        "P1_direction_components"
    ]
)

p3_component_rows = read_tsv(
    sources[
        "P3_rank_centroid_components"
    ]
)

primary_perm_rows = read_tsv(
    sources[
        "primary_permutation_distributions"
    ]
)


# ==========================================================================
# FIGURE 1
# ==========================================================================

fig = new_figure(
    148
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.07,
    right=0.98,
    bottom=0.07,
    top=0.95,
    wspace=0.22,
    hspace=0.28,
)

axA = fig.add_subplot(
    gs[0, 0]
)
axB = fig.add_subplot(
    gs[0, 1]
)
axC = fig.add_subplot(
    gs[1, 0]
)
axD = fig.add_subplot(
    gs[1, 1]
)

# A — study design

axA.axis(
    "off"
)

panel_label(
    axA,
    "A",
)

panel_heading(
    axA,
    "Source study and paired design",
)

def box(
    ax,
    xy,
    wh,
    text,
    edge,
    face="white",
    fontsize=ANNOT_FONT,
):

    x, y = xy
    w, h = wh

    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.02",
        facecolor=face,
        edgecolor=edge,
        linewidth=1.0,
    )

    ax.add_patch(
        patch
    )

    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
    )

    return patch


def arrow(
    ax,
    p1,
    p2,
):

    a = FancyArrowPatch(
        p1,
        p2,
        arrowstyle="-|>",
        mutation_scale=8,
        linewidth=0.9,
        color=NEUTRAL,
    )

    ax.add_patch(
        a
    )


box(
    axA,
    (0.32, 0.78),
    (0.36, 0.12),
    "20 paired participants",
    DARK,
)

box(
    axA,
    (0.06, 0.47),
    (0.38, 0.14),
    "Reconstructed\nhypercholesterolaemic block\nn = 10",
    HYPER,
)

box(
    axA,
    (0.56, 0.47),
    (0.38, 0.14),
    "Reconstructed\ncontrol block\nn = 10",
    CONTROL,
)

arrow(
    axA,
    (0.43, 0.78),
    (0.25, 0.62),
)

arrow(
    axA,
    (0.57, 0.78),
    (0.75, 0.62),
)

box(
    axA,
    (0.08, 0.19),
    (0.15, 0.10),
    "A\nbaseline",
    HYPER,
)

box(
    axA,
    (0.27, 0.19),
    (0.15, 0.10),
    "B\nfollow-up",
    HYPER,
)

box(
    axA,
    (0.58, 0.19),
    (0.15, 0.10),
    "A\nbaseline",
    CONTROL,
)

box(
    axA,
    (0.77, 0.19),
    (0.15, 0.10),
    "B\nfollow-up",
    CONTROL,
)

arrow(
    axA,
    (0.155, 0.46),
    (0.155, 0.30),
)

arrow(
    axA,
    (0.345, 0.46),
    (0.345, 0.30),
)

arrow(
    axA,
    (0.655, 0.46),
    (0.655, 0.30),
)

arrow(
    axA,
    (0.845, 0.46),
    (0.845, 0.30),
)

axA.text(
    0.5,
    0.04,
    "Group mapping: provisional indirect high-confidence",
    ha="center",
    va="bottom",
    fontsize=ANNOT_FONT,
    color=NEUTRAL,
)

axA.set_xlim(
    0,
    1,
)

axA.set_ylim(
    0,
    1,
)


# B — five analyte panel

axB.axis(
    "off"
)

panel_label(
    axB,
    "B",
)

panel_heading(
    axB,
    "Five-channel oxysterol panel",
)

xpos = np.linspace(
    0.10,
    0.90,
    5,
)

for x, label in zip(
    xpos,
    ANALYTE_LABELS,
):

    c = Circle(
        (x, 0.55),
        0.075,
        facecolor="#F5F5F5",
        edgecolor=DARK,
        linewidth=1.0,
    )

    axB.add_patch(
        c
    )

    axB.text(
        x,
        0.55,
        label,
        ha="center",
        va="center",
        fontsize=ANNOT_FONT,
        fontweight="bold",
    )

axB.text(
    0.5,
    0.28,
    "Longitudinal analysis preserves direction, ranks,\n"
    "dependence topology and multivariate architecture",
    ha="center",
    va="center",
    fontsize=ANNOT_FONT,
)

axB.set_xlim(
    0,
    1,
)

axB.set_ylim(
    0,
    1,
)


# C — measurement scale

axC.axis(
    "off"
)

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "Measurement-scale contract",
)

box(
    axC,
    (0.05, 0.65),
    (0.25, 0.14),
    "Instrument\nCPS",
    DARK,
)

box(
    axC,
    (0.38, 0.65),
    (0.30, 0.14),
    "Internal-standard\nnormalized response",
    DARK,
)

box(
    axC,
    (0.76, 0.65),
    (0.19, 0.14),
    "Analysis\nsurface",
    DARK,
)

arrow(
    axC,
    (0.30, 0.72),
    (0.38, 0.72),
)

arrow(
    axC,
    (0.68, 0.72),
    (0.76, 0.72),
)

axC.text(
    0.5,
    0.42,
    "Participant absolute authentic-analyte concentrations\n"
    "are not established from the deposited workbook.",
    ha="center",
    va="center",
    fontsize=ANNOT_FONT,
)

axC.text(
    0.5,
    0.20,
    "Response-scale sensitivities use ln(response B/A),\n"
    "not concentration fold-change.",
    ha="center",
    va="center",
    fontsize=ANNOT_FONT,
    color=NEUTRAL,
)

axC.set_xlim(
    0,
    1,
)

axC.set_ylim(
    0,
    1,
)


# D — hierarchy

axD.axis(
    "off"
)

panel_label(
    axD,
    "D",
)

panel_heading(
    axD,
    "Frozen inferential hierarchy",
)

levels = [
    (
        0.75,
        "PRIMARY",
        "P1 directional topology\n"
        "P2 rank-network rewiring\n"
        "P3 rank-centroid displacement",
        HYPER,
    ),
    (
        0.53,
        "HIERARCHICAL",
        "P2 localization\n"
        "27-OHC ↔ 7β-OHC",
        SUPPORTED,
    ),
    (
        0.31,
        "SECONDARY",
        "S2 between-block displacement\n"
        "S1 control-anchored departure",
        "#6A51A3",
    ),
    (
        0.09,
        "SENSITIVITY",
        "S3 hyper 5/5 | control 0/5\n"
        "S4 between-block 5/5",
        CONTROL,
    ),
]

for y, label, text, edge in levels:

    box(
        axD,
        (0.08, y),
        (0.84, 0.15),
        (
            f"{label}\n"
            f"{text}"
        ),
        edge,
        fontsize=ANNOT_FONT,
    )

axD.set_xlim(
    0,
    1,
)

axD.set_ylim(
    0,
    1,
)

outputs = []

outputs.extend(
    save_figure(
        fig,
        required_stems[0],
        148,
    )
)


# ==========================================================================
# FIGURE 2
# ==========================================================================

fig = new_figure(
    174
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.08,
    right=0.98,
    bottom=0.08,
    top=0.95,
    wspace=0.28,
    hspace=0.34,
)

axA = fig.add_subplot(
    gs[0, 0]
)
axB = fig.add_subplot(
    gs[0, 1]
)
axC = fig.add_subplot(
    gs[1, 0]
)
axD = fig.add_subplot(
    gs[1, 1]
)

# A P1 subject-direction matrix

p1_rows = read_tsv(
    p1_subject
)

P1 = np.asarray([
    [
        int(
            r[
                f"{a}_direction"
            ]
        )
        for a in ANALYTES
    ]
    for r in p1_rows
])

disc_cmap = ListedColormap([
    NEGATIVE,
    "#F2F2F2",
    POSITIVE,
])

im = axA.imshow(
    P1,
    cmap=disc_cmap,
    vmin=-1,
    vmax=1,
    aspect="auto",
    interpolation="nearest",
)

axA.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

axA.set_yticks(
    range(10),
    [
        f"H{i:02d}"
        for i in range(
            1,
            11,
        )
    ],
)

axA.set_xlabel(
    "Oxysterol channel"
)

axA.set_ylabel(
    "Hyper block participant"
)

panel_label(
    axA,
    "A",
)

panel_heading(
    axA,
    "Subject-level longitudinal direction",
)

cax = inset_axes(
    axA,
    width="4%",
    height="55%",
    loc="center right",
    borderpad=1.2,
)

cb = fig.colorbar(
    im,
    cax=cax,
    ticks=[
        -2 / 3,
        0,
        2 / 3,
    ],
)

cb.ax.set_yticklabels([
    "Decrease",
    "No change",
    "Increase",
])

cb.ax.tick_params(
    labelsize=TICK_FONT,
)


# B P1 mean sign vector

p1_by = {
    r["analyte"]: r
    for r in p1_component_rows
}

vals = np.asarray([
    float(
        p1_by[a][
            "mean_sign"
        ]
    )
    for a in ANALYTES
])

y = np.arange(
    5
)

for yi, value in zip(
    y,
    vals,
):

    axB.plot(
        [
            0,
            value,
        ],
        [
            yi,
            yi,
        ],
        color=NEUTRAL,
        linewidth=1.2,
    )

    axB.scatter(
        value,
        yi,
        s=30,
        color=(
            NEGATIVE
            if value < 0
            else POSITIVE
        ),
        zorder=3,
    )

axB.axvline(
    0,
    color=DARK,
    linewidth=0.8,
)

axB.set_xlim(
    -1.05,
    1.05,
)

axB.set_yticks(
    y,
    ANALYTE_LABELS,
)

axB.invert_yaxis()

axB.set_xlabel(
    "Mean paired direction sign"
)

panel_label(
    axB,
    "B",
)

panel_heading(
    axB,
    "P1 directional components",
)

clean_axis(
    axB,
    grid=True,
)


# C P3 subject delta-z

p3_rows = read_tsv(
    p3_subject
)

P3 = np.asarray([
    [
        float(
            r[
                f"{a}_delta_z_B_minus_A"
            ]
        )
        for a in ANALYTES
    ]
    for r in p3_rows
])

v = max(
    1e-12,
    float(
        np.max(
            np.abs(P3)
        )
    ),
)

im3 = axC.imshow(
    P3,
    cmap="RdBu_r",
    vmin=-v,
    vmax=v,
    aspect="auto",
    interpolation="nearest",
)

axC.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

axC.set_yticks(
    range(10),
    [
        f"H{i:02d}"
        for i in range(
            1,
            11,
        )
    ],
)

axC.set_xlabel(
    "Oxysterol channel"
)

axC.set_ylabel(
    "Hyper block participant"
)

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "P3 pooled rank-normal displacement",
)

cax3 = inset_axes(
    axC,
    width="4%",
    height="55%",
    loc="center right",
    borderpad=1.2,
)

cb3 = fig.colorbar(
    im3,
    cax=cax3,
)

cb3.set_label(
    "Δz (B − A)",
    fontsize=ANNOT_FONT,
)


# D P1/P3 exact nulls

axD.axis(
    "off"
)

panel_label(
    axD,
    "D",
)

panel_heading(
    axD,
    "Frozen exact null distributions",
)

p1_col = find_family_stat_col(
    primary_perm_rows,
    "P1",
)

p3_col = find_family_stat_col(
    primary_perm_rows,
    "P3",
)

p1_perm = np.asarray([
    float(
        r[p1_col]
    )
    for r in primary_perm_rows
])

p3_perm = np.asarray([
    float(
        r[p3_col]
    )
    for r in primary_perm_rows
])

axD1 = inset_axes(
    axD,
    width="94%",
    height="39%",
    loc="upper center",
    borderpad=1.25,
)

axD2 = inset_axes(
    axD,
    width="94%",
    height="39%",
    loc="lower center",
    borderpad=1.25,
)

p1_obs = float(
    primary_by_id[
        "P1"
    ][
        "observed_statistic"
    ]
)

p3_obs = float(
    primary_by_id[
        "P3"
    ][
        "observed_statistic"
    ]
)

p1_q = float(
    primary_by_id[
        "P1"
    ][
        "holm_p_adjusted"
    ]
)

p3_q = float(
    primary_by_id[
        "P3"
    ][
        "holm_p_adjusted"
    ]
)

histogram_panel(
    axD1,
    p1_perm,
    p1_obs,
    (
        "P1 Holm p="
        + fmt_p(
            p1_q
        )
    ),
    "P1 statistic",
    bins=28,
)

histogram_panel(
    axD2,
    p3_perm,
    p3_obs,
    (
        "P3 Holm p="
        + fmt_p(
            p3_q
        )
    ),
    "P3 statistic",
    bins=28,
)

outputs.extend(
    save_figure(
        fig,
        required_stems[1],
        174,
    )
)


# ==========================================================================
# FIGURE 3
# ==========================================================================

fig = new_figure(
    176
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.07,
    right=0.98,
    bottom=0.08,
    top=0.95,
    wspace=0.25,
    hspace=0.30,
)

axA = fig.add_subplot(
    gs[0, 0]
)
axB = fig.add_subplot(
    gs[0, 1]
)
axC = fig.add_subplot(
    gs[1, 0]
)
axD = fig.add_subplot(
    gs[1, 1]
)

labels_A, cols_A, net_A = (
    matrix_from_tsv(
        sources[
            "P2_rank_network_A"
        ]
    )
)

labels_B, cols_B, net_B = (
    matrix_from_tsv(
        sources[
            "P2_rank_network_B"
        ]
    )
)

if labels_A != labels_B:

    raise RuntimeError(
        "P2 network node ordering differs"
    )

draw_network(
    axA,
    net_A,
    labels_A,
    "Baseline rank network (A)",
)

panel_label(
    axA,
    "A",
)

draw_network(
    axB,
    net_B,
    labels_B,
    "Follow-up rank network (B)",
)

panel_label(
    axB,
    "B",
)

legend_handles = [
    Line2D(
        [0],
        [0],
        color=POSITIVE,
        linewidth=2,
        label="Positive ρ",
    ),
    Line2D(
        [0],
        [0],
        color=NEGATIVE,
        linewidth=2,
        label="Negative ρ",
    ),
]

axB.legend(
    handles=legend_handles,
    loc="lower center",
    bbox_to_anchor=(
        0.5,
        -0.08,
    ),
    frameon=False,
    fontsize=ANNOT_FONT,
    ncol=2,
)


# C delta-rho matrix from frozen edge table

delta_rows = read_tsv(
    sources[
        "P2_edge_deltas"
    ]
)

delta_matrix = np.zeros(
    (5, 5),
    dtype=float,
)

index = {
    a: i
    for i, a in enumerate(
        ANALYTES
    )
}

for row in delta_rows:

    a1 = first_existing(
        row,
        [
            "analyte_1",
            "analyte1",
        ],
    )

    a2 = first_existing(
        row,
        [
            "analyte_2",
            "analyte2",
        ],
    )

    delta = float(
        first_existing(
            row,
            [
                "delta_rho_B_minus_A",
                "delta_rho",
            ],
        )
    )

    i = index[
        a1
    ]
    j = index[
        a2
    ]

    delta_matrix[
        i,
        j
    ] = delta

    delta_matrix[
        j,
        i
    ] = delta

vdelta = max(
    1e-12,
    float(
        np.max(
            np.abs(
                delta_matrix
            )
        )
    ),
)

imd = axC.imshow(
    delta_matrix,
    cmap="RdBu_r",
    vmin=-vdelta,
    vmax=vdelta,
    interpolation="nearest",
)

axC.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

axC.set_yticks(
    range(5),
    ANALYTE_LABELS,
)

for i in range(5):

    for j in range(5):

        if i == j:
            continue

        axC.text(
            j,
            i,
            f"{delta_matrix[i,j]:.2f}",
            ha="center",
            va="center",
            fontsize=6.5,
            color=DARK,
        )

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "Change in rank dependence",
)

cax = inset_axes(
    axC,
    width="4%",
    height="60%",
    loc="center right",
    borderpad=1.2,
)

cb = fig.colorbar(
    imd,
    cax=cax,
)

cb.set_label(
    "Δρ = ρB − ρA",
    fontsize=ANNOT_FONT,
)


# D hierarchical maxT edges

maxt_rows = read_tsv(
    sources[
        "P2_exact_maxT_results"
    ]
)

def edge_delta(
    row,
):

    return float(
        first_existing(
            row,
            [
                "delta_rho_B_minus_A",
                "delta_rho",
            ],
        )
    )


ordered = sorted(
    maxt_rows,
    key=lambda r: abs(
        edge_delta(r)
    ),
    reverse=True,
)

yy = np.arange(
    len(ordered)
)

edge_labels = [
    (
        pretty_analyte(
            r["analyte_1"]
        )
        + "–"
        + pretty_analyte(
            r["analyte_2"]
        )
    )
    for r in ordered
]

edge_values = np.asarray([
    edge_delta(r)
    for r in ordered
])

support = [
    r[
        "maxT_support_0_05"
    ] == "YES"
    for r in ordered
]

for yi, value, supported in zip(
    yy,
    edge_values,
    support,
):

    axD.plot(
        [
            0,
            value,
        ],
        [
            yi,
            yi,
        ],
        color=LIGHT,
        linewidth=1.0,
    )

    axD.scatter(
        value,
        yi,
        s=(
            38
            if supported
            else 24
        ),
        color=(
            SUPPORTED
            if supported
            else NEUTRAL
        ),
        edgecolor=(
            DARK
            if supported
            else "none"
        ),
        linewidth=0.6,
        zorder=3,
    )

axD.axvline(
    0,
    color=DARK,
    linewidth=0.8,
)

axD.set_yticks(
    yy,
    edge_labels,
)

axD.invert_yaxis()

axD.set_xlabel(
    "Δρ (B − A)"
)

panel_label(
    axD,
    "D",
)

panel_heading(
    axD,
    "Hierarchical edge localization",
)

clean_axis(
    axD,
    grid=True,
)

xmin = min(
    -0.05,
    float(
        np.min(
            edge_values
        )
    ) - 0.10,
)

xmax = max(
    0.05,
    float(
        np.max(
            edge_values
        )
    ) + 0.32,
)

axD.set_xlim(
    xmin,
    xmax,
)

for yi, row, value in zip(
    yy,
    ordered,
    edge_values,
):

    p = float(
        row[
            "maxT_adjusted_p"
        ]
    )

    axD.text(
        xmax - 0.01,
        yi,
        (
            "maxT p="
            + fmt_p(p)
        ),
        ha="right",
        va="center",
        fontsize=6.2,
        color=(
            DARK
            if row[
                "maxT_support_0_05"
            ] == "YES"
            else NEUTRAL
        ),
    )

outputs.extend(
    save_figure(
        fig,
        required_stems[2],
        176,
    )
)


# ==========================================================================
# FIGURE 4
# ==========================================================================

fig = new_figure(
    154
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

s2_component_rows = read_tsv(
    sources[
        "S2_group_difference_components"
    ]
)

s2_by = {
    r["analyte"]: r
    for r in s2_component_rows
}

s2_values = np.asarray([
    float(
        s2_by[a][
            "hyper_minus_control"
        ]
    )
    for a in ANALYTES
])

y = np.arange(
    5
)

for yi, value in zip(
    y,
    s2_values,
):

    axA.plot(
        [
            0,
            value,
        ],
        [
            yi,
            yi,
        ],
        color=LIGHT,
        linewidth=1.3,
    )

    axA.scatter(
        value,
        yi,
        color=(
            NEGATIVE
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
    ANALYTE_LABELS,
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


# B subject heatmap

s2_subject_rows = read_tsv(
    sources[
        "S2_subject_rank_displacements"
    ]
)

S2mat = np.asarray([
    [
        float(
            row[
                f"{a}_delta_rank"
            ]
        )
        for a in ANALYTES
    ]
    for row in s2_subject_rows
])

v2 = max(
    1e-12,
    float(
        np.max(
            np.abs(
                S2mat
            )
        )
    ),
)

im2 = axB.imshow(
    S2mat,
    cmap="RdBu_r",
    vmin=-v2,
    vmax=v2,
    aspect="auto",
    interpolation="nearest",
)

axB.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

subject_labels = []

for row in s2_subject_rows:

    prefix = (
        "H"
        if row[
            "clinical_group"
        ] == "HYPERCHOLESTEROLAEMIC"
        else "C"
    )

    subject_labels.append(
        (
            prefix
            + str(
                row[
                    "subject_id"
                ]
            )
        )
    )

axB.set_yticks(
    range(
        len(
            s2_subject_rows
        )
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
    im2,
    cax=cax,
)

cb.set_label(
    "Δrank (B − A)",
    fontsize=ANNOT_FONT,
)


# C exact null

s2_global_rows = read_tsv(
    sources[
        "S2_global_result"
    ]
)

s2_global = s2_global_rows[
    0
]

s2_perm_rows = read_tsv(
    sources[
        "S2_permutation_distribution"
    ]
)

s2_perm_values = np.asarray([
    float(
        r[
            "global_statistic"
        ]
    )
    for r in s2_perm_rows
])

s2_obs = float(
    s2_global[
        "observed_statistic"
    ]
)

s2_p = float(
    s2_global[
        "exact_p"
    ]
)

histogram_panel(
    axC,
    s2_perm_values,
    s2_obs,
    (
        "Exact p="
        + fmt_p(
            s2_p
        )
    ),
    "S2 global rank-displacement statistic",
    bins=55,
)

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "Complete 10-of-20 label-permutation reference",
)

outputs.extend(
    save_figure(
        fig,
        required_stems[3],
        154,
    )
)


# ==========================================================================
# FIGURE 5
# ==========================================================================

fig = new_figure(
    176
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.08,
    right=0.98,
    bottom=0.09,
    top=0.95,
    wspace=0.30,
    hspace=0.36,
)

axA = fig.add_subplot(
    gs[0, 0]
)
axB = fig.add_subplot(
    gs[0, 1]
)
axC = fig.add_subplot(
    gs[1, 0]
)
axD = fig.add_subplot(
    gs[1, 1]
)

s1_subject_rows = read_tsv(
    sources[
        "S1_subject_departures"
    ]
)

# A paired D2

for row in s1_subject_rows:

    a = float(
        row[
            "D2_A"
        ]
    )

    b = float(
        row[
            "D2_B"
        ]
    )

    axA.plot(
        [
            0,
            1,
        ],
        [
            a,
            b,
        ],
        color=NEUTRAL,
        alpha=0.65,
        marker="o",
        markersize=3.5,
    )

axA.set_xticks(
    [
        0,
        1,
    ],
    [
        "A",
        "B",
    ],
)

axA.set_ylabel(
    "Squared Mahalanobis departure, D²"
)

panel_label(
    axA,
    "A",
)

panel_heading(
    axA,
    "Paired departure trajectories",
)

clean_axis(
    axA,
    grid=True,
)


# B delta D2 sorted

s1_sorted = sorted(
    s1_subject_rows,
    key=lambda r: float(
        r[
            "delta_D2_B_minus_A"
        ]
    ),
)

deltas = np.asarray([
    float(
        r[
            "delta_D2_B_minus_A"
        ]
    )
    for r in s1_sorted
])

yy = np.arange(
    len(
        s1_sorted
    )
)

for yi, value in zip(
    yy,
    deltas,
):

    axB.plot(
        [
            0,
            value,
        ],
        [
            yi,
            yi,
        ],
        color=LIGHT,
        linewidth=1.1,
    )

    axB.scatter(
        value,
        yi,
        s=28,
        color=(
            NEGATIVE
            if value < 0
            else POSITIVE
        ),
        zorder=3,
    )

axB.axvline(
    0,
    color=DARK,
    linewidth=0.8,
)

axB.set_yticks(
    yy,
    [
        "H"
        + str(
            r[
                "subject_id"
            ]
        )
        for r in s1_sorted
    ],
)

axB.set_xlabel(
    "ΔD² = D²B − D²A"
)

panel_label(
    axB,
    "B",
)

panel_heading(
    axB,
    "Subject movement relative to reference",
)

clean_axis(
    axB,
    grid=True,
)


# C control covariance

ref_rows = read_tsv(
    sources[
        "S1_control_reference"
    ]
)

cov_rows = [
    r
    for r in ref_rows
    if r[
        "object"
    ] == "covariance"
]

diag_rows = [
    r
    for r in ref_rows
    if r[
        "object"
    ] == "diagnostic"
]

cov = np.zeros(
    (5, 5),
    dtype=float,
)

for row in cov_rows:

    i = index[
        row[
            "row"
        ]
    ]

    j = index[
        row[
            "column"
        ]
    ]

    cov[
        i,
        j
    ] = float(
        row[
            "value"
        ]
    )

diag = {
    r["row"]: r["value"]
    for r in diag_rows
}

vcov = max(
    1e-12,
    float(
        np.max(
            np.abs(
                cov
            )
        )
    ),
)

imc = axC.imshow(
    cov,
    cmap="RdBu_r",
    vmin=-vcov,
    vmax=vcov,
    interpolation="nearest",
)

axC.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

axC.set_yticks(
    range(5),
    ANALYTE_LABELS,
)

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "OAS-shrunk control-baseline covariance",
)

shrink = float(
    diag[
        "oas_shrinkage"
    ]
)

cond = float(
    diag[
        "covariance_condition_number"
    ]
)

axC.text(
    0.02,
    -0.23,
    (
        f"OAS shrinkage={shrink:.3f}; "
        f"condition number={cond:.2f}"
    ),
    transform=axC.transAxes,
    ha="left",
    va="top",
    fontsize=ANNOT_FONT,
)

cax = inset_axes(
    axC,
    width="4%",
    height="58%",
    loc="center right",
    borderpad=1.1,
)

fig.colorbar(
    imc,
    cax=cax,
)


# D exact null

s1_global = read_tsv(
    sources[
        "S1_global_result"
    ]
)[0]

s1_perm_rows = read_tsv(
    sources[
        "S1_permutation_distribution"
    ]
)

s1_perm_values = np.asarray([
    float(
        r[
            "abs_mean_delta_D2"
        ]
    )
    for r in s1_perm_rows
])

s1_obs = float(
    s1_global[
        "observed_abs_mean_delta_D2"
    ]
)

s1_signed = float(
    s1_global[
        "signed_mean_delta_D2"
    ]
)

s1_p = float(
    s1_global[
        "exact_p"
    ]
)

histogram_panel(
    axD,
    s1_perm_values,
    s1_obs,
    (
        "Exact p="
        + fmt_p(
            s1_p
        )
        + "\nSigned mean ΔD²="
        + f"{s1_signed:.2f}"
    ),
    "|mean ΔD²|",
    bins=36,
)

panel_label(
    axD,
    "D",
)

panel_heading(
    axD,
    "Exact whole-subject A/B swap reference",
)

outputs.extend(
    save_figure(
        fig,
        required_stems[4],
        176,
    )
)


# ==========================================================================
# FIGURE 6
# ==========================================================================

fig = new_figure(
    174
)

gs = fig.add_gridspec(
    2,
    2,
    left=0.09,
    right=0.98,
    bottom=0.09,
    top=0.95,
    wspace=0.33,
    hspace=0.36,
)

axA = fig.add_subplot(
    gs[0, 0]
)
axB = fig.add_subplot(
    gs[0, 1]
)
axC = fig.add_subplot(
    gs[1, 0]
)
axD = fig.add_subplot(
    gs[1, 1]
)

s3h_rows = read_tsv(
    sources[
        "S3_hyper_results"
    ]
)

s3c_rows = read_tsv(
    sources[
        "S3_control_results"
    ]
)

s4_rows = read_tsv(
    sources[
        "S4_results"
    ]
)

s3h_by = {
    r["analyte"]: r
    for r in s3h_rows
}

s3c_by = {
    r["analyte"]: r
    for r in s3c_rows
}

s4_by = {
    r["analyte"]: r
    for r in s4_rows
}


def robust_response_panel(
    ax,
    rows_by,
    color,
    heading,
    panel,
):

    y = np.arange(
        5
    )

    means = np.asarray([
        float(
            rows_by[a][
                "mean_log_response"
            ]
        )
        for a in ANALYTES
    ])

    medians = np.asarray([
        float(
            rows_by[a][
                "median_log_response"
            ]
        )
        for a in ANALYTES
    ])

    q1 = np.asarray([
        float(
            rows_by[a][
                "q1_log_response"
            ]
        )
        for a in ANALYTES
    ])

    q3 = np.asarray([
        float(
            rows_by[a][
                "q3_log_response"
            ]
        )
        for a in ANALYTES
    ])

    for yi, m, med, lo, hi, a in zip(
        y,
        means,
        medians,
        q1,
        q3,
        ANALYTES,
    ):

        ax.plot(
            [
                lo,
                hi,
            ],
            [
                yi,
                yi,
            ],
            color=color,
            linewidth=2.0,
            alpha=0.75,
        )

        ax.scatter(
            med,
            yi,
            s=32,
            facecolor="white",
            edgecolor=color,
            linewidth=1.2,
            zorder=3,
        )

        ax.scatter(
            m,
            yi,
            s=24,
            marker="x",
            color=DARK,
            linewidth=1.0,
            zorder=4,
        )

        p = float(
            rows_by[a][
                "holm_p_adjusted"
            ]
        )

        ax.text(
            0.98,
            yi,
            (
                "Holm p="
                + fmt_p(p)
            ),
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=6.2,
        )

    ax.axvline(
        0,
        color=DARK,
        linewidth=0.8,
    )

    ax.set_yticks(
        y,
        ANALYTE_LABELS,
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "ln analytical response B/A"
    )

    panel_label(
        ax,
        panel,
    )

    panel_heading(
        ax,
        heading,
    )

    clean_axis(
        ax,
        grid=True,
    )


robust_response_panel(
    axA,
    s3h_by,
    HYPER,
    "S3 hyper-block analytical response",
    "A",
)

robust_response_panel(
    axB,
    s3c_by,
    CONTROL,
    "S3 control-block analytical response",
    "B",
)


# C S4

y = np.arange(
    5
)

s4_values = np.asarray([
    float(
        s4_by[a][
            "hyper_minus_control_mean_log_response"
        ]
    )
    for a in ANALYTES
])

for yi, value, a in zip(
    y,
    s4_values,
    ANALYTES,
):

    axC.plot(
        [
            0,
            value,
        ],
        [
            yi,
            yi,
        ],
        color=LIGHT,
        linewidth=1.2,
    )

    axC.scatter(
        value,
        yi,
        s=32,
        color=HYPER,
        zorder=3,
    )

    p = float(
        s4_by[a][
            "holm_p_adjusted"
        ]
    )

    axC.text(
        0.98,
        yi,
        (
            "Holm p="
            + fmt_p(p)
        ),
        transform=axC.get_yaxis_transform(),
        ha="right",
        va="center",
        fontsize=6.2,
    )

axC.axvline(
    0,
    color=DARK,
    linewidth=0.8,
)

axC.set_yticks(
    y,
    ANALYTE_LABELS,
)

axC.invert_yaxis()

axC.set_xlabel(
    "Hyper − control mean ln response"
)

panel_label(
    axC,
    "C",
)

panel_heading(
    axC,
    "S4 between-block contrast",
)

clean_axis(
    axC,
    grid=True,
)


# D support matrix

support_matrix = np.zeros(
    (3, 5),
    dtype=int,
)

for j, a in enumerate(
    ANALYTES
):

    support_matrix[
        0,
        j
    ] = int(
        s3h_by[a][
            "holm_support_0_05"
        ] == "YES"
    )

    support_matrix[
        1,
        j
    ] = int(
        s3c_by[a][
            "holm_support_0_05"
        ] == "YES"
    )

    support_matrix[
        2,
        j
    ] = int(
        s4_by[a][
            "holm_support_0_05"
        ] == "YES"
    )

support_cmap = ListedColormap([
    "#E6E6E6",
    SUPPORTED,
])

ims = axD.imshow(
    support_matrix,
    cmap=support_cmap,
    vmin=0,
    vmax=1,
    aspect="auto",
    interpolation="nearest",
)

axD.set_xticks(
    range(5),
    ANALYTE_LABELS,
    rotation=35,
    ha="right",
)

axD.set_yticks(
    range(3),
    [
        "S3 hyper",
        "S3 control",
        "S4 contrast",
    ],
)

for i in range(3):

    for j in range(5):

        axD.text(
            j,
            i,
            (
                "Supported"
                if support_matrix[
                    i,
                    j
                ]
                else "—"
            ),
            ha="center",
            va="center",
            fontsize=6.1,
            color=(
                "white"
                if support_matrix[
                    i,
                    j
                ]
                else NEUTRAL
            ),
            fontweight=(
                "bold"
                if support_matrix[
                    i,
                    j
                ]
                else "normal"
            ),
        )

panel_label(
    axD,
    "D",
)

panel_heading(
    axD,
    "Multiplicity-controlled support matrix",
)

outputs.extend(
    save_figure(
        fig,
        required_stems[5],
        174,
    )
)


# ==========================================================================
# 7. Output count / format validation
# ==========================================================================

if len(
    outputs
) != 24:

    print(
        f"REFUSE_OUTPUT_RECORD_COUNT="
        f"{len(outputs)}"
    )

    raise SystemExit(11)

format_counts = defaultdict(
    int
)

for fmt, *_ in outputs:

    format_counts[
        fmt
    ] += 1

required_counts = {
    "TIFF": 6,
    "SVG": 6,
    "PDF": 6,
    "PNG_QA": 6,
}

if dict(
    format_counts
) != required_counts:

    print(
        f"REFUSE_FORMAT_COUNTS="
        f"{dict(format_counts)}"
    )

    raise SystemExit(12)

print()
print("PUBLICATION_FIGURE_COUNT=6")
print("TIFF_COUNT=6")
print("SVG_COUNT=6")
print("PDF_COUNT=6")
print("QA_PNG_COUNT=6")
print("OUTPUT_COUNT_QC=PASS")
print("RASTER_DIMENSION_QC=PASS")


# ==========================================================================
# 8. Output manifest
# ==========================================================================

manifest_path = (
    out
    / "OX_G2A_R2Z_publication_figure_output_manifest.tsv"
)

with manifest_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([
        "format",
        "path",
        "bytes",
        "sha256",
        "pixel_width",
        "pixel_height",
        "dpi",
        "compression",
    ])

    for (
        fmt,
        path,
        width,
        height,
        dpi,
        compression,
    ) in outputs:

        w.writerow([
            fmt,
            str(path),
            path.stat().st_size,
            sha256(path),
            width,
            height,
            dpi,
            compression,
        ])


# ==========================================================================
# 9. Figure-level source provenance
# ==========================================================================

provenance_path = (
    out
    / "OX_G2A_R2Z_figure_provenance.tsv"
)

figure_sources = {
    "Figure_1": [
        "study_design_summary",
    ],
    "Figure_2": [
        "P1_subject_direction_matrix",
        "P1_direction_components",
        "P3_subject_rank_score_displacement",
        "P3_rank_centroid_components",
        "primary_global_results",
        "primary_permutation_distributions",
    ],
    "Figure_3": [
        "P2_rank_network_A",
        "P2_rank_network_B",
        "P2_edge_deltas",
        "P2_exact_maxT_results",
    ],
    "Figure_4": [
        "S2_global_result",
        "S2_group_difference_components",
        "S2_subject_rank_displacements",
        "S2_permutation_distribution",
    ],
    "Figure_5": [
        "S1_global_result",
        "S1_subject_departures",
        "S1_control_reference",
        "S1_permutation_distribution",
    ],
    "Figure_6": [
        "S3_hyper_results",
        "S3_control_results",
        "S4_results",
    ],
}

manifest_by_role = {
    r["source_role"]: r
    for r in manifest_rows
}

with provenance_path.open(
    "w",
    encoding="utf-8",
    newline="",
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([
        "figure",
        "source_role",
        "source_path",
        "source_sha256",
        "source_class",
    ])

    for figure, roles in figure_sources.items():

        for role in roles:

            row = manifest_by_role[
                role
            ]

            w.writerow([
                figure,
                role,
                row["path"],
                row["sha256"],
                row["source_class"],
            ])


# ==========================================================================
# 10. Visual-QA checklist
# ==========================================================================

qa_checklist = (
    out
    / "OX_G2A_R2Z_manual_visual_QA_checklist.tsv"
)

qa_rows = []

for stem in required_stems:

    qa_rows.extend([
        [
            stem,
            "NO_CLIPPED_TEXT",
            "PENDING_MANUAL_REVIEW",
        ],
        [
            stem,
            "NO_PANEL_OVERLAP",
            "PENDING_MANUAL_REVIEW",
        ],
        [
            stem,
            "LEGENDS_AND_LABELS_READABLE",
            "PENDING_MANUAL_REVIEW",
        ],
        [
            stem,
            "ANALYTE_NAMES_CONSISTENT",
            "PENDING_MANUAL_REVIEW",
        ],
        [
            stem,
            "PANEL_LABELS_CORRECT",
            "PENDING_MANUAL_REVIEW",
        ],
        [
            stem,
            "STATISTICAL_ANNOTATIONS_MATCH_FROZEN_RESULTS",
            "PENDING_MANUAL_REVIEW",
        ],
    ])

with qa_checklist.open(
    "w",
    encoding="utf-8",
    newline="",
) as f:

    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([
        "figure",
        "visual_check",
        "status",
    ])

    w.writerows(
        qa_rows
    )


# ==========================================================================
# 11. Governance receipt
# ==========================================================================

receipt = (
    gov
    / "OX_G2A_R2Z_publication_figure_generation_receipt_20260928.txt"
)

receipt_lines = [
    "OX_COMORBID_G2A_R2Z_PUBLICATION_FIGURE_GENERATION",
    "DATE=2026-09-28",
    f"R2Y_RECEIPT_SHA256={sha256(r2y_receipt)}",
    f"R2Y_FREEZE_SHA256={sha256(r2y_freeze)}",
    f"R2Y_SOURCE_MANIFEST_SHA256={sha256(source_manifest)}",
    f"R2Y_DESIGN_SYSTEM_SHA256={sha256(design_system)}",
    f"R2Y_PANEL_SPEC_SHA256={sha256(panel_spec)}",
    f"R2Y_NO_RECALC_SHA256={sha256(no_recalc)}",
    f"OUTPUT_MANIFEST_SHA256={sha256(manifest_path)}",
    f"FIGURE_PROVENANCE_SHA256={sha256(provenance_path)}",
    f"MANUAL_QA_CHECKLIST_SHA256={sha256(qa_checklist)}",
    "FIGURE_SOURCE_REVALIDATION=PASS",
    "EXPECTED_OUTPUT_IDENTITY=PASS",
    "PUBLICATION_FIGURE_COUNT=6",
    "TIFF_COUNT=6",
    "SVG_COUNT=6",
    "PDF_COUNT=6",
    "QA_PNG_COUNT=6",
    "PUBLICATION_TIFF_DPI=600",
    "QA_PNG_DPI=200",
    "RASTER_DIMENSION_QC=PASS",
    "FIXED_CANVAS_POLICY=PASS",
    "NEW_P_VALUES_COMPUTED=NO",
    "NEW_HYPOTHESIS_TESTS_EXECUTED=NO",
    "NEW_INFERENTIAL_FAMILIES=HOLD",
    "UPSTREAM_STATISTICAL_VALUES_REUSED=YES",
    "ABSOLUTE_CONCENTRATION_INTERPRETATION=NOT_AUTHORIZED",
    "CONCENTRATION_FOLD_CHANGE_INTERPRETATION=NOT_AUTHORIZED",
    "HOMEOSTASIS_RESTORATION_CLAIM=NOT_AUTHORIZED",
    "CAUSAL_TREATMENT_EFFECT_INTERPRETATION=NOT_AUTHORIZED",
    "VISUAL_QA_ADJUDICATION=PENDING",
    "PUBLICATION_FIGURES_FINAL_FREEZE=NOT_YET_AUTHORIZED",
]

for stem in required_stems:

    for ext in [
        ".tiff",
        ".svg",
        ".pdf",
        "_QA.png",
    ]:

        path = (
            figdir
            / (
                stem
                + ext
            )
        )

        key = (
            stem.upper()
            .replace(
                "-",
                "_",
            )
            .replace(
                " ",
                "_",
            )
        )

        ext_key = (
            ext
            .replace(
                ".",
                ""
            )
            .replace(
                "_",
                ""
            )
            .upper()
        )

        receipt_lines.append(
            f"{key}_{ext_key}_SHA256="
            f"{sha256(path)}"
        )

receipt.write_text(
    "\n".join(
        receipt_lines
    ) + "\n",
    encoding="utf-8",
)

freeze = (
    gov
    / "OX_G2A_R2Z_publication_figure_generation_receipt_20260928.sha256"
)

freeze_targets = [
    manifest_path,
    provenance_path,
    qa_checklist,
    receipt,
]

for stem in required_stems:

    freeze_targets.extend([
        figdir
        / (
            stem
            + ".tiff"
        ),
        figdir
        / (
            stem
            + ".svg"
        ),
        figdir
        / (
            stem
            + ".pdf"
        ),
        figdir
        / (
            stem
            + "_QA.png"
        ),
    ])

with freeze.open(
    "w",
    encoding="ascii",
    newline="\n",
) as f:

    for path in freeze_targets:

        f.write(
            f"{sha256(path)}  {path}\n"
        )


# ==========================================================================
# 12. Final terminal state
# ==========================================================================

print()
print("----- R2Z OUTPUT MANIFEST -----")

for row in read_tsv(
    manifest_path
):

    print(
        f"FORMAT={row['format']} "
        f"FILE={Path(row['path']).name} "
        f"BYTES={row['bytes']} "
        f"SHA256={row['sha256']}"
    )

print()
print("----- R2Z GOVERNANCE -----")
print(
    f"RECEIPT={receipt}"
)
print(
    f"RECEIPT_SHA256="
    f"{sha256(receipt)}"
)
print(
    f"FREEZE={freeze}"
)
print(
    f"FREEZE_SHA256="
    f"{sha256(freeze)}"
)

print()
print("=" * 112)
print("R2Z_EXECUTION=PASS")
print("PUBLICATION_FIGURES_GENERATED=YES")
print("PUBLICATION_FIGURE_COUNT=6")
print("TIFF_COUNT=6")
print("SVG_COUNT=6")
print("PDF_COUNT=6")
print("QA_PNG_COUNT=6")
print("FIGURE_SOURCE_REVALIDATION=PASS")
print("RASTER_DIMENSION_QC=PASS")
print("NEW_P_VALUES_COMPUTED=NO")
print("NEW_HYPOTHESIS_TESTS_EXECUTED=NO")
print("NEW_INFERENTIAL_FAMILIES=HOLD")
print("VISUAL_QA_ADJUDICATION=PENDING")
print("PUBLICATION_FIGURES_FINAL_FREEZE=NOT_YET_AUTHORIZED")
print("NEXT_GATE=R2Z_MANUAL_VISUAL_QA")
print("=" * 112)

