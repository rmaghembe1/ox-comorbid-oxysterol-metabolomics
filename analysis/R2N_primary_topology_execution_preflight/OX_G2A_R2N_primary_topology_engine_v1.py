
from statistics import NormalDist
import math
import numpy as np

ANALYTE_ORDER = (
    "24S_OHC",
    "25_OHC",
    "27_OHC",
    "7B_OHC",
    "7_KC",
)

FLOAT_TAIL_REL_TOL = 1e-12


def _as_matrix(x):
    x = np.asarray(x, dtype=float)

    if x.ndim != 2:
        raise ValueError("matrix must be two-dimensional")

    if not np.all(np.isfinite(x)):
        raise ValueError("matrix contains non-finite values")

    return x


def average_ranks_1d(x):
    x = np.asarray(x, dtype=float)

    if x.ndim != 1:
        raise ValueError("rank input must be one-dimensional")

    if not np.all(np.isfinite(x)):
        raise ValueError("rank input contains non-finite values")

    n = x.size
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(n, dtype=float)

    i = 0

    while i < n:
        j = i + 1

        while (
            j < n
            and x[order[j]] == x[order[i]]
        ):
            j += 1

        # Average 1-based ranks from i+1 through j.
        rank = ((i + 1) + j) / 2.0

        ranks[order[i:j]] = rank

        i = j

    return ranks


def rank_columns(x):
    x = _as_matrix(x)

    return np.column_stack([
        average_ranks_1d(x[:, j])
        for j in range(x.shape[1])
    ])


def spearman_matrix(x):
    x = _as_matrix(x)

    if x.shape[0] < 3:
        raise ValueError("at least three observations required")

    r = rank_columns(x)

    sds = r.std(axis=0, ddof=0)

    if np.any(sds == 0):
        raise ValueError("constant ranked variable")

    corr = np.corrcoef(
        r,
        rowvar=False
    )

    return corr


def p1_directional_topology(a, b):
    a = _as_matrix(a)
    b = _as_matrix(b)

    if a.shape != b.shape:
        raise ValueError("A and B shapes differ")

    signs = np.sign(b - a)

    mean_sign = signs.mean(axis=0)

    statistic = float(
        np.linalg.norm(
            mean_sign,
            ord=2
        )
    )

    return statistic, mean_sign, signs


def p2_rank_network_rewiring(a, b):
    a = _as_matrix(a)
    b = _as_matrix(b)

    if a.shape != b.shape:
        raise ValueError("A and B shapes differ")

    ra = spearman_matrix(a)
    rb = spearman_matrix(b)

    diff = rb - ra

    iu = np.triu_indices(
        a.shape[1],
        k=1
    )

    edge_delta = diff[iu]

    statistic = float(
        np.sqrt(
            np.sum(
                edge_delta ** 2
            )
        )
    )

    return statistic, ra, rb, edge_delta


def pooled_rank_normal_scores(a, b):
    a = _as_matrix(a)
    b = _as_matrix(b)

    if a.shape != b.shape:
        raise ValueError("A and B shapes differ")

    pooled = np.vstack([
        a,
        b,
    ])

    n_total = pooled.shape[0]

    ranks = rank_columns(
        pooled
    )

    probs = (
        ranks - 0.5
    ) / float(n_total)

    normal = NormalDist()

    z = np.empty_like(
        probs,
        dtype=float
    )

    for i in range(
        probs.shape[0]
    ):
        for j in range(
            probs.shape[1]
        ):
            p = float(
                probs[i, j]
            )

            if not (
                0.0 < p < 1.0
            ):
                raise ValueError(
                    "rank-normal probability outside (0,1)"
                )

            z[i, j] = (
                normal.inv_cdf(p)
            )

    n = a.shape[0]

    return (
        z[:n].copy(),
        z[n:].copy(),
        ranks[:n].copy(),
        ranks[n:].copy(),
    )


def p3_rank_centroid_displacement_from_scores(
    za,
    zb,
):
    za = _as_matrix(za)
    zb = _as_matrix(zb)

    if za.shape != zb.shape:
        raise ValueError("ZA and ZB shapes differ")

    delta = zb - za

    mean_delta = delta.mean(
        axis=0
    )

    statistic = float(
        np.linalg.norm(
            mean_delta,
            ord=2
        )
    )

    return (
        statistic,
        mean_delta,
        delta,
    )


def p3_rank_centroid_displacement(a, b):
    za, zb, ra, rb = (
        pooled_rank_normal_scores(
            a,
            b,
        )
    )

    statistic, mean_delta, delta = (
        p3_rank_centroid_displacement_from_scores(
            za,
            zb,
        )
    )

    return (
        statistic,
        mean_delta,
        delta,
        za,
        zb,
        ra,
        rb,
    )


def apply_time_swap(a, b, mask):
    a = _as_matrix(a)
    b = _as_matrix(b)

    if a.shape != b.shape:
        raise ValueError("A and B shapes differ")

    n = a.shape[0]

    if mask < 0 or mask >= (1 << n):
        raise ValueError("mask outside enumeration space")

    ap = a.copy()
    bp = b.copy()

    for i in range(n):

        if (
            mask >> i
        ) & 1:

            tmp = ap[i].copy()
            ap[i] = bp[i]
            bp[i] = tmp

    return ap, bp


def iter_time_swap_masks(n):
    if n < 0:
        raise ValueError("negative n")

    return range(
        1 << n
    )


def exact_upper_tail_p(
    observed,
    permutation_statistics,
):
    x = np.asarray(
        permutation_statistics,
        dtype=float,
    )

    if x.ndim != 1:
        raise ValueError(
            "permutation statistics must be one-dimensional"
        )

    if x.size == 0:
        raise ValueError(
            "empty permutation universe"
        )

    if not np.all(
        np.isfinite(x)
    ):
        raise ValueError(
            "non-finite permutation statistic"
        )

    observed = float(
        observed
    )

    tol = (
        FLOAT_TAIL_REL_TOL
        * max(
            1.0,
            abs(observed),
        )
    )

    count_ge = int(
        np.count_nonzero(
            x >= (
                observed - tol
            )
        )
    )

    # Exhaustive universe includes observed assignment.
    p = (
        count_ge
        / float(x.size)
    )

    return (
        p,
        count_ge,
        int(x.size),
        tol,
    )


def holm_adjust(p_values):
    p = np.asarray(
        p_values,
        dtype=float,
    )

    if p.ndim != 1:
        raise ValueError(
            "p-values must be one-dimensional"
        )

    if np.any(
        (p < 0.0)
        | (p > 1.0)
        | ~np.isfinite(p)
    ):
        raise ValueError(
            "invalid p-value"
        )

    m = p.size

    order = np.argsort(
        p,
        kind="mergesort"
    )

    adjusted_sorted = np.empty(
        m,
        dtype=float
    )

    running = 0.0

    for k, idx in enumerate(order):

        candidate = (
            (m - k)
            * p[idx]
        )

        running = max(
            running,
            candidate,
        )

        adjusted_sorted[k] = min(
            1.0,
            running,
        )

    adjusted = np.empty(
        m,
        dtype=float
    )

    for k, idx in enumerate(order):
        adjusted[idx] = (
            adjusted_sorted[k]
        )

    return adjusted


def unique_edge_labels(
    analyte_order=ANALYTE_ORDER,
):
    labels = []

    for j in range(
        len(analyte_order)
    ):
        for k in range(
            j + 1,
            len(analyte_order)
        ):
            labels.append(
                (
                    analyte_order[j],
                    analyte_order[k],
                )
            )

    return tuple(labels)
