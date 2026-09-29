
import itertools
import numpy as np


def _vector(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    if x.ndim != 1:
        raise ValueError(
            "input must be 1D"
        )

    if not np.all(
        np.isfinite(x)
    ):
        raise ValueError(
            "non-finite input"
        )

    return x


def subject_log_response(
    a,
    b,
):

    a = _vector(a)
    b = _vector(b)

    if a.shape != b.shape:
        raise ValueError(
            "A/B shape mismatch"
        )

    if np.any(a <= 0) or np.any(b <= 0):
        raise ValueError(
            "log response requires positive values"
        )

    return np.log(
        b / a
    )


def within_group_statistic(
    log_response,
):

    x = _vector(
        log_response
    )

    signed_mean = float(
        np.mean(x)
    )

    statistic = abs(
        signed_mean
    )

    return (
        float(statistic),
        signed_mean,
    )


def between_group_statistic(
    log_response,
    group1_indices,
):

    x = _vector(
        log_response
    )

    n = len(x)

    group1 = sorted(
        set(
            int(i)
            for i in group1_indices
        )
    )

    if len(group1) == 0 or len(group1) >= n:
        raise ValueError(
            "invalid group size"
        )

    group2 = [
        i
        for i in range(n)
        if i not in group1
    ]

    mean1 = float(
        np.mean(
            x[group1]
        )
    )

    mean2 = float(
        np.mean(
            x[group2]
        )
    )

    difference = (
        mean1
        - mean2
    )

    statistic = abs(
        difference
    )

    return (
        float(statistic),
        mean1,
        mean2,
        float(difference),
    )


def sign_flip_values(
    x,
    mask,
):

    x = _vector(x)

    n = len(x)

    if mask < 0 or mask >= (1 << n):
        raise ValueError(
            "mask outside space"
        )

    out = x.copy()

    for i in range(n):

        if (mask >> i) & 1:
            out[i] *= -1.0

    return out


def iter_sign_masks(n):

    return range(
        1 << n
    )


def iter_equal_group_combinations(
    n=20,
    group_size=10,
):

    return itertools.combinations(
        range(n),
        group_size,
    )


def exact_upper_tail_p(
    observed,
    values,
):

    values = _vector(
        values
    )

    observed = float(
        observed
    )

    tol = (
        1e-12
        * max(
            1.0,
            abs(observed),
        )
    )

    count = int(
        np.count_nonzero(
            values >= (
                observed - tol
            )
        )
    )

    return (
        count / float(len(values)),
        count,
        len(values),
        tol,
    )


def holm_adjust(
    p_values,
):

    p = _vector(
        p_values
    )

    if np.any(
        (p < 0.0)
        | (p > 1.0)
    ):
        raise ValueError(
            "invalid p-value"
        )

    m = len(p)

    order = np.argsort(
        p,
        kind="mergesort",
    )

    adjusted_sorted = np.empty(
        m,
        dtype=float,
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
        dtype=float,
    )

    for k, idx in enumerate(order):

        adjusted[idx] = (
            adjusted_sorted[k]
        )

    return adjusted


def robust_summary(
    x,
):

    x = _vector(x)

    q1 = float(
        np.quantile(
            x,
            0.25,
            method="linear",
        )
    )

    median = float(
        np.quantile(
            x,
            0.50,
            method="linear",
        )
    )

    q3 = float(
        np.quantile(
            x,
            0.75,
            method="linear",
        )
    )

    return (
        median,
        q1,
        q3,
        q3 - q1,
    )
