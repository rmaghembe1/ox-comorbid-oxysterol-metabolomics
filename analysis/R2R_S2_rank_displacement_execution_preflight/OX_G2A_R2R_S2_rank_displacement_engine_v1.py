
import itertools
import math
import numpy as np


def _matrix(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    if x.ndim != 2:
        raise ValueError(
            "matrix must be 2D"
        )

    if not np.all(
        np.isfinite(x)
    ):
        raise ValueError(
            "matrix contains non-finite values"
        )

    return x


def average_ranks_1d(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    if x.ndim != 1:
        raise ValueError(
            "rank input must be one-dimensional"
        )

    n = len(x)

    order = np.argsort(
        x,
        kind="mergesort",
    )

    ranks = np.empty(
        n,
        dtype=float,
    )

    i = 0

    while i < n:

        j = i + 1

        while (
            j < n
            and
            x[order[j]]
            == x[order[i]]
        ):
            j += 1

        rank = (
            (i + 1 + j)
            / 2.0
        )

        ranks[
            order[i:j]
        ] = rank

        i = j

    return ranks


def rank_columns(x):

    x = _matrix(x)

    return np.column_stack([
        average_ranks_1d(
            x[:, j]
        )
        for j in range(
            x.shape[1]
        )
    ])


def longitudinal_rank_displacement(
    a,
    b,
):

    a = _matrix(a)
    b = _matrix(b)

    if a.shape != b.shape:
        raise ValueError(
            "A/B shapes differ"
        )

    rank_a = rank_columns(a)
    rank_b = rank_columns(b)

    delta = (
        rank_b
        - rank_a
    )

    return (
        delta,
        rank_a,
        rank_b,
    )


def between_group_statistic(
    displacement,
    group1_indices,
):

    displacement = _matrix(
        displacement
    )

    n = displacement.shape[0]

    group1 = set(
        int(x)
        for x in group1_indices
    )

    if len(group1) == 0:
        raise ValueError(
            "empty group1"
        )

    if len(group1) >= n:
        raise ValueError(
            "group1 contains all subjects"
        )

    if any(
        x < 0 or x >= n
        for x in group1
    ):
        raise ValueError(
            "group index out of range"
        )

    group2 = [
        i
        for i in range(n)
        if i not in group1
    ]

    group1 = sorted(
        group1
    )

    mean1 = displacement[
        group1,
        :
    ].mean(
        axis=0
    )

    mean2 = displacement[
        group2,
        :
    ].mean(
        axis=0
    )

    difference = (
        mean1 - mean2
    )

    statistic = float(
        np.sum(
            difference ** 2
        )
    )

    return (
        statistic,
        mean1,
        mean2,
        difference,
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

    values = np.asarray(
        values,
        dtype=float,
    )

    if values.ndim != 1:
        raise ValueError(
            "permutation values must be 1D"
        )

    if not np.all(
        np.isfinite(values)
    ):
        raise ValueError(
            "non-finite permutation values"
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
