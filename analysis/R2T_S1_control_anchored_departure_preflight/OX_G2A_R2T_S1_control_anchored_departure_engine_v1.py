
import numpy as np


def _matrix(x):

    x = np.asarray(
        x,
        dtype=float,
    )

    if x.ndim != 2:
        raise ValueError("matrix must be 2D")

    if not np.all(np.isfinite(x)):
        raise ValueError("matrix contains non-finite values")

    return x


def oas_covariance(x):

    x = _matrix(x)

    n, p = x.shape

    if n < 2:
        raise ValueError("at least two observations required")

    center = x.mean(axis=0)

    xc = x - center

    empirical = (
        xc.T @ xc
    ) / float(n)

    mu = float(
        np.trace(empirical)
        / p
    )

    alpha = float(
        np.mean(
            empirical ** 2
        )
    )

    denominator = (
        (n + 1.0)
        * (
            alpha
            - (
                mu ** 2
                / p
            )
        )
    )

    numerator = (
        alpha
        + mu ** 2
    )

    if denominator <= 0.0:
        shrinkage = 1.0
    else:
        shrinkage = min(
            numerator
            / denominator,
            1.0,
        )

    covariance = (
        (1.0 - shrinkage)
        * empirical
    )

    covariance = covariance.copy()

    covariance.flat[
        :: p + 1
    ] += (
        shrinkage
        * mu
    )

    eig = np.linalg.eigvalsh(
        covariance
    )

    if not np.all(
        eig > 0.0
    ):
        raise ValueError(
            "OAS covariance not positive definite"
        )

    precision = np.linalg.inv(
        covariance
    )

    return (
        center,
        covariance,
        precision,
        float(shrinkage),
    )


def mahalanobis_squared(
    x,
    center,
    precision,
):

    x = _matrix(x)

    center = np.asarray(
        center,
        dtype=float,
    )

    precision = np.asarray(
        precision,
        dtype=float,
    )

    if center.ndim != 1:
        raise ValueError(
            "center must be 1D"
        )

    p = x.shape[1]

    if center.shape != (p,):
        raise ValueError(
            "center dimension mismatch"
        )

    if precision.shape != (p, p):
        raise ValueError(
            "precision dimension mismatch"
        )

    delta = (
        x - center
    )

    d2 = np.einsum(
        "ij,jk,ik->i",
        delta,
        precision,
        delta,
    )

    return d2


def departure_change_statistic(
    za,
    zb,
    center,
    precision,
):

    za = _matrix(za)
    zb = _matrix(zb)

    if za.shape != zb.shape:
        raise ValueError(
            "A/B shapes differ"
        )

    d2_a = mahalanobis_squared(
        za,
        center,
        precision,
    )

    d2_b = mahalanobis_squared(
        zb,
        center,
        precision,
    )

    delta_d2 = (
        d2_b - d2_a
    )

    mean_delta = float(
        np.mean(
            delta_d2
        )
    )

    statistic = abs(
        mean_delta
    )

    return (
        float(statistic),
        mean_delta,
        d2_a,
        d2_b,
        delta_d2,
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
            "values must be 1D"
        )

    if not np.all(
        np.isfinite(values)
    ):
        raise ValueError(
            "non-finite permutation statistic"
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
