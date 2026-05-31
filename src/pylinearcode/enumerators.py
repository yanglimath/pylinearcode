"""Weight-enumerator utilities."""

from __future__ import annotations

from math import comb
from typing import Mapping


def krawtchouk(q: int, length: int, j: int, i: int) -> int:
    """Return the q-ary Krawtchouk polynomial ``K_j(i)``."""

    if q < 2:
        raise ValueError("q must be at least 2")
    if not 0 <= i <= length or not 0 <= j <= length:
        raise ValueError("indices must lie between 0 and length")
    total = 0
    for ell in range(j + 1):
        if ell <= i and j - ell <= length - i:
            total += (-1) ** ell * (q - 1) ** (j - ell) * comb(i, ell) * comb(length - i, j - ell)
    return total


def macwilliams_transform(
    distribution: Mapping[int, int], *, q: int, length: int, dimension: int
) -> dict[int, int]:
    """Return the dual weight distribution using the MacWilliams identity."""

    size = q**dimension
    dual: dict[int, int] = {}
    for j in range(length + 1):
        numerator = sum(
            count * krawtchouk(q, length, j, i) for i, count in distribution.items()
        )
        if numerator % size != 0:
            raise ArithmeticError("MacWilliams transform did not produce integral coefficients")
        dual[j] = numerator // size
    return dual

