"""Classical bounds for q-ary linear codes."""

from __future__ import annotations

from math import ceil, comb


def singleton_bound(length: int, dimension: int) -> int:
    if not 0 <= dimension <= length:
        raise ValueError("dimension must satisfy 0 <= k <= n")
    if dimension == 0:
        return length + 1
    return length - dimension + 1


def hamming_sphere_volume(q: int, length: int, radius: int) -> int:
    if q < 2:
        raise ValueError("q must be at least 2")
    if not 0 <= radius <= length:
        raise ValueError("radius must satisfy 0 <= t <= n")
    return sum(comb(length, i) * (q - 1) ** i for i in range(radius + 1))


def satisfies_hamming_bound(q: int, length: int, dimension: int, distance: int) -> bool:
    if not 0 <= dimension <= length:
        raise ValueError("dimension must satisfy 0 <= k <= n")
    t = (distance - 1) // 2
    return q**dimension * hamming_sphere_volume(q, length, t) <= q**length


def griesmer_bound_length(q: int, dimension: int, distance: int) -> int:
    if q < 2:
        raise ValueError("q must be at least 2")
    if dimension < 0:
        raise ValueError("dimension must be nonnegative")
    if distance < 0:
        raise ValueError("distance must be nonnegative")
    return sum(ceil(distance / (q**i)) for i in range(dimension))

