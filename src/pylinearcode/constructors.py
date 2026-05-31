"""Constructors for common linear codes."""

from __future__ import annotations

import random
from itertools import combinations, product
from typing import Sequence

from .fields import FiniteField
from .linear_code import LinearCode
from .matrix import Matrix, identity, mat_vec_left, nullspace, rank
from .polynomials import lcm as poly_lcm
from .polynomials import monic as poly_monic
from .polynomials import mul as poly_mul


def zero_code(field: FiniteField, length: int) -> LinearCode:
    return LinearCode(field, [], length=length)


def whole_space_code(field: FiniteField, length: int) -> LinearCode:
    return LinearCode(field, identity(field, length), length=length, reduce=False)


def repetition_code(field: FiniteField, length: int) -> LinearCode:
    if length < 0:
        raise ValueError("length must be nonnegative")
    if length == 0:
        return zero_code(field, 0)
    return LinearCode(field, [[field.one] * length], length=length, reduce=False)


def parity_check_code(field: FiniteField, length: int) -> LinearCode:
    """Single parity-check code ``sum_i c_i = 0``."""

    if length < 1:
        raise ValueError("length must be positive")
    return LinearCode.from_parity_check(field, [[field.one] * length], length=length)


def hamming_code(field: FiniteField, r: int) -> LinearCode:
    """q-ary Hamming code with an ``r x n`` parity-check matrix."""

    if r < 2:
        raise ValueError("r must be at least 2")
    columns = _projective_points(field, r)
    h = [[col[i] for col in columns] for i in range(r)]
    return LinearCode.from_parity_check(field, h, length=len(columns))


def extended_hamming_code(field: FiniteField, r: int) -> LinearCode:
    """Return the Hamming code with one overall parity coordinate appended."""

    return hamming_code(field, r).extend_with_parity()


def simplex_code(field: FiniteField, r: int) -> LinearCode:
    """Return the dual of the q-ary Hamming code."""

    return hamming_code(field, r).dual()


def reed_solomon_code(
    field: FiniteField,
    *,
    length: int,
    dimension: int,
    evaluation_points: Sequence[int] | None = None,
    multipliers: Sequence[int] | None = None,
) -> LinearCode:
    """Generalized Reed-Solomon code.

    The generator row ``j`` is ``(v_i * alpha_i**j)_i`` for ``0 <= j < k``.
    """

    if not 0 <= dimension <= length:
        raise ValueError("dimension must satisfy 0 <= k <= n")
    if length > field.q:
        raise ValueError("this constructor requires length <= q")
    if evaluation_points is None:
        evaluation_points = list(field.elements())[:length]
    if len(set(evaluation_points)) != length:
        raise ValueError("evaluation points must be distinct")
    points = [field.coerce(a) for a in evaluation_points]

    if multipliers is None:
        multipliers = [field.one] * length
    if len(multipliers) != length:
        raise ValueError("multipliers must have length n")
    multipliers = [field.coerce(v) for v in multipliers]
    if any(v == 0 for v in multipliers):
        raise ValueError("multipliers must be nonzero")

    rows: Matrix = []
    for j in range(dimension):
        rows.append([field.mul(v, field.pow(alpha, j)) for alpha, v in zip(points, multipliers)])
    return LinearCode(field, rows, length=length, reduce=False)


def reed_muller_code(order: int, variables: int) -> LinearCode:
    """Binary Reed-Muller code RM(order, variables)."""

    if variables < 0:
        raise ValueError("variables must be nonnegative")
    field = FiniteField(2)
    length = 2**variables
    if order < 0:
        return zero_code(field, length)
    if order >= variables:
        return whole_space_code(field, length)

    points = list(product((0, 1), repeat=variables))
    rows: Matrix = []
    for degree in range(order + 1):
        for support in combinations(range(variables), degree):
            row = []
            for point in points:
                value = 1
                for idx in support:
                    value &= point[idx]
                row.append(value)
            rows.append(row)
    return LinearCode(field, rows, length=length, reduce=False)


def bch_generator_polynomial(
    field: FiniteField,
    *,
    extension_degree: int,
    designed_distance: int,
    offset: int = 1,
) -> list[int]:
    """Generator polynomial for a primitive narrow-sense BCH code.

    The current implementation supports prime base fields ``GF(p)`` and length
    ``p^m - 1``. The roots are ``alpha^offset`` through
    ``alpha^(offset + designed_distance - 2)`` in ``GF(p^m)``.
    """

    if field.degree != 1:
        raise NotImplementedError("BCH construction currently supports prime base fields GF(p)")
    if extension_degree < 1:
        raise ValueError("extension_degree must be positive")
    if designed_distance < 2:
        raise ValueError("designed_distance must be at least 2")

    extension = FiniteField(field.p, extension_degree, table_limit=field.table_limit)
    length = extension.q - 1
    alpha = extension.primitive_element()
    generator = [field.one]
    used_exponents: set[int] = set()

    for exponent in range(offset, offset + designed_distance - 1):
        exponent %= length
        if exponent in used_exponents:
            continue
        coset = _cyclotomic_coset(field.q, length, exponent)
        used_exponents.update(coset)
        minimal = _minimal_polynomial_from_coset(field, extension, alpha, coset)
        generator = poly_lcm(field, generator, minimal)
    return poly_monic(field, generator)


def bch_code(
    field: FiniteField,
    *,
    extension_degree: int,
    designed_distance: int,
    length: int | None = None,
    offset: int = 1,
) -> LinearCode:
    """Primitive narrow-sense BCH code over a prime field."""

    full_length = field.q**extension_degree - 1
    if length is None:
        length = full_length
    if length != full_length:
        raise NotImplementedError("only primitive BCH length q^m - 1 is currently supported")
    generator = bch_generator_polynomial(
        field,
        extension_degree=extension_degree,
        designed_distance=designed_distance,
        offset=offset,
    )
    return cyclic_code_from_generator(field, length=length, generator_polynomial=generator)


def random_linear_code(
    field: FiniteField,
    *,
    length: int,
    dimension: int,
    seed: int | None = None,
    systematic: bool = False,
    max_attempts: int = 10_000,
) -> LinearCode:
    if not 0 <= dimension <= length:
        raise ValueError("dimension must satisfy 0 <= k <= n")
    rng = random.Random(seed)
    if systematic:
        rows = identity(field, dimension)
        for row in rows:
            row.extend(rng.randrange(field.q) for _ in range(length - dimension))
        return LinearCode(field, rows, length=length, reduce=False)

    rows: Matrix = []
    attempts = 0
    while len(rows) < dimension and attempts < max_attempts:
        attempts += 1
        candidate = [rng.randrange(field.q) for _ in range(length)]
        if rank(field, [*rows, candidate], ncols=length) > len(rows):
            rows.append(candidate)
    if len(rows) != dimension:
        raise RuntimeError("failed to generate full-rank random generator matrix")
    return LinearCode(field, rows, length=length, reduce=False)


def cyclic_code_from_generator(
    field: FiniteField, *, length: int, generator_polynomial: Sequence[int]
) -> LinearCode:
    """Cyclic code generated by a polynomial ``g(x)``.

    Coefficients are ascending: ``[g_0, g_1, ...]``. The caller is responsible
    for choosing a generator that divides ``x^n - 1`` when a mathematically
    canonical cyclic code is required.
    """

    if length < 0:
        raise ValueError("length must be nonnegative")
    coeffs = [field.coerce(c) for c in generator_polynomial]
    while coeffs and coeffs[-1] == 0:
        coeffs.pop()
    if not coeffs:
        return zero_code(field, length)
    degree = len(coeffs) - 1
    if degree > length:
        raise ValueError("generator degree must not exceed length")
    dimension = length - degree
    rows: Matrix = []
    for shift in range(dimension):
        row = [0] * length
        for i, coeff in enumerate(coeffs):
            row[(i + shift) % length] = field.add(row[(i + shift) % length], coeff)
        rows.append(row)
    return LinearCode(field, rows, length=length)


def _cyclotomic_coset(q: int, modulus: int, start: int) -> list[int]:
    coset: list[int] = []
    value = start % modulus
    while value not in coset:
        coset.append(value)
        value = (value * q) % modulus
    return coset


def _minimal_polynomial_from_coset(
    base_field: FiniteField, extension: FiniteField, alpha: int, coset: Sequence[int]
) -> list[int]:
    poly = [extension.one]
    for exponent in coset:
        root = extension.pow(alpha, exponent)
        poly = poly_mul(extension, poly, [extension.neg(root), extension.one])

    coeffs: list[int] = []
    for coeff in poly:
        if coeff >= base_field.q:
            raise ArithmeticError("minimal polynomial coefficient did not descend to base field")
        coeffs.append(base_field.coerce(coeff))
    return poly_monic(base_field, coeffs)


def _projective_points(field: FiniteField, r: int) -> list[tuple[int, ...]]:
    seen: set[tuple[int, ...]] = set()
    points: list[tuple[int, ...]] = []
    for vector in product(range(field.q), repeat=r):
        if all(x == 0 for x in vector):
            continue
        first = next(i for i, x in enumerate(vector) if x != 0)
        inv = field.inv(vector[first])
        canonical = tuple(field.mul(inv, x) for x in vector)
        if canonical not in seen:
            seen.add(canonical)
            points.append(canonical)
    return points


def generator_from_parity_check(field: FiniteField, parity_check_matrix: Sequence[Sequence[int]]) -> Matrix:
    if not parity_check_matrix:
        return []
    return nullspace(field, parity_check_matrix, ncols=len(parity_check_matrix[0]))


def encode_rows(field: FiniteField, messages: Sequence[Sequence[int]], generator: Matrix) -> Matrix:
    return [mat_vec_left(field, message, generator) for message in messages]
