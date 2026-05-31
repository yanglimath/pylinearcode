"""Linear algebra over :mod:`pylinearcode.fields` finite fields."""

from __future__ import annotations

from typing import Sequence

from .fields import FiniteField

Matrix = list[list[int]]
Vector = list[int]


def normalize_vector(field: FiniteField, vector: Sequence[int], *, length: int | None = None) -> Vector:
    out = [field.coerce(x) for x in vector]
    if length is not None and len(out) != length:
        raise ValueError(f"expected vector of length {length}, got {len(out)}")
    return out


def normalize_matrix(
    field: FiniteField, matrix: Sequence[Sequence[int]], *, ncols: int | None = None
) -> Matrix:
    rows = [normalize_vector(field, row) for row in matrix]
    if ncols is None:
        if rows:
            ncols = len(rows[0])
        else:
            return []
    for row in rows:
        if len(row) != ncols:
            raise ValueError("matrix rows must all have the same length")
    return rows


def zeros(nrows: int, ncols: int) -> Matrix:
    return [[0 for _ in range(ncols)] for _ in range(nrows)]


def identity(field: FiniteField, size: int) -> Matrix:
    return [[field.one if i == j else field.zero for j in range(size)] for i in range(size)]


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    if not matrix:
        return []
    return [list(col) for col in zip(*matrix)]


def dot(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != len(right):
        raise ValueError("dot product vectors must have the same length")
    total = field.zero
    for a, b in zip(left, right):
        total = field.add(total, field.mul(a, b))
    return total


def vector_add(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Vector:
    if len(left) != len(right):
        raise ValueError("vectors must have the same length")
    return [field.add(a, b) for a, b in zip(left, right)]


def vector_sub(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Vector:
    if len(left) != len(right):
        raise ValueError("vectors must have the same length")
    return [field.sub(a, b) for a, b in zip(left, right)]


def scalar_mul(field: FiniteField, scalar: int, vector: Sequence[int]) -> Vector:
    scalar = field.coerce(scalar)
    if scalar == 0:
        return [0] * len(vector)
    return [field.mul(scalar, x) for x in vector]


def matmul(field: FiniteField, left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> Matrix:
    if not left:
        return []
    if not right:
        raise ValueError("right matrix must not be empty")
    left_cols = len(left[0])
    right_rows = len(right)
    if left_cols != right_rows:
        raise ValueError("incompatible matrix dimensions")
    right_t = transpose(right)
    return [[dot(field, row, col) for col in right_t] for row in left]


def mat_vec_left(field: FiniteField, vector: Sequence[int], matrix: Sequence[Sequence[int]]) -> Vector:
    """Return ``vector * matrix``."""

    if not matrix:
        return []
    if len(vector) != len(matrix):
        raise ValueError("left vector length must equal matrix row count")
    ncols = len(matrix[0])
    out = [field.zero] * ncols
    for coeff, row in zip(vector, matrix):
        coeff = field.coerce(coeff)
        if coeff == 0:
            continue
        for j, value in enumerate(row):
            out[j] = field.add(out[j], field.mul(coeff, value))
    return out


def rref(field: FiniteField, matrix: Sequence[Sequence[int]], *, ncols: int | None = None) -> tuple[Matrix, list[int]]:
    """Reduced row-echelon form and pivot columns."""

    a = normalize_matrix(field, matrix, ncols=ncols)
    if not a:
        return [], []
    row_count = len(a)
    col_count = len(a[0])
    pivot_row = 0
    pivots: list[int] = []

    for col in range(col_count):
        pivot = None
        for r in range(pivot_row, row_count):
            if a[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue

        if pivot != pivot_row:
            a[pivot_row], a[pivot] = a[pivot], a[pivot_row]

        inv = field.inv(a[pivot_row][col])
        a[pivot_row] = [field.mul(inv, x) for x in a[pivot_row]]

        for r in range(row_count):
            if r == pivot_row or a[r][col] == 0:
                continue
            factor = a[r][col]
            a[r] = [field.sub(x, field.mul(factor, y)) for x, y in zip(a[r], a[pivot_row])]

        pivots.append(col)
        pivot_row += 1
        if pivot_row == row_count:
            break

    return a, pivots


def rank(field: FiniteField, matrix: Sequence[Sequence[int]], *, ncols: int | None = None) -> int:
    return len(rref(field, matrix, ncols=ncols)[1])


def row_space_basis(
    field: FiniteField, matrix: Sequence[Sequence[int]], *, ncols: int | None = None
) -> Matrix:
    reduced, _ = rref(field, matrix, ncols=ncols)
    return [row for row in reduced if any(x != 0 for x in row)]


def nullspace(
    field: FiniteField, matrix: Sequence[Sequence[int]], *, ncols: int | None = None
) -> Matrix:
    """Basis for the right nullspace ``{x : matrix * x^T = 0}``."""

    if ncols is None:
        if matrix:
            ncols = len(matrix[0])
        else:
            raise ValueError("ncols is required for the nullspace of an empty matrix")

    reduced, pivots = rref(field, matrix, ncols=ncols)
    pivot_set = set(pivots)
    free_cols = [col for col in range(ncols) if col not in pivot_set]
    if not free_cols:
        return []

    basis: Matrix = []
    for free_col in free_cols:
        vector = [field.zero] * ncols
        vector[free_col] = field.one
        for row_index, pivot_col in enumerate(pivots):
            vector[pivot_col] = field.neg(reduced[row_index][free_col])
        basis.append(vector)
    return basis


def orthogonal_complement(
    field: FiniteField, rows: Sequence[Sequence[int]], *, ncols: int | None = None
) -> Matrix:
    return nullspace(field, rows, ncols=ncols)


def hamming_weight(vector: Sequence[int]) -> int:
    return sum(1 for value in vector if value != 0)


def hamming_distance(left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != len(right):
        raise ValueError("vectors must have the same length")
    return sum(1 for a, b in zip(left, right) if a != b)

