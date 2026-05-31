"""Linear codes over finite fields."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, product
from math import prod
from typing import Iterator, Sequence

from .bounds import singleton_bound
from .enumerators import macwilliams_transform
from .exceptions import ComputationLimitError, LinearCodeError
from .fields import FiniteField
from .matrix import (
    Matrix,
    Vector,
    dot,
    hamming_distance,
    hamming_weight,
    mat_vec_left,
    matmul,
    normalize_matrix,
    normalize_vector,
    nullspace,
    rank,
    row_space_basis,
    transpose,
    vector_sub,
)


@dataclass(frozen=True)
class SystematicForm:
    """A systematic generator matrix plus its column permutation.

    ``permutation[j]`` is the original column index now located at column ``j``.
    """

    generator_matrix: Matrix
    permutation: list[int]
    information_set: list[int]


class LinearCode:
    """A linear code with an integer-encoded finite field."""

    def __init__(
        self,
        field: FiniteField,
        generator_matrix: Sequence[Sequence[int]],
        *,
        length: int | None = None,
        reduce: bool = True,
    ) -> None:
        self.field = field
        if length is not None and length < 0:
            raise LinearCodeError("length must be nonnegative")
        matrix = normalize_matrix(field, generator_matrix, ncols=length)
        if length is None:
            length = len(matrix[0]) if matrix else 0
        self._length = length
        self._generator_matrix = (
            row_space_basis(field, matrix, ncols=length) if reduce else [row[:] for row in matrix]
        )
        self._parity_check_matrix: Matrix | None = None

    @classmethod
    def from_parity_check(
        cls, field: FiniteField, parity_check_matrix: Sequence[Sequence[int]], *, length: int | None = None
    ) -> "LinearCode":
        h = normalize_matrix(field, parity_check_matrix, ncols=length)
        if length is None:
            length = len(h[0]) if h else 0
        return cls(field, nullspace(field, h, ncols=length), length=length)

    @property
    def generator_matrix(self) -> Matrix:
        return [row[:] for row in self._generator_matrix]

    @property
    def parity_check_matrix(self) -> Matrix:
        if self._parity_check_matrix is None:
            self._parity_check_matrix = nullspace(
                self.field, self._generator_matrix, ncols=self.length
            )
        return [row[:] for row in self._parity_check_matrix]

    @property
    def length(self) -> int:
        return self._length

    @property
    def dimension(self) -> int:
        return len(self._generator_matrix)

    @property
    def redundancy(self) -> int:
        return self.length - self.dimension

    @property
    def cardinality(self) -> int:
        return self.field.q**self.dimension

    def __repr__(self) -> str:
        return f"LinearCode([{self.length}, {self.dimension}] over {self.field!r})"

    def parameters(self, *, max_codewords: int = 1_000_000) -> tuple[int, int, int]:
        return (self.length, self.dimension, self.minimum_distance(max_codewords=max_codewords))

    def encode(self, message: Sequence[int]) -> Vector:
        msg = normalize_vector(self.field, message, length=self.dimension)
        return mat_vec_left(self.field, msg, self._generator_matrix)

    def syndrome(self, word: Sequence[int]) -> Vector:
        w = normalize_vector(self.field, word, length=self.length)
        return [dot(self.field, row, w) for row in self.parity_check_matrix]

    def contains(self, word: Sequence[int]) -> bool:
        return all(value == 0 for value in self.syndrome(word))

    def dual(self) -> "LinearCode":
        return LinearCode(self.field, self.parity_check_matrix, length=self.length, reduce=False)

    def hull(self) -> "LinearCode":
        if self.dimension == 0:
            return LinearCode(self.field, [], length=self.length)
        gram = matmul(self.field, self._generator_matrix, transpose(self._generator_matrix))
        message_basis = nullspace(self.field, gram, ncols=self.dimension)
        rows = [mat_vec_left(self.field, message, self._generator_matrix) for message in message_basis]
        return LinearCode(self.field, rows, length=self.length)

    def is_self_orthogonal(self) -> bool:
        for row in self._generator_matrix:
            for other in self._generator_matrix:
                if dot(self.field, row, other) != 0:
                    return False
        return True

    def is_self_dual(self) -> bool:
        return 2 * self.dimension == self.length and self.is_self_orthogonal()

    def is_mds(self, *, max_codewords: int = 1_000_000) -> bool:
        if self.dimension == 0:
            return False
        return self.minimum_distance(max_codewords=max_codewords) == singleton_bound(
            self.length, self.dimension
        )

    def codewords(self, *, max_codewords: int | None = 1_000_000) -> Iterator[Vector]:
        total = self.cardinality
        if max_codewords is not None and total > max_codewords:
            raise ComputationLimitError(
                f"exact enumeration would visit {total} codewords; raise max_codewords to allow it"
            )
        for message in product(range(self.field.q), repeat=self.dimension):
            yield self.encode(message)

    def minimum_weight_codeword(
        self,
        *,
        max_codewords: int = 1_000_000,
        algorithm: str = "auto",
        max_supports: int = 1_000_000,
    ) -> tuple[int, Vector]:
        if algorithm not in {"auto", "codewords", "parity_check"}:
            raise ValueError("algorithm must be 'auto', 'codewords', or 'parity_check'")
        if algorithm == "parity_check":
            return self.minimum_weight_codeword_by_parity_check(max_supports=max_supports)
        if algorithm == "auto" and self.cardinality - 1 > max_codewords:
            return self.minimum_weight_codeword_by_parity_check(max_supports=max_supports)

        if self.dimension == 0:
            return 0, [0] * self.length
        total = self.cardinality - 1
        if total > max_codewords:
            raise ComputationLimitError(
                f"exact minimum-distance search would visit {total} nonzero words; "
                "use a higher max_codewords or add a specialized algorithm"
            )

        best_weight = self.length + 1
        best_word: Vector | None = None
        for message in product(range(self.field.q), repeat=self.dimension):
            if all(x == 0 for x in message):
                continue
            word = self.encode(message)
            weight = hamming_weight(word)
            if 0 < weight < best_weight:
                best_weight = weight
                best_word = word
                if best_weight == 1:
                    break
        assert best_word is not None
        return best_weight, best_word

    def minimum_weight_codeword_by_parity_check(
        self, *, max_weight: int | None = None, max_supports: int = 1_000_000
    ) -> tuple[int, Vector]:
        """Find a minimum word via dependent column sets of a parity-check matrix.

        This is exact and can be dramatically faster than enumerating all
        codewords when the dimension is large and the target distance is small.
        """

        if self.dimension == 0:
            return 0, [0] * self.length
        if max_weight is None:
            max_weight = singleton_bound(self.length, self.dimension)
        h = self.parity_check_matrix
        checked = 0
        for weight in range(1, max_weight + 1):
            for support in combinations(range(self.length), weight):
                checked += 1
                if checked > max_supports:
                    raise ComputationLimitError(
                        f"parity-check search exceeded {max_supports} supports"
                    )
                submatrix = [[row[j] for j in support] for row in h]
                if rank(self.field, submatrix, ncols=weight) == weight:
                    continue
                dependencies = nullspace(self.field, submatrix, ncols=weight)
                if not dependencies:
                    continue
                word = [0] * self.length
                for index, value in zip(support, dependencies[0]):
                    word[index] = value
                actual_weight = hamming_weight(word)
                if actual_weight:
                    return actual_weight, word
        raise ComputationLimitError(f"no nonzero codeword found up to weight {max_weight}")

    def minimum_distance(
        self,
        *,
        max_codewords: int = 1_000_000,
        algorithm: str = "auto",
        max_supports: int = 1_000_000,
    ) -> int:
        return self.minimum_weight_codeword(
            max_codewords=max_codewords, algorithm=algorithm, max_supports=max_supports
        )[0]

    def weight_distribution(self, *, max_codewords: int = 1_000_000) -> dict[int, int]:
        distribution = {weight: 0 for weight in range(self.length + 1)}
        for word in self.codewords(max_codewords=max_codewords):
            distribution[hamming_weight(word)] += 1
        return distribution

    def dual_weight_distribution(self, *, max_codewords: int = 1_000_000) -> dict[int, int]:
        """Return the dual distribution using MacWilliams, not dual enumeration."""

        return macwilliams_transform(
            self.weight_distribution(max_codewords=max_codewords),
            q=self.field.q,
            length=self.length,
            dimension=self.dimension,
        )

    def distance(self, left: Sequence[int], right: Sequence[int]) -> int:
        return hamming_distance(
            normalize_vector(self.field, left, length=self.length),
            normalize_vector(self.field, right, length=self.length),
        )

    def puncture(self, positions: Sequence[int]) -> "LinearCode":
        pos = _validated_positions(positions, self.length)
        keep = [j for j in range(self.length) if j not in pos]
        rows = [[row[j] for j in keep] for row in self._generator_matrix]
        return LinearCode(self.field, rows, length=len(keep))

    def shorten(self, positions: Sequence[int]) -> "LinearCode":
        pos = _validated_positions(positions, self.length)
        if not pos:
            return LinearCode(self.field, self._generator_matrix, length=self.length)
        if self.dimension == 0:
            return LinearCode(self.field, [], length=self.length - len(pos))

        constraints = [[self._generator_matrix[row][col] for row in range(self.dimension)] for col in pos]
        message_basis = nullspace(self.field, constraints, ncols=self.dimension)
        keep = [j for j in range(self.length) if j not in pos]
        rows = []
        for message in message_basis:
            word = mat_vec_left(self.field, message, self._generator_matrix)
            rows.append([word[j] for j in keep])
        return LinearCode(self.field, rows, length=len(keep))

    def extend_with_parity(self) -> "LinearCode":
        rows: Matrix = []
        for row in self._generator_matrix:
            parity = self.field.neg(sum_field(self.field, row))
            rows.append([*row, parity])
        return LinearCode(self.field, rows, length=self.length + 1, reduce=False)

    def systematic_form(self) -> SystematicForm:
        if self.dimension == 0:
            return SystematicForm([], list(range(self.length)), [])

        a = self.generator_matrix
        k = self.dimension
        n = self.length
        permutation = list(range(n))
        pivot_row = 0

        for pivot_col in range(k):
            pivot = None
            for r in range(pivot_row, k):
                for c in range(pivot_col, n):
                    if a[r][c] != 0:
                        pivot = (r, c)
                        break
                if pivot is not None:
                    break
            if pivot is None:
                break

            r, c = pivot
            if r != pivot_row:
                a[pivot_row], a[r] = a[r], a[pivot_row]
            if c != pivot_col:
                for row in a:
                    row[pivot_col], row[c] = row[c], row[pivot_col]
                permutation[pivot_col], permutation[c] = permutation[c], permutation[pivot_col]

            inv = self.field.inv(a[pivot_row][pivot_col])
            a[pivot_row] = [self.field.mul(inv, x) for x in a[pivot_row]]
            for r2 in range(k):
                if r2 == pivot_row or a[r2][pivot_col] == 0:
                    continue
                factor = a[r2][pivot_col]
                a[r2] = [
                    self.field.sub(x, self.field.mul(factor, y))
                    for x, y in zip(a[r2], a[pivot_row])
                ]
            pivot_row += 1

        if pivot_row != k:
            raise LinearCodeError("generator matrix unexpectedly lost rank")
        return SystematicForm(a, permutation, permutation[:k])

    def coset_leader(self, syndrome: Sequence[int], *, max_weight: int | None = None) -> Vector:
        target = normalize_vector(self.field, syndrome, length=self.redundancy)
        if all(x == 0 for x in target):
            return [0] * self.length
        if max_weight is None:
            max_weight = self.length
        nonzero = list(self.field.nonzero_elements())
        for weight in range(1, max_weight + 1):
            for support in combinations(range(self.length), weight):
                for values in product(nonzero, repeat=weight):
                    error = [0] * self.length
                    for pos, value in zip(support, values):
                        error[pos] = value
                    if self.syndrome(error) == target:
                        return error
        raise ComputationLimitError(f"no coset leader found up to weight {max_weight}")

    def syndrome_decode(
        self, received: Sequence[int], *, max_weight: int | None = None
    ) -> tuple[Vector, Vector]:
        r = normalize_vector(self.field, received, length=self.length)
        error = self.coset_leader(self.syndrome(r), max_weight=max_weight)
        return vector_sub(self.field, r, error), error

    def syndrome_table(
        self, *, max_weight: int | None = None, max_errors: int = 1_000_000
    ) -> dict[tuple[int, ...], Vector]:
        """Return syndrome-to-coset-leader entries in increasing error weight."""

        if max_weight is None:
            max_weight = self.length
        target_size = self.field.q**self.redundancy
        table: dict[tuple[int, ...], Vector] = {tuple([0] * self.redundancy): [0] * self.length}
        nonzero = list(self.field.nonzero_elements())
        visited = 0
        for weight in range(1, max_weight + 1):
            for support in combinations(range(self.length), weight):
                for values in product(nonzero, repeat=weight):
                    visited += 1
                    if visited > max_errors:
                        raise ComputationLimitError(
                            f"syndrome-table construction exceeded {max_errors} errors"
                        )
                    error = [0] * self.length
                    for pos, value in zip(support, values):
                        error[pos] = value
                    table.setdefault(tuple(self.syndrome(error)), error)
                    if len(table) == target_size:
                        return table
        return table

    def covering_radius(self, *, max_weight: int | None = None, max_errors: int = 1_000_000) -> int:
        table = self.syndrome_table(max_weight=max_weight, max_errors=max_errors)
        if len(table) != self.field.q**self.redundancy:
            raise ComputationLimitError("syndrome table did not cover all cosets")
        return max(hamming_weight(error) for error in table.values())

    # Magma/Sage-style aliases for users porting notebooks or scripts.
    GeneratorMatrix = generator_matrix
    ParityCheckMatrix = parity_check_matrix
    Dual = dual
    Hull = hull
    MinimumDistance = minimum_distance
    WeightDistribution = weight_distribution
    DualWeightDistribution = dual_weight_distribution
    Syndrome = syndrome


def sum_field(field: FiniteField, values: Sequence[int]) -> int:
    total = field.zero
    for value in values:
        total = field.add(total, value)
    return total


def _validated_positions(positions: Sequence[int], length: int) -> set[int]:
    pos = {int(p) for p in positions}
    if any(p < 0 or p >= length for p in pos):
        raise LinearCodeError(f"positions must be in range(0, {length})")
    return pos


def direct_sum(left: LinearCode, right: LinearCode) -> LinearCode:
    if left.field != right.field:
        raise LinearCodeError("direct sum requires both codes to use the same field")
    rows: Matrix = []
    for row in left.generator_matrix:
        rows.append([*row, *([0] * right.length)])
    for row in right.generator_matrix:
        rows.append([*([0] * left.length), *row])
    return LinearCode(left.field, rows, length=left.length + right.length)


def expected_codewords(code: LinearCode) -> int:
    return prod([code.field.q] * code.dimension)
