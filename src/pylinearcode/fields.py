"""Finite fields represented by compact integers.

Elements of GF(p^m) are integers ``0 <= a < p**m``. For extension fields the
integer is interpreted as coefficients in base ``p``:

``a = c_0 + c_1*p + ... + c_{m-1}*p**(m-1)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from math import prod
from typing import Iterable, Sequence


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def _trim(poly: Sequence[int]) -> list[int]:
    result = [int(c) for c in poly]
    while result and result[-1] == 0:
        result.pop()
    return result or [0]


def _poly_add(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
    return _trim(out)


def _poly_sub(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p
    return _trim(out)


def _poly_mul(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    if a == [0] or b == [0]:
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % p
    return _trim(out)


def _poly_divmod(a: Sequence[int], b: Sequence[int], p: int) -> tuple[list[int], list[int]]:
    dividend = _trim(a)
    divisor = _trim(b)
    if divisor == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(dividend) < len(divisor):
        return [0], dividend

    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = dividend[:]
    inv_lead = pow(divisor[-1], -1, p)
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        coeff = remainder[-1] * inv_lead % p
        quotient[shift] = coeff
        for i, di in enumerate(divisor):
            remainder[shift + i] = (remainder[shift + i] - coeff * di) % p
        remainder = _trim(remainder)
    return _trim(quotient), remainder


def _poly_mod(a: Sequence[int], modulus: Sequence[int], p: int) -> list[int]:
    return _poly_divmod(a, modulus, p)[1]


def _poly_gcd(a: Sequence[int], b: Sequence[int], p: int) -> list[int]:
    left = _trim(a)
    right = _trim(b)
    while right != [0]:
        _, rem = _poly_divmod(left, right, p)
        left, right = right, rem
    inv = pow(left[-1], -1, p)
    return [(c * inv) % p for c in left]


def _poly_pow_mod(base: Sequence[int], exponent: int, modulus: Sequence[int], p: int) -> list[int]:
    result = [1]
    power = _poly_mod(base, modulus, p)
    e = exponent
    while e:
        if e & 1:
            result = _poly_mod(_poly_mul(result, power, p), modulus, p)
        power = _poly_mod(_poly_mul(power, power, p), modulus, p)
        e >>= 1
    return result


def _prime_factors(n: int) -> set[int]:
    factors: set[int] = set()
    d = 2
    remaining = n
    while d * d <= remaining:
        while remaining % d == 0:
            factors.add(d)
            remaining //= d
        d += 1 if d == 2 else 2
    if remaining > 1:
        factors.add(remaining)
    return factors


def _is_irreducible(poly: Sequence[int], p: int) -> bool:
    f = _trim([c % p for c in poly])
    degree = len(f) - 1
    if degree <= 0:
        return False
    if f[-1] != 1:
        return False
    x = [0, 1]
    for factor in _prime_factors(degree):
        test = _poly_sub(_poly_pow_mod(x, p ** (degree // factor), f, p), x, p)
        if _poly_gcd(f, test, p) != [1]:
            return False
    final = _poly_sub(_poly_pow_mod(x, p**degree, f, p), x, p)
    return _poly_mod(final, f, p) == [0]


def _find_irreducible_polynomial(p: int, degree: int) -> tuple[int, ...]:
    if degree == 1:
        return (0, 1)
    # Search monic polynomials c_0 + ... + c_{m-1}x^{m-1} + x^m.
    # Constant term 0 is reducible for degree > 1.
    for c0 in range(1, p):
        for rest in product(range(p), repeat=degree - 1):
            candidate = (c0, *rest, 1)
            if _is_irreducible(candidate, p):
                return candidate
    raise ValueError(f"could not find irreducible polynomial over GF({p}) of degree {degree}")


@dataclass(slots=True)
class FiniteField:
    """A finite field GF(p^m).

    Parameters
    ----------
    p:
        Prime characteristic.
    degree:
        Extension degree. ``1`` gives a prime field.
    modulus:
        Optional monic irreducible polynomial coefficients in ascending order.
    table_limit:
        Precompute operation tables when ``q <= table_limit``.
    """

    p: int
    degree: int = 1
    modulus: Sequence[int] | None = None
    table_limit: int = 2048

    def __post_init__(self) -> None:
        self.p = int(self.p)
        self.degree = int(self.degree)
        if not _is_prime(self.p):
            raise ValueError("field characteristic must be prime")
        if self.degree < 1:
            raise ValueError("field degree must be positive")
        self.q = self.p**self.degree
        if self.degree == 1:
            self.modulus = (0, 1)
        elif self.modulus is None:
            self.modulus = _find_irreducible_polynomial(self.p, self.degree)
        else:
            normalized = tuple(int(c) % self.p for c in self.modulus)
            if len(normalized) != self.degree + 1 or normalized[-1] != 1:
                raise ValueError("modulus must be monic of length degree + 1")
            if not _is_irreducible(normalized, self.p):
                raise ValueError("modulus must be irreducible over the prime field")
            self.modulus = normalized

        self._add_table: list[list[int]] | None = None
        self._mul_table: list[list[int]] | None = None
        self._neg_table: list[int] | None = None
        self._inv_table: list[int | None] | None = None
        if self.q <= self.table_limit:
            self._build_tables()

    q: int = 0
    _add_table: list[list[int]] | None = None
    _mul_table: list[list[int]] | None = None
    _neg_table: list[int] | None = None
    _inv_table: list[int | None] | None = None

    @property
    def zero(self) -> int:
        return 0

    @property
    def one(self) -> int:
        return 1

    def __call__(self, value: int) -> int:
        return self.coerce(value)

    def __repr__(self) -> str:
        if self.degree == 1:
            return f"GF({self.p})"
        coeffs = ",".join(str(c) for c in self.modulus)
        return f"GF({self.p}^{self.degree}, modulus=({coeffs}))"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, FiniteField):
            return NotImplemented
        return self.p == other.p and self.degree == other.degree and tuple(self.modulus) == tuple(
            other.modulus
        )

    def __hash__(self) -> int:
        return hash((self.p, self.degree, tuple(self.modulus)))

    def coerce(self, value: int) -> int:
        value = int(value)
        if self.degree == 1:
            return value % self.p
        if not 0 <= value < self.q:
            raise ValueError(f"extension-field element must be in range(0, {self.q})")
        return value

    def elements(self) -> range:
        return range(self.q)

    def nonzero_elements(self) -> range:
        return range(1, self.q)

    def add(self, a: int, b: int) -> int:
        a = self.coerce(a)
        b = self.coerce(b)
        if self._add_table is not None:
            return self._add_table[a][b]
        return self._add_raw(a, b)

    def sub(self, a: int, b: int) -> int:
        return self.add(a, self.neg(b))

    def neg(self, a: int) -> int:
        a = self.coerce(a)
        if self._neg_table is not None:
            return self._neg_table[a]
        return self._neg_raw(a)

    def mul(self, a: int, b: int) -> int:
        a = self.coerce(a)
        b = self.coerce(b)
        if self._mul_table is not None:
            return self._mul_table[a][b]
        return self._mul_raw(a, b)

    def inv(self, a: int) -> int:
        a = self.coerce(a)
        if a == 0:
            raise ZeroDivisionError("0 has no multiplicative inverse")
        if self._inv_table is not None:
            inverse = self._inv_table[a]
            assert inverse is not None
            return inverse
        return self.pow(a, self.q - 2)

    def div(self, a: int, b: int) -> int:
        return self.mul(a, self.inv(b))

    def pow(self, a: int, exponent: int) -> int:
        a = self.coerce(a)
        if exponent < 0:
            return self.pow(self.inv(a), -exponent)
        result = 1
        power = a
        e = exponent
        while e:
            if e & 1:
                result = self.mul(result, power)
            power = self.mul(power, power)
            e >>= 1
        return result

    def format(self, a: int, variable: str = "a") -> str:
        a = self.coerce(a)
        if self.degree == 1:
            return str(a)
        coeffs = self._digits(a)
        terms: list[str] = []
        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue
            if i == 0:
                term = str(coeff)
            elif i == 1:
                term = variable if coeff == 1 else f"{coeff}{variable}"
            else:
                term = f"{variable}^{i}" if coeff == 1 else f"{coeff}{variable}^{i}"
            terms.append(term)
        return "0" if not terms else " + ".join(terms)

    def primitive_element(self) -> int:
        """Return a generator of the multiplicative group."""

        if self.q == 2:
            return 1
        order = self.q - 1
        factors = _prime_factors(order)
        for candidate in self.nonzero_elements():
            if all(self.pow(candidate, order // factor) != 1 for factor in factors):
                return candidate
        raise RuntimeError("finite field multiplicative group has no generator")

    def _build_tables(self) -> None:
        self._add_table = [[self._add_raw(a, b) for b in range(self.q)] for a in range(self.q)]
        self._mul_table = [[self._mul_raw(a, b) for b in range(self.q)] for a in range(self.q)]
        self._neg_table = [self._neg_raw(a) for a in range(self.q)]
        inv_table: list[int | None] = [None] * self.q
        for a in range(1, self.q):
            inv_table[a] = self._pow_raw(a, self.q - 2)
        self._inv_table = inv_table

    def _digits(self, value: int) -> list[int]:
        digits = []
        x = value
        for _ in range(self.degree):
            digits.append(x % self.p)
            x //= self.p
        return digits

    def _from_digits(self, digits: Iterable[int]) -> int:
        value = 0
        factor = 1
        for coeff in digits:
            value += (coeff % self.p) * factor
            factor *= self.p
        return value

    def _add_raw(self, a: int, b: int) -> int:
        if self.degree == 1:
            return (a + b) % self.p
        return self._from_digits((x + y) % self.p for x, y in zip(self._digits(a), self._digits(b)))

    def _neg_raw(self, a: int) -> int:
        if self.degree == 1:
            return (-a) % self.p
        return self._from_digits((-x) % self.p for x in self._digits(a))

    def _mul_raw(self, a: int, b: int) -> int:
        if self.degree == 1:
            return (a * b) % self.p
        if a == 0 or b == 0:
            return 0

        left = self._digits(a)
        right = self._digits(b)
        tmp = [0] * (2 * self.degree - 1)
        for i, ai in enumerate(left):
            for j, bj in enumerate(right):
                tmp[i + j] = (tmp[i + j] + ai * bj) % self.p

        assert self.modulus is not None
        for d in range(len(tmp) - 1, self.degree - 1, -1):
            coeff = tmp[d] % self.p
            if coeff == 0:
                continue
            offset = d - self.degree
            for i in range(self.degree):
                tmp[offset + i] = (tmp[offset + i] - coeff * self.modulus[i]) % self.p

        return self._from_digits(tmp[: self.degree])

    def _pow_raw(self, a: int, exponent: int) -> int:
        result = 1
        power = a
        e = exponent
        while e:
            if e & 1:
                result = self._mul_raw(result, power)
            power = self._mul_raw(power, power)
            e >>= 1
        return result


def GF(q: int, *, modulus: Sequence[int] | None = None, table_limit: int = 2048) -> FiniteField:
    """Construct ``GF(q)``.

    ``q`` must be a prime power. For extension fields, pass ``modulus`` when a
    specific representation is needed; otherwise the first monic irreducible
    polynomial found by deterministic search is used.
    """

    q = int(q)
    if q < 2:
        raise ValueError("field order must be at least 2")
    if _is_prime(q):
        return FiniteField(q, 1, table_limit=table_limit)

    for p in range(2, q + 1):
        if not _is_prime(p):
            continue
        power = p
        degree = 1
        while power < q:
            power *= p
            degree += 1
        if power == q:
            return FiniteField(p, degree, modulus=modulus, table_limit=table_limit)
    raise ValueError("field order must be a prime power")


def product_of(values: Iterable[int]) -> int:
    """Small compatibility helper for callers that need a public product."""

    return prod(values)

