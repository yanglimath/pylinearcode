"""Polynomial arithmetic over finite fields.

Polynomials use ascending coefficient order:

``[a0, a1, a2]`` means ``a0 + a1*x + a2*x^2``.
"""

from __future__ import annotations

from typing import Sequence

from .fields import FiniteField

Polynomial = list[int]


def trim(poly: Sequence[int]) -> Polynomial:
    result = [int(c) for c in poly]
    while result and result[-1] == 0:
        result.pop()
    return result or [0]


def degree(poly: Sequence[int]) -> int:
    poly = trim(poly)
    return -1 if poly == [0] else len(poly) - 1


def normalize(field: FiniteField, poly: Sequence[int]) -> Polynomial:
    return trim([field.coerce(c) for c in poly])


def add(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Polynomial:
    n = max(len(left), len(right))
    out = [0] * n
    for i in range(n):
        a = left[i] if i < len(left) else 0
        b = right[i] if i < len(right) else 0
        out[i] = field.add(a, b)
    return trim(out)


def sub(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Polynomial:
    n = max(len(left), len(right))
    out = [0] * n
    for i in range(n):
        a = left[i] if i < len(left) else 0
        b = right[i] if i < len(right) else 0
        out[i] = field.sub(a, b)
    return trim(out)


def neg(field: FiniteField, poly: Sequence[int]) -> Polynomial:
    return trim([field.neg(c) for c in poly])


def scalar_mul(field: FiniteField, scalar: int, poly: Sequence[int]) -> Polynomial:
    scalar = field.coerce(scalar)
    if scalar == 0:
        return [0]
    return trim([field.mul(scalar, c) for c in poly])


def mul(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Polynomial:
    left = normalize(field, left)
    right = normalize(field, right)
    if left == [0] or right == [0]:
        return [0]
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        if a == 0:
            continue
        for j, b in enumerate(right):
            if b != 0:
                out[i + j] = field.add(out[i + j], field.mul(a, b))
    return trim(out)


def divmod(field: FiniteField, dividend: Sequence[int], divisor: Sequence[int]) -> tuple[Polynomial, Polynomial]:
    dividend = normalize(field, dividend)
    divisor = normalize(field, divisor)
    if divisor == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(dividend) < len(divisor):
        return [0], dividend

    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = dividend[:]
    inv_lead = field.inv(divisor[-1])
    while len(remainder) >= len(divisor) and remainder != [0]:
        shift = len(remainder) - len(divisor)
        coeff = field.mul(remainder[-1], inv_lead)
        quotient[shift] = coeff
        for i, di in enumerate(divisor):
            remainder[shift + i] = field.sub(remainder[shift + i], field.mul(coeff, di))
        remainder = trim(remainder)
    return trim(quotient), remainder


def mod(field: FiniteField, dividend: Sequence[int], modulus: Sequence[int]) -> Polynomial:
    return divmod(field, dividend, modulus)[1]


def monic(field: FiniteField, poly: Sequence[int]) -> Polynomial:
    poly = normalize(field, poly)
    if poly == [0]:
        return [0]
    inv_lead = field.inv(poly[-1])
    return scalar_mul(field, inv_lead, poly)


def gcd(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Polynomial:
    a = normalize(field, left)
    b = normalize(field, right)
    while b != [0]:
        _, rem = divmod(field, a, b)
        a, b = b, rem
    return monic(field, a)


def lcm(field: FiniteField, left: Sequence[int], right: Sequence[int]) -> Polynomial:
    left = normalize(field, left)
    right = normalize(field, right)
    if left == [0] or right == [0]:
        return [0]
    common = gcd(field, left, right)
    quotient, remainder = divmod(field, left, common)
    if remainder != [0]:
        raise ArithmeticError("internal polynomial division error")
    return monic(field, mul(field, quotient, right))


def powmod(field: FiniteField, base: Sequence[int], exponent: int, modulus: Sequence[int]) -> Polynomial:
    if exponent < 0:
        raise ValueError("negative polynomial exponents are not supported")
    result = [field.one]
    power = mod(field, base, modulus)
    e = exponent
    while e:
        if e & 1:
            result = mod(field, mul(field, result, power), modulus)
        power = mod(field, mul(field, power, power), modulus)
        e >>= 1
    return result


def derivative(field: FiniteField, poly: Sequence[int]) -> Polynomial:
    poly = normalize(field, poly)
    if len(poly) <= 1:
        return [0]
    return trim([field.mul(i % field.p, poly[i]) for i in range(1, len(poly))])


def evaluate(field: FiniteField, poly: Sequence[int], point: int) -> int:
    point = field.coerce(point)
    total = field.zero
    for coeff in reversed(normalize(field, poly)):
        total = field.add(field.mul(total, point), coeff)
    return total


def format_polynomial(field: FiniteField, poly: Sequence[int], variable: str = "x") -> str:
    poly = normalize(field, poly)
    if poly == [0]:
        return "0"
    terms: list[str] = []
    for power, coeff in enumerate(poly):
        if coeff == 0:
            continue
        coeff_text = field.format(coeff)
        if power == 0:
            terms.append(coeff_text)
        elif power == 1:
            terms.append(variable if coeff == 1 else f"{coeff_text}*{variable}")
        else:
            terms.append(f"{variable}^{power}" if coeff == 1 else f"{coeff_text}*{variable}^{power}")
    return " + ".join(terms)

