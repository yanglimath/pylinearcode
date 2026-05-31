# PyLinearCode

PyLinearCode is a pure-Python toolkit for finite-field linear codes. It is
designed around a Magma-inspired coding-theory API, while keeping the core small,
auditable, dependency-free, and easy to embed in ordinary Python projects.

The project is not yet a full SageMath or Magma replacement. Its current goal is
to provide a clean open-source kernel for linear-code computations and a stable
place to add faster backends later.

## Highlights

- Finite fields `GF(p)` and `GF(p^m)` with compact integer-encoded elements.
- Public polynomial arithmetic over finite fields.
- Finite-field matrix row reduction, rank, nullspace, and orthogonal complement.
- Linear-code operations:
  - generator and parity-check matrices
  - encoding and membership checks
  - syndrome computation and bounded syndrome decoding
  - syndrome table and covering radius
  - dual code, hull, self-orthogonality, and self-duality
  - puncturing, shortening, and parity extension
  - systematic form
  - exact minimum distance by codeword enumeration or parity-check column search
  - weight distribution and dual weight distribution via MacWilliams
- Code constructors:
  - zero code, whole-space code, repetition code, and single parity-check code
  - Hamming, extended Hamming, and simplex codes
  - generalized Reed-Solomon codes
  - primitive narrow-sense BCH codes over prime fields
  - binary Reed-Muller codes
  - random linear codes
  - cyclic codes from generator polynomials
- Classical bounds:
  - Singleton
  - Hamming sphere volume
  - Griesmer length bound

## Installation

From the repository root:

```bash
python -m pip install -e .
```

For development tools:

```bash
python -m pip install -e ".[dev]"
```

## Quick Start

```python
from pylinearcode import GF, hamming_code, reed_solomon_code

F2 = GF(2)
C = hamming_code(F2, r=3)

print(C.parameters())          # (7, 4, 3)
print(C.minimum_distance())    # 3

word = C.encode([1, 0, 1, 1])
received = word[:]
received[2] ^= 1

decoded, error = C.syndrome_decode(received, max_weight=1)
assert decoded == word

F5 = GF(5)
RS = reed_solomon_code(F5, length=5, dimension=3)
print(RS.parameters())         # (5, 3, 3)
```

## BCH Codes

```python
from pylinearcode import GF, bch_code, bch_generator_polynomial

F2 = GF(2)
g = bch_generator_polynomial(F2, extension_degree=4, designed_distance=5)
print(g)                       # coefficients in ascending order

C = bch_code(F2, extension_degree=4, designed_distance=5)
print(C.length, C.dimension)   # 15 7
print(C.minimum_distance())    # 5
```

The BCH constructor currently supports primitive narrow-sense BCH codes over
prime fields `GF(p)` with length `p^m - 1`.

## Reed-Muller Codes

```python
from pylinearcode import reed_muller_code

RM = reed_muller_code(order=1, variables=3)
print(RM.parameters())          # (8, 4, 4)
print(RM.weight_distribution()) # A_4 = 14, A_8 = 1
```

## Polynomial Arithmetic

```python
from pylinearcode import GF
from pylinearcode.polynomials import divmod, evaluate, gcd, mul

F2 = GF(2)
f = mul(F2, [1, 1], [1, 0, 1])
q, r = divmod(F2, f, [1, 1])

assert f == [1, 1, 1, 1]
assert q == [1, 0, 1]
assert r == [0]
assert gcd(F2, f, [1, 1]) == [1, 1]
assert evaluate(F2, [1, 1, 1], 1) == 1
```

## Command Line

```bash
python -m pylinearcode hamming --q 2 --r 3
python -m pylinearcode rs --q 5 --n 5 --k 3
python -m pylinearcode bch --q 2 --m 4 --delta 5
python -m pylinearcode rm --r 1 --m 3
```

## Verification

```bash
PYTHONPATH=src python -m unittest discover -s tests
```

## Roadmap

Near-term algorithmic targets:

1. Information-set decoding and Lee-Brickell-style improvements.
2. Brouwer-Zimmermann and Leon minimum-distance algorithms.
3. BCH variants, Goppa codes, alternant codes, and AG-code scaffolding.
4. Code equivalence and automorphism-group routines.
5. Optional Numba, Cython, Rust, or GPU backends for hot integer-matrix paths.

The design principle is simple: keep the public Python API stable and
mathematically explicit, then accelerate the internals without changing user
code.
