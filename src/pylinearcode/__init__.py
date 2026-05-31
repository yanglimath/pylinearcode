"""Finite-field linear-code algorithms.

The public API intentionally uses plain Python containers. Field elements are
encoded as integers in ``range(q)`` for speed and for easy interoperability with
array backends.
"""

from .bounds import griesmer_bound_length, hamming_sphere_volume, singleton_bound
from .constructors import (
    bch_code,
    bch_generator_polynomial,
    cyclic_code_from_generator,
    extended_hamming_code,
    hamming_code,
    parity_check_code,
    random_linear_code,
    reed_muller_code,
    reed_solomon_code,
    repetition_code,
    simplex_code,
    whole_space_code,
    zero_code,
)
from .enumerators import krawtchouk, macwilliams_transform
from .exceptions import ComputationLimitError, LinearCodeError
from .fields import FiniteField, GF
from .linear_code import LinearCode, SystematicForm

__all__ = [
    "ComputationLimitError",
    "FiniteField",
    "GF",
    "LinearCode",
    "LinearCodeError",
    "SystematicForm",
    "bch_code",
    "bch_generator_polynomial",
    "cyclic_code_from_generator",
    "extended_hamming_code",
    "griesmer_bound_length",
    "hamming_code",
    "hamming_sphere_volume",
    "krawtchouk",
    "macwilliams_transform",
    "parity_check_code",
    "random_linear_code",
    "reed_muller_code",
    "reed_solomon_code",
    "repetition_code",
    "simplex_code",
    "singleton_bound",
    "whole_space_code",
    "zero_code",
]
