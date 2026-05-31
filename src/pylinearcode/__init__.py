"""Finite-field linear-code algorithms.

The public API intentionally uses plain Python containers. Field elements are
encoded as integers in ``range(q)`` for speed and for easy interoperability with
array backends.
"""

from .bounds import griesmer_bound_length, hamming_sphere_volume, singleton_bound
from .constructors import (
    cyclic_code_from_generator,
    hamming_code,
    parity_check_code,
    random_linear_code,
    reed_solomon_code,
    repetition_code,
    whole_space_code,
    zero_code,
)
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
    "cyclic_code_from_generator",
    "griesmer_bound_length",
    "hamming_code",
    "hamming_sphere_volume",
    "parity_check_code",
    "random_linear_code",
    "reed_solomon_code",
    "repetition_code",
    "singleton_bound",
    "whole_space_code",
    "zero_code",
]

