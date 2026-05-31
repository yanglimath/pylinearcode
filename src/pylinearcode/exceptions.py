"""Project-specific exceptions."""


class LinearCodeError(ValueError):
    """Raised when a linear-code object or operation is invalid."""


class ComputationLimitError(RuntimeError):
    """Raised when an exact algorithm exceeds its configured computation limit."""

