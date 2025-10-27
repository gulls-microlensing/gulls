"""Custom exceptions used by the smoke test utilities."""


class SmokeTestError(RuntimeError):
    """Raised when the smoke test encounters a setup or runtime failure."""


__all__ = ["SmokeTestError"]
