"""Helpers for reproducible condition-specific seed derivation."""


def derive_condition_seed(base_seed, *offsets):
    """Derive a reproducible condition-specific seed from one base seed."""
    value = int(base_seed)
    for offset in offsets:
        value = (value * 1009) + int(offset)
    return value
