TRUE_VALUES = {"1", "true", "yes", "y"}


def parse_bool_arg(value: str) -> bool:
    """Parse common CLI truthy values."""
    return value.strip().lower() in TRUE_VALUES
