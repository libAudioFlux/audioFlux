from enum import Enum

__all__ = [
    'PitchType'
]


class PitchType(Enum):
    """
    Pitch Type
    """
    YIN = 0
    NCF = 1
    PEF = 2
