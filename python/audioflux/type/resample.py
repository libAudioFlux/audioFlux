from enum import Enum

__all__ = [
    "ResampleAlgType",
    "ResampleQualityType"
]


class ResampleAlgType(Enum):
    """
    Resample Alg Type
    """
    POLYPHASE = 0  # stand
    BANDLIMITED = 1  # ccrma


class ResampleQualityType(Enum):
    """
    Resample Quality Type
    """
    BEST = 0
    MID = 1
    FAST = 2
