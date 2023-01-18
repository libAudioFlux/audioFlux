from enum import Enum

__all__ = [
    "ReduceType",
    "NoveltyType"
]


class ReduceType(Enum):
    """
    Reduce Type
    """
    MEAN = 0  # Spectral flux
    SUM = 1  # Spectral flux

    LOG = 2  # Modified Kullback-Leibler mean
    # MUL = 3 # High Frequency Content mean


class NoveltyType(Enum):
    """
    Novelty Type
    """
    FLUX = 0

    HFC = 1
    SD = 2
    SF = 3
    MKL = 4

    PD = 5
    WPD = 6
    NWPD = 7

    CD = 8
    RCD = 9

    BROADBAND = 10
